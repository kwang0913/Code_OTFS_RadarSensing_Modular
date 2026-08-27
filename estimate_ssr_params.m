function final_ests = estimate_ssr_params(YTF, XTF, coarse, params, n0)
%ESTIMATE_SSR_PARAMS Direct OMP on fine dictionary for monostatic DFRC.
%   Builds a fine (angle, delay, Doppler) dictionary from coarse candidate
%   neighborhoods, then runs OMP with residual-reduction threshold stopping.
%   No binary refinement — single-stage, grid-independent.

    K  = length(coarse.AoAs);
    pilot_blocks = params.pilot_blocks;
    Nr = params.Nr;
    M  = params.M;
    N  = params.N;

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr);

    % --- Fine grid config ---
    step_a = 1;     % angle step (degrees)
    step_r = 1;     % delay step (integer taps)
    step_v = 0.1;   % Doppler step (taps)
    range_a = params.ssr.range.a;
    range_r = params.ssr.range.r;
    range_v = params.ssr.range.v;

    %% --- Build fine dictionary from all coarse candidates ---
    % Collect all grid points, then deduplicate
    all_a = [];  all_r = [];  all_v = [];
    for ci = 1:K
        if params.skip_aoa
            grid_a = 0;  % fixed angle, no search
        else
            grid_a = coarse.AoAs(ci) - range_a : step_a : coarse.AoAs(ci) + range_a;
        end
        grid_r = round(coarse.pathDelays(ci) - range_r) : step_r : round(coarse.pathDelays(ci) + range_r);
        grid_r = grid_r(grid_r >= 0);  % delay >= 0
        grid_v = coarse.pathDopplers(ci) - range_v : step_v : coarse.pathDopplers(ci) + range_v;

        [A, R, V] = ndgrid(grid_a, grid_r, grid_v);
        all_a = [all_a; A(:)]; %#ok<AGROW>
        all_r = [all_r; R(:)]; %#ok<AGROW>
        all_v = [all_v; V(:)]; %#ok<AGROW>
    end

    % Deduplicate: round to grid precision (angle 0.1°, delay integer, Doppler 0.01)
    arv = [round(all_a, 1), round(all_r), round(all_v, 2)];
    [arv_unique, ~] = unique(arv, 'rows', 'stable');
    dict_a = arv_unique(:, 1)';
    dict_r = arv_unique(:, 2)';
    dict_v = arv_unique(:, 3)';
    n_dict = length(dict_a);

    if isfield(params, 'verbose_output') && params.verbose_output
        fprintf('  SSR dictionary: %d atoms from %d coarse candidates\n', n_dict, K);
    end

    % Build dictionary matrix: vectorize over (delay, Doppler) per angle
    Phi_dict = zeros(n_meas, n_dict);
    unique_angles = unique(dict_a);

    for ai = 1:length(unique_angles)
        theta = unique_angles(ai);
        col_mask = (dict_a == theta);
        col_idx = find(col_mask);
        local_r = dict_r(col_mask);
        local_v = dict_v(col_mask);

        a_steer = exp(1j*2*pi*params.ddlr*(0:Nr-1).'*sind(theta));

        % Group by unique delay values for vectorized compute_xi_H
        unique_r = unique(local_r);
        nu_all = local_v(:);
        k_all  = floor(nu_all);
        kappa_all = nu_all - k_all;

        idx = 0;
        for b = 1:length(pilot_blocks)
            blk = pilot_blocks(b);
            d_tx = exp(-1j*2*pi*params.ddlt*(blk.tx-1)*sind(theta));

            for pp = 1:blk.n_pilots
                m0 = blk.m(pp) - 1;
                n0_pos = blk.n(pp) - 1;

                % Compute xi_H for each atom individually (paired mode)
                xi_vals = compute_xi_H(nu_all, k_all, kappa_all, local_r(:), n0_pos, m0, N, M);
                Phi_dict(idx+(1:Nr), col_idx) = (a_steer * d_tx) * xi_vals.';
                idx = idx + Nr;
            end
        end
    end

    %% --- OMP with residual-reduction threshold ---
    max_targets = params.targets_num * 2;  % upper bound: let OMP find more than coarse

    r = y;
    active = [];
    prev_res = real(y' * y);

    for omp_iter = 1:max_targets
        % Correlate residual with all dictionary columns
        cn = vecnorm(Phi_dict);
        cn(cn == 0) = Inf;
        corr_vals = abs(Phi_dict' * r) ./ cn';

        % Zero out already-selected columns
        for si = 1:length(active)
            corr_vals(active(si)) = 0;
        end

        [~, best] = max(corr_vals);
        active = [active, best]; %#ok<AGROW>

        % Joint LS for active set
        Phi_active = Phi_dict(:, active);
        gains_active = (Phi_active' * Phi_active + n0 * eye(length(active))) \ (Phi_active' * y);
        r = y - Phi_active * gains_active;

        cur_res = real(r' * r);
        reduction = (prev_res - cur_res) / prev_res;

        if isfield(params, 'verbose_output') && params.verbose_output
            fprintf('  OMP iter %d: atom %d (a=%.0f r=%d v=%.1f), ||r||²=%.2f, reduction=%.1f%%\n', ...
                omp_iter, best, dict_a(best), dict_r(best), dict_v(best), ...
                cur_res, reduction*100);
        end

        % Stop when adding an atom barely reduces residual (keep at least 1)
        if reduction < 0.1 && length(active) > 1
            active(end) = [];
            break;
        end

        prev_res = cur_res;
    end

    %% --- Package output ---
    Phi_active = Phi_dict(:, active);
    gains_active = (Phi_active' * Phi_active + n0 * eye(length(active))) \ (Phi_active' * y);
    r = y - Phi_active * gains_active;

    final_ests.AoAs         = dict_a(active);
    final_ests.pathDelays   = dict_r(active);
    final_ests.pathDopplers = dict_v(active);
    final_ests.pathGains    = gains_active;
    final_ests.residual_nmse = real(r' * r) / real(y' * y);

    if isfield(params, 'verbose_output') && params.verbose_output
        fprintf('--- SSR Results ---\n');
        if params.skip_aoa
            fprintf('AoA:     unresolvable (Nr=1)\n');
        else
            fprintf('AoA:     %s\n', sprintf('%5.2f ', final_ests.AoAs));
        end
        fprintf('Delay:   %s\n', sprintf('%5.2f ', final_ests.pathDelays));
        fprintf('Doppler: %s\n', sprintf('%5.2f ', final_ests.pathDopplers));
        fprintf('--------------------------------------------\n');
    end

    % Optional: SSR result stem plots
    if isfield(params, 'show_figures') && params.show_figures
        gain_mag = abs(final_ests.pathGains(:)');
        labels = {'Angle [°]', 'Delay [tap]', 'Doppler [tap]'};
        fields = {'AoAs', 'pathDelays', 'pathDopplers'};

        figure;
        for fi = 1:3
            subplot(1, 3, fi); hold on;
            stem(final_ests.(fields{fi}), gain_mag, 'b', 'filled', ...
                'LineWidth', 1.5, 'DisplayName', 'SSR Est.');
            if isfield(params, 'chanParams_for_plot')
                cp = params.chanParams_for_plot;
                for j = 1:length(cp.(fields{fi}))
                    xline(cp.(fields{fi})(j), 'r--', 'LineWidth', 2, ...
                        'DisplayName', sprintf('True %.1f', cp.(fields{fi})(j)));
                end
            end
            hold off;
            xlabel(labels{fi}); ylabel('|Gain|');
            title(sprintf('SSR: %s', labels{fi}));
            set(gca, 'FontSize', 11, 'FontWeight', 'Bold');
            legend('Location', 'best');
        end
    end
end
