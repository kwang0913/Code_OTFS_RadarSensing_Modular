function final_ests = estimate_ssr_params(YTF, XTF, coarse, params, chanParams, n0, derivedParams)
%ESTIMATE_SSR_PARAMS Direct OMP on fine dictionary for monostatic DFRC.
%   Stage 1: fine (angle, delay, Doppler) dictionary from coarse candidate
%   neighborhoods (angle step ssr.length.a, delay/Doppler steps at the min
%   steps), OMP with residual-reduction threshold stopping.
%   Stage 2: angle of each selected path refined within +-ssr.length.a at the
%   angle min step (delay/Doppler held), then gains refit.
%   Gains use the pilot model y = c*Phi*gains + noise (pilot_noise_model).
%   params.ssr.max_J optionally caps the number of recovered scatterers.

    K  = length(coarse.AoAs);
    pilot_blocks = params.pilot_blocks;
    Nr = params.Nr;
    M  = params.M;
    N  = params.N;

    max_targets = params.targets_num;
    if ~isempty(params.ssr.max_J)
        if ~isscalar(params.ssr.max_J) || ~isreal(params.ssr.max_J) || ...
                ~isfinite(params.ssr.max_J) || params.ssr.max_J < 1 || ...
                params.ssr.max_J ~= floor(params.ssr.max_J)
            error('ssr.max_J must be a positive integer.');
        end
        max_targets = params.ssr.max_J;
    end

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr, n0);
    [c, s2] = pilot_noise_model(XTF, params, n0);  % y = c*Phi*gains + noise(s2)

    % --- Fine grid config ---
    min_step = ssr_min_steps(params);
    step_a = params.ssr.length.a;  % stage-1 angle step (degrees)
    step_r = min_step.r;           % delay step (taps)
    step_v = min_step.v;           % Doppler step (taps)
    steps = [step_a, step_r, step_v];
    if numel(steps) ~= 3 || ~isreal(steps) || any(~isfinite(steps) | steps <= 0)
        error('SSR grid steps must be positive finite scalars.');
    end
    range_a = params.ssr.range.a;
    range_r = params.ssr.range.r;
    range_v = params.ssr.range.v;

    %% --- Build fine dictionary from all coarse candidates ---
    % Collect all grid points, then deduplicate
    all_a = zeros(0, 1);  all_r = zeros(0, 1);  all_v = zeros(0, 1);
    for ci = 1:K
        if params.skip_aoa
            grid_a = 0;  % fixed angle, no search
        else
            grid_a = coarse.AoAs(ci) - range_a : step_a : coarse.AoAs(ci) + range_a;
        end
        delay_bounds = round((coarse.pathDelays(ci) + [-range_r, range_r]) / step_r, 12);
        grid_r = (ceil(delay_bounds(1)):floor(delay_bounds(2))) * step_r;
        grid_r = grid_r(grid_r >= 0 & grid_r < M);
        doppler_bounds = round((coarse.pathDopplers(ci) + [-range_v, range_v]) / step_v, 12);
        grid_v = (ceil(doppler_bounds(1)):floor(doppler_bounds(2))) * step_v;

        [A, R, V] = ndgrid(grid_a, grid_r, grid_v);
        all_a = [all_a; A(:)]; %#ok<AGROW>
        all_r = [all_r; R(:)]; %#ok<AGROW>
        all_v = [all_v; V(:)]; %#ok<AGROW>
    end

    % Remove floating-point duplicates while preserving fractional taps.
    arv = round([all_a, all_r, all_v], 12);
    arv_unique = unique(arv, 'rows', 'stable');
    dict_a = arv_unique(:, 1)';
    dict_r = arv_unique(:, 2)';
    dict_v = arv_unique(:, 3)';
    n_dict = length(dict_a);

    if isfield(params, 'verbose_output') && params.verbose_output
        fprintf('  SSR dictionary: %d atoms from %d coarse candidates\n', n_dict, K);
    end

    % Build dictionary matrix: vectorize over (delay, Doppler) per angle.
    % Each angle is built in a small local block and written with ONE full-height assignment:
    % per-site row-block writes into scattered columns of the large matrix were pathologically slow.
    Phi_dict = zeros(n_meas, n_dict);
    sites = pilot_sites(pilot_blocks);
    unique_angles = unique(dict_a);

    for ai = 1:length(unique_angles)
        theta = unique_angles(ai);
        col_mask = (dict_a == theta);
        col_idx = find(col_mask);
        local_r = dict_r(col_mask);
        local_v = dict_v(col_mask);

        a_steer = exp(1j*2*pi*params.ddlr*(0:Nr-1).'*sind(theta));

        nu_all = local_v(:);

        Phi_blk = zeros(n_meas, numel(col_idx));
        for s = 1:sites.n_sites
            d_tx = exp(-1j*2*pi*params.ddlt*(sites.tx(s)-1)*sind(theta));

            % xi_H for each atom (delay and Doppler paired elementwise)
            xi_vals = compute_xi_H(nu_all, local_r(:), sites.n(s)-1, sites.m(s)-1, N, M, derivedParams, params.padType);
            Phi_blk((s-1)*Nr+(1:Nr), :) = (a_steer * d_tx) * xi_vals.';
        end
        Phi_dict(:, col_idx) = c * Phi_blk;
    end

    %% --- OMP with residual-reduction threshold ---
    max_targets = min(max_targets, n_dict);  % cannot select more distinct atoms

    r = y;
    active = [];
    prev_res = real(y' * y);
    cn = vecnorm(Phi_dict);
    cn(cn == 0) = Inf;

    for omp_iter = 1:max_targets
        % Correlate residual with all dictionary columns
        corr_vals = abs(Phi_dict' * r) ./ cn';

        % Zero out already-selected columns
        corr_vals(active) = 0;

        [~, best] = max(corr_vals);
        active = [active, best]; %#ok<AGROW>

        % Joint LS for active set
        Phi_active = Phi_dict(:, active);
        gains_active = (Phi_active' * Phi_active + s2 * eye(length(active))) \ (Phi_active' * y);
        r = y - Phi_active * gains_active;

        cur_res = real(r' * r);
        reduction = (prev_res - cur_res) / prev_res;

        if isfield(params, 'verbose_output') && params.verbose_output
            fprintf('  OMP iter %d: atom %d (a=%.1f r=%.1f v=%.1f), ||r||²=%.2f, reduction=%.1f%%\n', ...
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

    %% --- Stage 2: refine each selected angle within +-step_a at min_step.a ---
    % Delay/Doppler held from stage 1; one pass, each path against the others' fit.
    sel_a = dict_a(active);  sel_r = dict_r(active);  sel_v = dict_v(active);
    n_ref = round(step_a / min_step.a);
    if ~params.skip_aoa && min_step.a < step_a  % stage 2 only if strictly finer than stage 1
        n_sel = length(active);
        for i = 1:n_sel
            others = [1:i-1, i+1:n_sel];
            Phi_sel = atoms(sel_a, sel_r, sel_v, c, params, derivedParams);
            g = (Phi_sel' * Phi_sel + s2 * eye(n_sel)) \ (Phi_sel' * y);
            r_i = y - Phi_sel(:, others) * reshape(g(others), [], 1);  % reshape: n_sel = 1 gives a 0x1 empty
            cand_a = round(sel_a(i) + (-n_ref:n_ref) * min_step.a, 12);
            Phi_c = atoms(cand_a, repmat(sel_r(i), size(cand_a)), repmat(sel_v(i), size(cand_a)), c, params, derivedParams);
            [~, best_c] = max(abs(Phi_c' * r_i) ./ vecnorm(Phi_c)');
            sel_a(i) = cand_a(best_c);
        end
    end

    %% --- Package output ---
    Phi_active = atoms(sel_a, sel_r, sel_v, c, params, derivedParams);
    gains_active = (Phi_active' * Phi_active + s2 * eye(length(active))) \ (Phi_active' * y);
    r = y - Phi_active * gains_active;

    final_ests.AoAs         = sel_a;
    final_ests.AoDs         = sel_a;  % monostatic
    final_ests.pathDelays   = sel_r;
    final_ests.pathDopplers = sel_v;
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
            for j = 1:length(chanParams.(fields{fi}))
                xline(chanParams.(fields{fi})(j), 'r--', 'LineWidth', 2, ...
                    'DisplayName', sprintf('True %.1f', chanParams.(fields{fi})(j)));
            end
            hold off;
            xlabel(labels{fi}); ylabel('|Gain|');
            title(sprintf('SSR: %s', labels{fi}));
            set(gca, 'FontSize', 11, 'FontWeight', 'Bold');
            legend('Location', 'best');
        end
    end
end

function Phi = atoms(angles, delays, dopplers, c, params, derivedParams)
%ATOMS Scaled response columns c*Phi for the given (angle, delay, Doppler) triples (monostatic).
    ests = struct('AoAs', angles, 'AoDs', angles, 'pathDelays', delays, 'pathDopplers', dopplers);
    Phi = c * build_response_matrix(ests, params, derivedParams);
end
