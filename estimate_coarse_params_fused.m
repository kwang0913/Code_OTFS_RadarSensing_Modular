function coarse = estimate_coarse_params_fused(YTF, YDD, XTF, XDD, params, chanParams, n0, derivedParams)
%ESTIMATE_COARSE_PARAMS_FUSED Coarse DFRC parameter estimation using pilot arms.
%   Stage 1:  SIC AoA detection
%   Stage 2:  FFT xcorr detection + pairing
%   Stage 2b: pilot-arm CLEAN: xcorr candidates snapped on the residual arm
%             spectra, refined on a local coherent LS grid, selected by LS
%             improvement, joint gain refit
%   Stage 3:  ridge LS path gains, model y = c*Phi*gains + noise (pilot_noise_model)

    M  = params.M;
    N  = params.N;
    Nr = params.Nr;
    Nt = params.Nt;
    max_targets = params.targets_num;
    interpolate_len = 64;  % FFT zero-pad length for angular spectrum interpolation

    if params.skip_aoa
        %% --- Nr=1: skip SIC, fix AoA=0 ---
        detected_angles = 0;
        n_detected = 1;
        ang = 0;
    else
        %% -------------------------------------------------------------------
        %  Stage 1: SIC AoA Detection
        % --------------------------------------------------------------------
        rx_OTFS_resid = reshape(YTF, [], Nr);
        x_ang = asind(((0:interpolate_len-1) - interpolate_len/2) / (params.ddlr * interpolate_len));

        detected_angles = [];
        original_peak   = 0;
        Cor_angle_init  = [];
        sic_spectra   = {};
        sic_stop_reason = '';

        for iter = 1:max_targets
            rx_idft = fft(rx_OTFS_resid, interpolate_len, 2);
            Cor_angle = fftshift(median(abs(rx_idft), 1));

            if iter == 1
                original_peak  = max(Cor_angle);
                Cor_angle_init = Cor_angle / original_peak;
            end

            sic_spectra{iter} = Cor_angle / original_peak; %#ok<AGROW>

            [pks, locs, ~, proms] = findpeaks(Cor_angle, 'SortStr', 'descend');
            if isempty(pks)
                sic_stop_reason = 'no peaks found';
                break;
            end
            if proms(1) / original_peak < 0.1 && iter > 1
                sic_stop_reason = sprintf('prominence %.3f < 0.1', proms(1)/original_peak);
                break;
            end

            spec      = Cor_angle / original_peak;
            spec_mean = mean(spec);
            spec_med  = median(spec);
            if spec_mean < spec_med && ...
               (pks(1)/original_peak - spec_med) < (spec_med - spec_mean) && ...
               iter > 1
                sic_stop_reason = sprintf('mean/med: peak-med=%.3f < med-mean=%.3f', ...
                    pks(1)/original_peak - spec_med, spec_med - spec_mean);
                break;
            end

            theta_k = x_ang(locs(1));
            detected_angles(end+1) = theta_k; %#ok<AGROW>

            a_k     = exp(1j * 2*pi * (0:Nr-1)' * params.ddlr * sind(theta_k));
            alpha_k = (a_k' * a_k) \ (a_k' * rx_OTFS_resid.');
            rx_OTFS_resid = rx_OTFS_resid - (a_k * alpha_k).';
        end

        if isempty(sic_stop_reason) && length(detected_angles) == max_targets
            sic_stop_reason = 'reached max_targets';
        end

        n_detected = length(detected_angles);
        if n_detected == 0
            rx_idft = fft(reshape(YTF, [], Nr), interpolate_len, 2);
            Cor_angle = fftshift(median(abs(rx_idft), 1));
            [~, fb_idx] = max(Cor_angle);
            detected_angles = x_ang(fb_idx);
            n_detected = 1;
            Cor_angle_init = Cor_angle / max(Cor_angle);
        end

        ang = round(detected_angles);

        % Plot SIC step-by-step debug figure
        if params.show_figures
            n_iters = length(sic_spectra);
            figure;
            for it = 1:n_iters
                subplot(n_iters, 1, it); hold on;
                plot(x_ang, sic_spectra{it}, '-b', 'LineWidth', 1.5);
                for j = 1:min(max_targets, length(chanParams.AoAs))
                    xline(chanParams.AoAs(j), 'r--', 'LineWidth', 1);
                end
                if it <= n_detected
                    xline(detected_angles(it), 'g-', 'LineWidth', 2);
                    title_str = sprintf('Iter %d: detected %.1f°', it, detected_angles(it));
                else
                    title_str = sprintf('Iter %d: STOPPED — %s', it, sic_stop_reason);
                end
                yline(0.1, 'k:', 'LineWidth', 1);
                hold off;
                ylabel('Norm. Mag.');
                title(title_str, 'FontSize', 10);
                set(gca, 'FontSize', 9, 'FontWeight', 'Bold');
                if it == n_iters, xlabel('Angle [°]'); end
            end
            sgtitle('SIC AoA Detection (step-by-step)');

            figure; hold on;
            plot(x_ang, Cor_angle_init, '-', 'DisplayName', 'SIC Spectrum', 'LineWidth', 2);
            for j = 1:min(max_targets, length(chanParams.AoAs))
                xline(chanParams.AoAs(j), 'r--', 'LineWidth', 2, ...
                    'DisplayName', sprintf('True %.1f°', chanParams.AoAs(j)));
            end
            for j = 1:n_detected
                xline(ang(j), 'g:', 'LineWidth', 2, ...
                    'DisplayName', sprintf('Est %.0f°', ang(j)));
            end
            hold off;
            legend('Location', 'southwest');
            xlabel('Angle [\circ]'); ylabel('Normalized Magnitude');
            set(gca, 'FontSize', 11, 'FontWeight', 'Bold');
            title('Coarse: SIC Angular Spectrum');
        end
    end

    %% -------------------------------------------------------------------
    %  Stage 2: FFT circular xcorr for delay-Doppler (with beamforming when Nr > 1)
    % --------------------------------------------------------------------
    xcorr2_thresh = 0.1;  % fraction of per-bin peak max for local maxima detection
    doppler_taps  = 0:N-1;  % circular Doppler lags; wrapN maps to (-N/2, N/2]
    wrapN = @(v) N/2 - mod(N/2 - v, N);  % principal interval (-N/2, N/2]
    n_ang = length(ang);
    XDD_flat = reshape(permute(XDD, [3,1,2]), Nt, []);  % Nt x (M*N)

    if params.skip_aoa
        % Nr=1: no beamforming, single DD grid
        ADD = reshape(squeeze(YDD), 1, M, N);
        XDD_ref = reshape(sum(XDD, 3), 1, M, N);
    else
        % MMSE beamforming weights
        Ang_mtx = exp(1j * 2*pi * (0:Nr-1)' * params.ddlr * sind(ang));
        W_bf = (Ang_mtx' * Ang_mtx + n0 * eye(n_ang)) \ Ang_mtx';

        % Per-angle DD grids via MMSE beamforming (vectorized)
        ADD = reshape(W_bf * reshape(permute(YDD, [3,1,2]), Nr, []), n_ang, M, N);

        % Tx-steered DD reference per angle (vectorized)
        XDD_ref = zeros(n_ang, M, N);
        for p = 1:n_ang
            d_tx = exp(-1j*2*pi*(0:Nt-1)'*params.ddlt*sind(ang(p)));
            XDD_ref(p,:,:) = reshape(d_tx.' * XDD_flat, M, N);
        end
    end

    % Circular xcorr per angle with local maxima detection
    all_angles    = [];
    all_delays    = [];
    all_dopplers  = [];
    all_peak_vals = [];

    for p = 1:n_ang
        Y_bf  = squeeze(ADD(p, :, :));
        X_ref = squeeze(XDD_ref(p, :, :));

        Hdd_pos = abs(ifft2(fft2(Y_bf) .* conj(fft2(X_ref))));  % M x N circular xcorr: row r -> delay r-1, col c -> Doppler lag c-1

        peak_max = max(Hdd_pos(:));
        if peak_max == 0, continue; end

        thresh = xcorr2_thresh * peak_max;

        % Local maxima: >= all 8 circular neighbors and above threshold
        is_peak = true(size(Hdd_pos));
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0, continue; end
                is_peak = is_peak & (Hdd_pos >= circshift(Hdd_pos, [dr dc]));
            end
        end
        is_peak = is_peak & (Hdd_pos >= thresh);

        [r_idx, c_idx] = find(is_peak);
        for k = 1:length(r_idx)
            all_angles(end+1)    = ang(p); %#ok<AGROW>
            all_delays(end+1)    = r_idx(k) - 1; %#ok<AGROW>
            all_dopplers(end+1)  = wrapN(doppler_taps(c_idx(k))); %#ok<AGROW>
            all_peak_vals(end+1) = Hdd_pos(r_idx(k), c_idx(k)); %#ok<AGROW>
        end
    end

    % Fallback: if no peaks found, use global max of first angle bin
    if isempty(all_delays) && n_ang > 0
        Y_bf  = squeeze(ADD(1, :, :));
        X_ref = squeeze(XDD_ref(1, :, :));
        Hdd_pos = abs(ifft2(fft2(Y_bf) .* conj(fft2(X_ref))));
        [~, pk_lin] = max(Hdd_pos(:));
        [row_pk, col_pk] = ind2sub(size(Hdd_pos), pk_lin);
        all_angles    = ang(1);
        all_delays    = row_pk - 1;
        all_dopplers  = wrapN(doppler_taps(col_pk));
        all_peak_vals = Hdd_pos(row_pk, col_pk);
    end

    % --- Deduplicate: for shared (delay, Doppler) peaks at close angles, keep strongest ---
    % Only treat as leakage if angles are within ~1 beamwidth
    beamwidth_deg = asind(1 / (Nr * params.ddlr));  % ULA 3dB beamwidth ≈ asin(1/(Nr·d/λ))
    ang_dedup_thresh = beamwidth_deg;

    n_peaks = length(all_delays);
    keep = true(1, n_peaks);
    for i = 1:n_peaks
        if ~keep(i), continue; end
        for j = i+1:n_peaks
            if ~keep(j), continue; end
            if all_delays(i) == all_delays(j) && all_dopplers(i) == all_dopplers(j) ...
                    && abs(all_angles(i) - all_angles(j)) < ang_dedup_thresh
                % Same (delay, Doppler) at close angles → leakage
                if all_peak_vals(i) >= all_peak_vals(j)
                    keep(j) = false;
                else
                    keep(i) = false;
                    break;
                end
            end
        end
    end

    if params.verbose_output && sum(~keep) > 0
        fprintf('  Dedup: removed %d leaked peak(s)\n', sum(~keep));
    end

    all_angles    = all_angles(keep);
    all_delays    = all_delays(keep);
    all_dopplers  = all_dopplers(keep);
    all_peak_vals = all_peak_vals(keep);

    n_detected   = length(all_delays);
    ang          = all_angles;
    pathDelays   = all_delays;
    pathDopplers = all_dopplers;

    % Plot xcorr2 surfaces
    if params.show_figures
        unique_ang = unique(ang, 'stable');
        n_ang_bins = length(unique_ang);
        for ab = 1:min(n_ang_bins, 4)
            ang_val = unique_ang(ab);
            tgt_mask = (ang == ang_val);

            [~, bf_idx] = min(abs(detected_angles - ang_val));
            Y_bf_p  = squeeze(ADD(bf_idx, :, :));
            X_ref_p = squeeze(XDD_ref(bf_idx, :, :));
            Hdd_pos_plot = abs(ifft2(fft2(Y_bf_p) .* conj(fft2(X_ref_p))));
            Hdd_pos_plot = Hdd_pos_plot / max(Hdd_pos_plot(:));

            delay_axis = 0:(M-1);

            figure;
            mesh(-N/2:N/2-1, delay_axis, fftshift(Hdd_pos_plot, 2), 'EdgeColor', 'interp', 'FaceColor', 'interp');
            hold on;
            for tidx = 1:length(chanParams.AoAs)
                [~, nearest_bin] = min(abs(unique_ang - chanParams.AoAs(tidx)));
                if nearest_bin == ab
                    stem3(chanParams.pathDopplers(tidx), chanParams.pathDelays(tidx), 1, ...
                        'go', 'filled', 'MarkerSize', 10, 'LineWidth', 2, ...
                        'DisplayName', sprintf('True tgt %d (%.0f°)', tidx, chanParams.AoAs(tidx)));
                end
            end
            % Find all local maxima above threshold on this surface
            is_pk = true(size(Hdd_pos_plot));
            for dr = -1:1
                for dc = -1:1
                    if dr == 0 && dc == 0, continue; end
                    is_pk = is_pk & (Hdd_pos_plot >= circshift(Hdd_pos_plot, [dr dc]));
                end
            end
            is_pk = is_pk & (Hdd_pos_plot >= xcorr2_thresh);
            [pk_r, pk_c] = find(is_pk);

            % Plot all peaks above threshold
            est_set_delays   = pathDelays(tgt_mask);
            est_set_dopplers = pathDopplers(tgt_mask);
            for qi = 1:length(pk_r)
                pk_delay = pk_r(qi) - 1;
                pk_dopp  = wrapN(doppler_taps(pk_c(qi)));
                pk_mag   = Hdd_pos_plot(pk_r(qi), pk_c(qi));

                % Check if this peak is a selected estimate
                is_est = any(est_set_delays == pk_delay & est_set_dopplers == pk_dopp);
                if is_est
                    marker = 'r^';  face = 'r';  clr = 'r';
                    lbl = sprintf('Est (l=%d, v=%.0f) mag=%.2f', pk_delay, pk_dopp, pk_mag);
                else
                    marker = 'bs';  face = 'b';  clr = 'b';
                    lbl = sprintf('Peak (l=%d, v=%.0f) mag=%.2f', pk_delay, pk_dopp, pk_mag);
                end
                stem3(pk_dopp, pk_delay, 1, marker, 'MarkerSize', 10, ...
                    'LineWidth', 2, 'MarkerFaceColor', face, 'DisplayName', lbl);
                text(pk_dopp, pk_delay, 1.05, sprintf('%.2f', pk_mag), ...
                    'FontSize', 10, 'FontWeight', 'bold', 'Color', clr);
            end
            legend('Location', 'best'); hold off;
            xlabel('Doppler [tap]'); ylabel('Delay [tap]'); zlabel('Normalized Magnitude');
            title(sprintf('Coarse: Range-Velocity @ %.0f°', ang_val));
            set(gca, 'FontSize', 11, 'FontWeight', 'Bold');
        end
    end

    %% -------------------------------------------------------------------
    %  Stage 2b: pilot-arm CLEAN (snap on arm spectra, coherent local LS grid, select by LS improvement)
    % --------------------------------------------------------------------
    blocks = params.pilot_blocks;
    if ~all(ismember({'delay', 'doppler'}, {blocks.dim}))
        error('estimate_coarse_params_fused: requires delay and doppler blocks in params.pilot_blocks_def.');
    end
    sites = pilot_sites(blocks);

    % Comb ambiguity periods from the pilot spacing (contiguous arms: M, N)
    K_r = 1;  K_v = 1;
    sp_r = [];  sp_v = [];
    for b = 1:length(blocks)
        if strcmp(blocks(b).dim, 'delay'),   sp_r = [sp_r, diff(blocks(b).m(:)')]; end %#ok<AGROW>
        if strcmp(blocks(b).dim, 'doppler'), sp_v = [sp_v, diff(blocks(b).n(:)')]; end %#ok<AGROW>
    end
    if ~isempty(sp_r), K_r = sp_r(1); for d = sp_r(2:end), K_r = gcd(K_r, d); end, end
    if ~isempty(sp_v), K_v = sp_v(1); for d = sp_v(2:end), K_v = gcd(K_v, d); end, end
    per_r = M / K_r;  per_v = N / K_v;
    nfft_r = 16*M;  nfft_v = 16*N;
    tau_axis = (0:nfft_r-1) * M / nfft_r;
    nu_axis  = (-nfft_v/2:nfft_v/2-1) * N / nfft_v;

    [c, s2] = pilot_noise_model(XTF, params, n0);  % y = c*Phi*gains + noise(s2)
    y = build_measurement_vector(YTF, XTF, blocks, Nr, n0);
    min_relative_gain = 0.01;
    if isfield(params, 'coarse') && isfield(params.coarse, 'min_relative_gain')
        min_relative_gain = params.coarse.min_relative_gain;
    end
    minimum_gain = min_relative_gain * real(y' * y);

    % CLEAN: snap the remaining candidates on the residual arm spectra, refine each
    % on a local grid by coherent LS improvement, accept the best, refit all gains
    [dr, dv] = ndgrid(-0.5:1/16:0.5);  % local refinement grid offsets
    remaining = 1:n_detected;
    ang = [];  pathDelays = [];  pathDopplers = [];  history = {};
    A = zeros(length(y), 0);
    residual = y;
    while ~isempty(remaining)
        [P_r, P_v] = arm_spectra(residual, blocks, sites, M, N, Nr, nfft_r, nfft_v);
        pk_r = (P_r >= circshift(P_r, 1)) & (P_r >= circshift(P_r, -1));
        pk_v = (P_v >= circshift(P_v, 1)) & (P_v >= circshift(P_v, -1));

        % Snap each candidate to the strongest residual-spectrum peak within +-1 tap
        % (cand_r/cand_v then hold the coherently refined delay/Doppler)
        cand_r = all_delays(remaining);  cand_v = all_dopplers(remaining);
        best_gain = -inf;  best = 1;
        for i = 1:length(remaining)
            off_r = mod(tau_axis(:) - cand_r(i) + per_r/2, per_r) - per_r/2;
            [~, k] = max(P_r .* (pk_r & abs(off_r) <= 1));
            if pk_r(k) && abs(off_r(k)) <= 1
                cand_r(i) = max(cand_r(i) + off_r(k), 0);
            end
            off_v = mod(nu_axis(:) - cand_v(i) + per_v/2, per_v) - per_v/2;
            [~, k] = max(P_v .* (pk_v & abs(off_v) <= 1));
            if pk_v(k) && abs(off_v(k)) <= 1
                cand_v(i) = cand_v(i) + off_v(k);
            end

            % Coherent refinement: Schur improvement on a local delay/Doppler grid
            g = unique([max(cand_r(i) + dr(:), 0), cand_v(i) + dv(:)], 'rows');
            grid = struct('AoAs', repmat(all_angles(remaining(i)), 1, size(g,1)), ...
                'pathDelays', g(:,1).', 'pathDopplers', g(:,2).');
            grid.AoDs = grid.AoAs;  % monostatic
            Phi = c * build_response_matrix(grid, params, derivedParams);
            improvement = schur_gain(A, residual, Phi, s2);
            % Ridge regularization must not reward duplicating an existing path.
            for j = 1:length(ang)
                duplicate = abs(grid.AoAs - ang(j)) < 1e-8 & ...
                    abs(grid.pathDelays - pathDelays(j)) < 1e-8 & ...
                    abs(grid.pathDopplers - pathDopplers(j)) < 1e-8;
                improvement(duplicate) = -inf;
            end
            [gain, k] = max(improvement);
            cand_r(i) = g(k,1);  cand_v(i) = g(k,2);
            if gain > best_gain, best_gain = gain;  best = i;  phi_best = Phi(:,k); end
        end
        if best_gain <= minimum_gain, break; end

        ang(end+1)          = all_angles(remaining(best)); %#ok<AGROW>
        pathDelays(end+1)   = cand_r(best); %#ok<AGROW>
        pathDopplers(end+1) = cand_v(best); %#ok<AGROW>
        history{end+1} = struct('P_r', P_r, 'P_v', P_v, ...
            'delay', pathDelays(end), 'doppler', pathDopplers(end)); %#ok<AGROW>
        remaining(best) = [];
        A = [A, phi_best]; %#ok<AGROW>
        beta = (A'*A + s2*eye(size(A,2))) \ (A'*y);
        residual = y - A*beta;
    end

    % Nothing supported: keep the strongest xcorr candidate
    if isempty(ang)
        [~, i] = max(all_peak_vals);
        ang = all_angles(i);  pathDelays = all_delays(i);  pathDopplers = all_dopplers(i);
    end
    n_detected = length(ang);

    if params.show_figures && ~isempty(history)
        n_it = length(history);
        figure('Position', [100, 100, 1000, max(300, 200*n_it)]);
        for it = 1:n_it
            h = history{it};
            subplot(n_it, 2, 2*it-1); hold on;
            plot(tau_axis, h.P_r / max(h.P_r), '-b', 'LineWidth', 1.5);
            xline(chanParams.pathDelays, 'r--');  xline(h.delay, 'g--', 'LineWidth', 2);
            hold off; xlabel('Delay [tap]'); ylabel('Norm. Power');
            title(sprintf('CLEAN %d: delay', it));
            subplot(n_it, 2, 2*it); hold on;
            plot(nu_axis, h.P_v / max(h.P_v), '-b', 'LineWidth', 1.5);
            xline(chanParams.pathDopplers, 'r--');  xline(h.doppler, 'g--', 'LineWidth', 2);
            hold off; xlabel('Doppler [tap]'); ylabel('Norm. Power');
            title(sprintf('CLEAN %d: Doppler', it));
        end
        sgtitle('Coarse (fused): pilot-arm CLEAN spectra');
    end

    %% -------------------------------------------------------------------
    %  Stage 3: ridge LS path gains (c*Phi, ridge s2)
    % --------------------------------------------------------------------
    coarse_for_phi.AoAs         = ang;
    coarse_for_phi.AoDs         = ang;  % monostatic
    coarse_for_phi.pathDelays   = pathDelays;
    coarse_for_phi.pathDopplers = pathDopplers;

    Phi = c * build_response_matrix(coarse_for_phi, params, derivedParams);
    pathGains = (Phi' * Phi + s2 * eye(n_detected)) \ (Phi' * y);

    %% -------------------------------------------------------------------
    %  Output
    % --------------------------------------------------------------------
    coarse.AoAs         = ang;
    coarse.AoDs         = ang;  % monostatic
    coarse.pathDelays   = pathDelays;
    coarse.pathDopplers = pathDopplers;
    coarse.pathGains    = pathGains;
    coarse.n_detected   = n_detected;

    if params.verbose_output
        fprintf('--- Coarse Estimation Results ---\n');
        if params.skip_aoa
            fprintf('Cor. AoA:     unresolvable (Nr=1)\n');
        else
            fprintf('Cor. AoA:     %s\n', sprintf('%5.2f ', coarse.AoAs));
        end
        fprintf('Cor. Delay:   %s\n', sprintf('%5.2f ', coarse.pathDelays));
        fprintf('Cor. Doppler: %s\n', sprintf('%5.2f ', coarse.pathDopplers));
        fprintf('---------------------------------\n');
    end
end

function improvement = schur_gain(A, residual, Phi, ridge)
%SCHUR_GAIN Regularized LS improvement from adding each column of Phi to A and refitting ALL gains.
    G = A'*A + ridge*eye(size(A,2));
    cross = A'*Phi;
    norms = sum(abs(Phi).^2, 1);
    innovation = norms + ridge - real(sum(conj(cross).*(G\cross), 1));
    valid = innovation > 1e-12*max(1, norms);
    improvement = -inf(1, size(Phi,2));
    improvement(valid) = abs(Phi(:,valid)'*residual).'.^2 ./ innovation(valid);
end

function [P_r, P_v] = arm_spectra(residual, blocks, sites, M, N, Nr, nfft_r, nfft_v)
%ARM_SPECTRA Delay/Doppler spectra of the pilot arms, summed incoherently over Rx and blocks.
%   residual has the build_measurement_vector layout (Nr rows per unique site).
    P_r = zeros(nfft_r, 1);  P_v = zeros(nfft_v, 1);
    for b = 1:length(blocks)
        blk = blocks(b);
        values = residual((sites.pilot_site{b}-1)*Nr + (1:Nr));  % n_pilots x Nr
        if strcmp(blk.dim, 'delay')
            s = zeros(M, Nr);  s(blk.m, :) = values;
            P_r = P_r + sum(abs(ifft(s, nfft_r)).^2, 2);
        elseif strcmp(blk.dim, 'doppler')
            s = zeros(N, Nr);  s(blk.n, :) = values;
            P_v = P_v + sum(abs(fftshift(fft(s, nfft_v), 1)).^2, 2);
        end
    end
end
