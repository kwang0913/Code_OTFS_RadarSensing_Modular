function coarse = estimate_coarse_params_new(YTF, YDD, XTF, XDD, params, chanParams, n0)
%ESTIMATE_COARSE_PARAMS Coarse DFRC parameter estimation for monostatic radar.
%   Stage 1: SIC AoA detection
%   Stage 2: MMSE beamforming + DFT-based delay-Doppler estimation
%   Stage 3: LS path gains

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
    %  Stage 2: TF-grid DFT delay-Doppler estimation
    % --------------------------------------------------------------------
    peak_thresh = 0.1;   % fraction of per-spectrum peak for peak detection
    dftM        = 4 * M; % zero-padded delay DFT size
    dftN        = 4 * N; % zero-padded Doppler DFT size
    n_ang       = length(ang);
    XTF_flat    = reshape(permute(XTF, [3,1,2]), Nt, []);  % Nt x (M*N)

    if params.skip_aoa
        % Nr=1: no beamforming, single TF grid
        ATF = reshape(squeeze(YTF), 1, M, N);
        XTF_ref = reshape(sum(XTF, 3), 1, M, N);
    else
        % MMSE beamforming weights
        Ang_mtx = exp(1j * 2*pi * (0:Nr-1)' * params.ddlr * sind(ang));
        W_bf = (Ang_mtx' * Ang_mtx + n0 * eye(n_ang)) \ Ang_mtx';

        % Per-angle TF grids via MMSE beamforming
        ATF = reshape(W_bf * reshape(permute(YTF, [3,1,2]), Nr, []), n_ang, M, N);

        % Tx-steered TF reference per angle
        XTF_ref = zeros(n_ang, M, N);
        for p = 1:n_ang
            d_tx = exp(-1j*2*pi*(0:Nt-1)'*params.ddlt*sind(ang(p)));
            XTF_ref(p,:,:) = reshape(d_tx.' * XTF_flat, M, N);
        end
    end

    % --- Build private-bin mask and skip full rows/columns that touch it ---
    private_mask = false(M, N);
    for tx = 1:Nt
        private_mask = private_mask | (XDD(:, :, tx) == 0);
    end
    for b = 1:length(params.pilot_blocks)
        blk = params.pilot_blocks(b);
        for k = 1:blk.n_pilots
            private_mask(blk.m(k), blk.n(k)) = true;
        end
    end
    valid_delay_cols = ~any(private_mask, 1);
    valid_doppler_rows = ~any(private_mask, 2);
    if ~any(valid_delay_cols), valid_delay_cols(:) = true; end
    if ~any(valid_doppler_rows), valid_doppler_rows(:) = true; end

    all_angles    = [];
    all_delays    = [];
    all_dopplers  = [];
    all_peak_vals = [];
    delay_spectra = zeros(n_ang, dftM);
    doppler_spectra = zeros(n_ang, dftN);

    reg_eps = 0.01;
    max_per_angle = max(1, ceil(max_targets / max(1, n_ang)));
    min_delay_sep = max(1, round(dftM / M));
    min_doppler_sep = max(1, round(dftN / N));

    for p = 1:n_ang
        X_ref = squeeze(XTF_ref(p, :, :));
        H = squeeze(ATF(p, :, :)) .* conj(X_ref) ./ (abs(X_ref).^2 + reg_eps);
        H(private_mask) = 0;

        delay_map = abs(ifft(H(:, valid_delay_cols), dftM, 1));
        delay_spec = mean(delay_map, 2);
        doppler_map = abs(fft(H(valid_doppler_rows, :), dftN, 2));
        doppler_spec = mean(doppler_map, 1).';

        if max(delay_spec) > 0
            delay_spectra(p, :) = (delay_spec / max(delay_spec)).';
        end
        if max(doppler_spec) > 0
            doppler_spectra(p, :) = (doppler_spec / max(doppler_spec)).';
        end

        [d_pk, d_loc] = circular_top_peaks(delay_spec, max_per_angle, peak_thresh, min_delay_sep);
        [v_pk, v_loc] = circular_top_peaks(doppler_spec, max_per_angle, peak_thresh, min_doppler_sep);
        n_pair = min([length(d_loc), length(v_loc), max_per_angle]);

        for q = 1:n_pair
            delay_tap = round((d_loc(q) - 1) * M / dftM);
            delay_tap = min(max(delay_tap, 0), M - 1);

            doppler_bin = v_loc(q) - 1;
            if doppler_bin >= dftN/2
                doppler_bin = doppler_bin - dftN;
            end
            doppler_tap = round(doppler_bin * N / dftN);

            all_angles(end+1)    = ang(p); %#ok<AGROW>
            all_delays(end+1)    = delay_tap; %#ok<AGROW>
            all_dopplers(end+1)  = doppler_tap; %#ok<AGROW>
            all_peak_vals(end+1) = sqrt(d_pk(q) * v_pk(q)); %#ok<AGROW>
        end
    end

    % --- Deduplicate: for shared (delay, Doppler) peaks at close angles, keep strongest ---
    if params.skip_aoa
        beamwidth_deg = inf;
    else
        beamwidth_deg = asind(min(1, 1 / (Nr * params.ddlr)));
    end
    ang_dedup_thresh = beamwidth_deg;

    n_peaks = length(all_delays);
    keep = true(1, n_peaks);
    for i = 1:n_peaks
        if ~keep(i), continue; end
        for j = i+1:n_peaks
            if ~keep(j), continue; end
            if all_delays(i) == all_delays(j) && all_dopplers(i) == all_dopplers(j) ...
                    && abs(all_angles(i) - all_angles(j)) < ang_dedup_thresh
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

    if length(all_peak_vals) > max_targets
        [~, ord] = sort(all_peak_vals, 'descend');
        ord = ord(1:max_targets);
        all_angles    = all_angles(ord);
        all_delays    = all_delays(ord);
        all_dopplers  = all_dopplers(ord);
        all_peak_vals = all_peak_vals(ord);
    end

    n_detected   = length(all_delays);
    ang          = all_angles;
    pathDelays   = all_delays;
    pathDopplers = all_dopplers;

    % Plot TF DFT delay and Doppler spectra
    if params.show_figures
        delay_axis = (0:dftM-1) * M / dftM;
        doppler_axis = (0:dftN-1) * N / dftN;
        doppler_axis(doppler_axis >= N/2) = doppler_axis(doppler_axis >= N/2) - N;

        unique_ang = unique(ang, 'stable');
        for ab = 1:min(length(unique_ang), 4)
            ang_val = unique_ang(ab);
            [~, bf_idx] = min(abs(detected_angles - ang_val));
            tgt_mask = (ang == ang_val);

            figure; hold on;
            plot(delay_axis, delay_spectra(bf_idx, :), '-b', 'LineWidth', 1.5);
            for tidx = 1:length(chanParams.AoAs)
                [~, nearest_bin] = min(abs(unique_ang - chanParams.AoAs(tidx)));
                if nearest_bin == ab
                    xline(chanParams.pathDelays(tidx), 'r--', 'LineWidth', 2);
                end
            end
            est_set = unique(pathDelays(tgt_mask));
            for qi = 1:length(est_set)
                xline(est_set(qi), 'g:', 'LineWidth', 2);
            end
            hold off;
            xlabel('Delay [tap]'); ylabel('Normalized Magnitude');
            title(sprintf('Coarse: TF Delay Spectrum @ %.0f°', ang_val));
            set(gca, 'FontSize', 11, 'FontWeight', 'Bold');

            figure; hold on;
            plot(doppler_axis, doppler_spectra(bf_idx, :), '-b', 'LineWidth', 1.5);
            for tidx = 1:length(chanParams.AoAs)
                [~, nearest_bin] = min(abs(unique_ang - chanParams.AoAs(tidx)));
                if nearest_bin == ab
                    xline(chanParams.pathDopplers(tidx), 'r--', 'LineWidth', 2);
                end
            end
            est_set = unique(pathDopplers(tgt_mask));
            for qi = 1:length(est_set)
                xline(est_set(qi), 'g:', 'LineWidth', 2);
            end
            hold off;
            xlabel('Doppler [tap]'); ylabel('Normalized Magnitude');
            title(sprintf('Coarse: TF Doppler Spectrum @ %.0f°', ang_val));
            set(gca, 'FontSize', 11, 'FontWeight', 'Bold');
        end
    end

    %% -------------------------------------------------------------------
    %  Stage 3: LS path gains
    % --------------------------------------------------------------------
    coarse_for_phi.AoAs         = ang;
    coarse_for_phi.pathDelays   = pathDelays;
    coarse_for_phi.pathDopplers = pathDopplers;

    Phi = build_response_matrix(coarse_for_phi, params);
    y   = build_measurement_vector(YTF, XTF, params.pilot_blocks, Nr);
    pathGains = (Phi' * Phi + n0 * eye(n_detected)) \ (Phi' * y);

    %% -------------------------------------------------------------------
    %  Output
    % --------------------------------------------------------------------
    coarse.AoAs         = ang;
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

function [pk_vals, pk_locs] = circular_top_peaks(spec, max_count, thresh_frac, min_sep)
%CIRCULAR_TOP_PEAKS Top circular local maxima with endpoint support.
    spec = spec(:);
    if isempty(spec) || max_count <= 0 || all(spec == 0)
        pk_vals = [];
        pk_locs = [];
        return;
    end

    thresh = thresh_frac * max(spec);
    is_peak = spec >= circshift(spec, 1) & ...
              spec >= circshift(spec, -1) & ...
              spec >= thresh;
    cand_locs = find(is_peak);

    if isempty(cand_locs)
        [pk_vals, pk_locs] = max(spec);
        return;
    end

    cand_vals = spec(cand_locs);
    [cand_vals, ord] = sort(cand_vals, 'descend');
    cand_locs = cand_locs(ord);

    pk_vals = [];
    pk_locs = [];
    L = length(spec);
    for ii = 1:length(cand_locs)
        loc = cand_locs(ii);
        if isempty(pk_locs)
            keep = true;
        else
            dist = abs(pk_locs - loc);
            dist = min(dist, L - dist);
            keep = all(dist >= min_sep);
        end
        if keep
            pk_vals(end+1) = cand_vals(ii); %#ok<AGROW>
            pk_locs(end+1) = loc; %#ok<AGROW>
            if length(pk_locs) >= max_count
                break;
            end
        end
    end
end
