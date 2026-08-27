function coarse = estimate_coarse_params(YTF, YDD, XTF, XDD, params, chanParams, n0)
%ESTIMATE_COARSE_PARAMS Coarse DFRC parameter estimation for monostatic radar.
%   Stage 1: SIC AoA detection
%   Stage 2: MMSE beamforming + xcorr2 in DD domain for delay-Doppler
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
    %  Stage 2: xcorr2 for delay-Doppler (with beamforming when Nr > 1)
    % --------------------------------------------------------------------
    xcorr2_thresh = 0.1;  % fraction of per-bin peak max for local maxima detection
    doppler_taps  = -(N-1):(N-1);
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

    % xcorr2 per angle with local maxima detection
    all_angles    = [];
    all_delays    = [];
    all_dopplers  = [];
    all_peak_vals = [];

    for p = 1:n_ang
        Y_bf  = squeeze(ADD(p, :, :));
        X_ref = squeeze(XDD_ref(p, :, :));

        Hdd = xcorr2(Y_bf, X_ref);
        Hdd_pos = abs(Hdd(M:end, :));

        peak_max = max(Hdd_pos(:));
        if peak_max == 0, continue; end

        thresh = xcorr2_thresh * peak_max;
        [nR, nC] = size(Hdd_pos);

        % Local maxima: >= all 8 neighbors and above threshold
        Hpad = zeros(nR+2, nC+2);
        Hpad(2:end-1, 2:end-1) = Hdd_pos;
        is_peak = true(nR, nC);
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0, continue; end
                is_peak = is_peak & (Hdd_pos >= Hpad((2:end-1)+dr, (2:end-1)+dc));
            end
        end
        is_peak = is_peak & (Hdd_pos >= thresh);

        [r_idx, c_idx] = find(is_peak);
        for k = 1:length(r_idx)
            all_angles(end+1)    = ang(p); %#ok<AGROW>
            all_delays(end+1)    = r_idx(k) - 1; %#ok<AGROW>
            all_dopplers(end+1)  = doppler_taps(c_idx(k)); %#ok<AGROW>
            all_peak_vals(end+1) = Hdd_pos(r_idx(k), c_idx(k)); %#ok<AGROW>
        end
    end

    % Fallback: if no peaks found, use global max of first angle bin
    if isempty(all_delays) && n_ang > 0
        Y_bf  = squeeze(ADD(1, :, :));
        X_ref = squeeze(XDD_ref(1, :, :));
        Hdd = xcorr2(Y_bf, X_ref);
        Hdd_pos = abs(Hdd(M:end, :));
        [~, pk_lin] = max(Hdd_pos(:));
        [row_pk, col_pk] = ind2sub(size(Hdd_pos), pk_lin);
        all_angles    = ang(1);
        all_delays    = row_pk - 1;
        all_dopplers  = doppler_taps(col_pk);
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
            Hdd_p = xcorr2(Y_bf_p, X_ref_p);

            Hdd_pos_plot = abs(Hdd_p(M:end, :));
            Hdd_pos_plot = Hdd_pos_plot / max(Hdd_pos_plot(:));

            delay_axis = 0:(M-1);

            figure;
            mesh(doppler_taps, delay_axis, Hdd_pos_plot, 'EdgeColor', 'interp', 'FaceColor', 'interp');
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
            [nR_p, nC_p] = size(Hdd_pos_plot);
            Hpad_p = zeros(nR_p+2, nC_p+2);
            Hpad_p(2:end-1, 2:end-1) = Hdd_pos_plot;
            is_pk = true(nR_p, nC_p);
            for dr = -1:1
                for dc = -1:1
                    if dr == 0 && dc == 0, continue; end
                    is_pk = is_pk & (Hdd_pos_plot >= Hpad_p((2:end-1)+dr, (2:end-1)+dc));
                end
            end
            is_pk = is_pk & (Hdd_pos_plot >= xcorr2_thresh);
            [pk_r, pk_c] = find(is_pk);

            % Plot all peaks above threshold
            est_set_delays   = pathDelays(tgt_mask);
            est_set_dopplers = pathDopplers(tgt_mask);
            for qi = 1:length(pk_r)
                pk_delay = pk_r(qi) - 1;
                pk_dopp  = doppler_taps(pk_c(qi));
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
