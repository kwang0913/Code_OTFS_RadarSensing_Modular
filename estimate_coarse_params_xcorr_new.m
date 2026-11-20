function coarse = estimate_coarse_params_xcorr_new(YTF, YDD, XTF, XDD, params, chanParams, n0, derivedParams)
%ESTIMATE_COARSE_PARAMS_XCORR_NEW Coarse DFRC parameter estimation, joint virtual-array angle.
%   Stage 1: per-(Rx,Tx) circular DD xcorr C_{r,p}; angle-independent map P = sum |C_{r,p}|^2
%   Stage 2: (delay, Doppler) cells = circular 8-neighbour local maxima of P above a
%            dynamic-range and a noise-floor threshold (a neighbour only masks a cell if it
%            carries the same angular signature), strongest first, at most targets_num
%   Stage 3: per cell, angles from the Nr*Nt values C_{r,p}(cell) on the virtual array
%            (position ddlr*r - ddlt*p) by grid matched filter + OMP (joint LS refit)
%   Stage 4: ridge LS path gains, model y = c*Phi*gains + noise (pilot_noise_model)
%   No Rx beamforming, no SIC on the Rx array. Angles snap to the SSR angle min step.

    M  = params.M;
    N  = params.N;
    Nr = params.Nr;
    Nt = params.Nt;
    L  = Nr * Nt;                  % virtual array size

    % --- Named constants, calibrated on separate seeds (101-105), single and 3-target scenes, 20 and -10 dB, ddlt 8 and 0.5 ---
    PEAK_REL      = 0.01;   % cell power >= 0.01*max(P)  (= amplitude >= 0.1 of the peak)
    SAME_PATH     = 0.95;   % |cos| between the (r,p) vectors of neighbouring cells above which they are one path:
                            % skirts of one path >= 0.988; adjacent cells of paths at other angles <= 0.887
                            % (2 deg apart at ddlt=0.5, the unresolvable worst case), <= 0.05 at ddlt=8
    MAX_CAND      = 400;    % strongest cells above both thresholds that are examined
    NOISE_K       = 20;     % cell power >= NOISE_K*median(P): worst non-target local max 9.3x (31.6x only as
                            % target skirts, < 0.003*max(P)); weakest target cell >= 346x at -10 dB
    ANGLE_GRID    = 4096;   % angle grid points, uniform in sin(theta) over [-1,1) (~0.03 deg)
    MAX_ANGLES    = 3;      % angles per cell
    BACKFIT       = 3;      % re-pick sweeps over the accepted angles (greedy picks of close angles are biased)
    ANGLE_REL     = 0.1;    % further angle needs residual amplitude >= ANGLE_REL * first amplitude
                            % (single-target LS residual: <= 0.019 over all calibration runs)
    ANGLE_NOISE_K = 25;     % ... and residual power >= ANGLE_NOISE_K * matched-filter noise power
                            % (pure-noise residual peaks: <= 9.6x at -10 dB)
    max_cells     = params.targets_num;

    wrapN = @(v) N/2 - mod(N/2 - v, N);  % circular Doppler lag -> (-N/2, N/2]

    %% Stage 1: per-pair xcorr and angle-independent map
    FY = fft2(YDD);                      % M x N x Nr
    FX = fft2(XDD);                      % M x N x Nt
    C  = zeros(M, N, Nr, Nt);            % row -> delay row-1, column -> Doppler lag col-1
    P  = zeros(M, N);
    for p = 1:Nt
        Cp = ifft2(FY .* conj(FX(:, :, p)));
        C(:, :, :, p) = Cp;
        P = P + sum(abs(Cp).^2, 3);
    end

    %% Stage 2: (delay, Doppler) cells
    %  Peak = above both thresholds and not dominated by a larger 8-neighbour of the SAME path
    %  (coherent angular signature); a larger neighbour at another angle does not mask it.
    Pmed = median(P(:));
    cand = find((P >= PEAK_REL * max(P(:))) & (P >= NOISE_K * Pmed));
    [~, order] = sort(P(cand), 'descend');
    cand = cand(order(1:min(numel(order), MAX_CAND)));
    Vc = reshape(C, M * N, L);           % row = cell, column = (r,p) pair
    is_peak = true(size(cand));
    for i = 1:numel(cand)
        [l, k] = ind2sub([M N], cand(i));
        for dr = -1:1
            for dc = -1:1
                nb = sub2ind([M N], mod(l-1+dr, M) + 1, mod(k-1+dc, N) + 1);
                if P(nb) > P(cand(i)) && ...
                        abs(Vc(cand(i), :) * Vc(nb, :)') / (norm(Vc(cand(i), :)) * norm(Vc(nb, :))) > SAME_PATH
                    is_peak(i) = false;
                end
            end
        end
    end
    cells = cand(is_peak);               % strongest first
    cells = cells(1:min(numel(cells), max_cells));
    if isempty(cells)
        [~, cells] = max(P(:));          % nothing passes: global maximum
    end

    cell_delay   = mod(cells - 1, M).';
    cell_doppler = wrapN(floor((cells - 1) / M)).';

    %% Stage 3: angle per cell on the virtual array
    min_step = ssr_min_steps(params).a;
    if params.skip_aoa
        A = [];
    else
        [R, T] = ndgrid(0:Nr-1, 0:Nt-1);
        pos = params.ddlr * R(:) - params.ddlt * T(:);          % wavelengths, matches the channel signs
        sin_grid = -1 + 2 * (0:ANGLE_GRID-1).' / ANGLE_GRID;
        A = exp(1j * 2*pi * pos * sin_grid.');                  % L x G
        sigma2 = Pmed / L;                                       % per-pair noise power
        noise_pow = sigma2 / L;                                  % matched-filter noise power
    end

    est = zeros(0, 3);                   % [angle delay doppler]
    spec1 = cell(numel(cells), 1);
    for c = 1:numel(cells)
        [l, k] = ind2sub([M N], cells(c));
        delay   = cell_delay(c);
        doppler = cell_doppler(c);
        if params.skip_aoa
            angles = 0;
        else
            v = reshape(C(l, k, :, :), [], 1);
            sel = [];  resid = v;
            for it = 1:MAX_ANGLES
                spec = abs(A' * resid) / L;
                [s, g] = max(spec);
                if it == 1
                    s1 = s;  spec1{c} = spec;
                elseif s < ANGLE_REL * s1 || s^2 < ANGLE_NOISE_K * noise_pow
                    break;
                end
                sel(end+1) = g; %#ok<AGROW>
                for sweep = 1:BACKFIT * (numel(sel) > 1)
                    for j = 1:numel(sel)
                        rest = sel([1:j-1, j+1:end]);
                        [~, sel(j)] = max(abs(A' * (v - A(:, rest) * (A(:, rest) \ v))));
                    end
                end
                resid = v - A(:, sel) * (A(:, sel) \ v);
            end
            angles = round(round(asind(sin_grid(sel)) / min_step) * min_step, 12);
        end
        est = [est; angles(:), repmat([delay doppler], numel(angles), 1)]; %#ok<AGROW>
    end
    est = unique(est, 'rows', 'stable');
    ang          = est(:, 1).';
    pathDelays   = est(:, 2).';
    pathDopplers = est(:, 3).';
    n_detected   = numel(ang);

    if params.show_figures
        figure;
        imagesc(-N/2:N/2-1, 0:M-1, 10*log10(fftshift(P, 2) + eps)); axis xy; hold on;
        plot(cell_doppler, cell_delay, 'wo', 'MarkerSize', 12, 'LineWidth', 2);
        colorbar; xlabel('Doppler [tap]'); ylabel('Delay [tap]');
        title('Coarse (xcorr\_new): P = \Sigma|C_{r,p}|^2 [dB], detected cells');
        if ~params.skip_aoa
            figure;
            for c = 1:numel(cells)
                subplot(numel(cells), 1, c); hold on;
                plot(asind(sin_grid), spec1{c} / max(spec1{c}), '-b', 'LineWidth', 1.5);
                for j = 1:length(chanParams.AoAs), xline(chanParams.AoAs(j), 'r--'); end
                for a = ang(pathDelays == cell_delay(c) & pathDopplers == cell_doppler(c))
                    xline(a, 'g-', 'LineWidth', 2);
                end
                ylabel(sprintf('cell (%d,%g)', cell_delay(c), cell_doppler(c)));
            end
            xlabel('Angle [\circ]'); sgtitle('Coarse (xcorr\_new): per-cell angle spectrum');
        end
    end

    %% Stage 4: ridge LS path gains (c*Phi, ridge s2)
    coarse_for_phi.AoAs         = ang;
    coarse_for_phi.AoDs         = ang;  % monostatic
    coarse_for_phi.pathDelays   = pathDelays;
    coarse_for_phi.pathDopplers = pathDopplers;

    [c, s2] = pilot_noise_model(XTF, params, n0);
    Phi = c * build_response_matrix(coarse_for_phi, params, derivedParams);
    y   = build_measurement_vector(YTF, XTF, params.pilot_blocks, Nr, n0);
    pathGains = (Phi' * Phi + s2 * eye(n_detected)) \ (Phi' * y);

    %% Output
    coarse.AoAs         = ang;
    coarse.AoDs         = ang;  % monostatic
    coarse.pathDelays   = pathDelays;
    coarse.pathDopplers = pathDopplers;
    coarse.pathGains    = pathGains;
    coarse.n_detected   = n_detected;

    if params.verbose_output
        fprintf('--- Coarse Estimation Results (xcorr_new) ---\n');
        if params.skip_aoa
            fprintf('Cor. AoA:     unresolvable (Nr=1)\n');
        else
            fprintf('Cor. AoA:     %s\n', sprintf('%5.2f ', coarse.AoAs));
        end
        fprintf('Cor. Delay:   %s\n', sprintf('%5.2f ', coarse.pathDelays));
        fprintf('Cor. Doppler: %s\n', sprintf('%5.2f ', coarse.pathDopplers));
        fprintf('---------------------------------------------\n');
    end
end
