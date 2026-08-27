function chanParams = generate_channel_params(params)
%GENERATE_CHANNEL_PARAMS Generate monostatic target parameters.
%   Handles three target selection modes:
%     Mode 3 (deterministic): params.targets defines exact targets
%     Mode 2 (random):        params.n_targets + aoa/delay/doppler ranges
%     Mode 1 (pool):          params.selected_targets from params.pool
%   kappa and gain are applied as post-processing on any mode.
%   Path gains are random complex Rayleigh unless params.gain is set.

    c = params.c;
    lambda = params.lambda;
    fsamp = params.fsamp;

    %% Target selection
    if ~isempty(params.targets)
        % --- Mode 3: deterministic ---
        angles   = [params.targets.angle];
        delays   = [params.targets.delay];
        dopplers = [params.targets.doppler];
        K = numel(params.targets);

    elseif ~isempty(params.n_targets)
        % --- Mode 2: random n_targets ---
        K = params.n_targets;

        % Select from pool, then override parameters if ranges are specified
        sel = sort(randperm(numel(params.pool), K));
        target_list = params.pool(sel);
        angles   = [target_list.angle];
        delays   = [target_list.delay];
        dopplers = [target_list.doppler];

        if ~isempty(params.aoa_range)
            sel_aoa = randperm(length(params.aoa_range), K);
            angles = params.aoa_range(sel_aoa);
        end
        if ~isempty(params.delay_range)
            delays = randi(params.delay_range, 1, K);  % delay_range must be scalar (max) or [lo, hi]
        end
        if ~isempty(params.doppler_range)
            dop_mag  = params.doppler_range(1) + diff(params.doppler_range) * rand(1, K);
            dop_sign = 2*randi([0,1], 1, K) - 1;
            dopplers = dop_mag .* dop_sign;
        end

    else
        % --- Mode 1: pool selection ---
        target_list = params.pool(params.selected_targets);
        angles   = [target_list.angle];
        delays   = [target_list.delay];
        dopplers = [target_list.doppler];
        K = numel(params.selected_targets);
    end

    %% Post-processing: override fractional Doppler part (for controlled experiments)
    if ~isempty(params.kappa)
        dopplers = floor(dopplers) + params.kappa;
    end

    %% Path gains
    if isfield(params, 'gain') && ~isempty(params.gain)
        gains = params.gain;
        if isscalar(gains)
            gains = repmat(gains, 1, K);
        elseif numel(gains) == K
            gains = reshape(gains, 1, []);
        else
            error('generate_channel_params:badGainLength', ...
                'params.gain must be scalar or have one entry per target.');
        end
    else
        gains = sqrt(1/2) * (randn(1, K) + 1j * randn(1, K));
    end

    %% Assemble chanParams
    chanParams.AoAs         = angles;
    chanParams.pathDelays   = delays;
    chanParams.pathDopplers = dopplers;
    chanParams.pathGains    = gains;

    %% Derived quantities
    chanParams.padLen = max(delays);
    [chanParams.Meff, chanParams.numSamps, chanParams.offset] = ...
        compute_padding_params(params.M, params.N, chanParams.padLen, params.padType);
    chanParams.fsamp = fsamp;
    chanParams.Tsamp = chanParams.Meff / fsamp;

    chanParams.ranges           = delays * c / (2 * fsamp);
    chanParams.pathDopplerFreqs = dopplers / (params.N * chanParams.Tsamp);
    chanParams.velocities       = chanParams.pathDopplerFreqs * lambda / 2;

    %% Verbose: print ground truth
    if isfield(params, 'verbose_output') && params.verbose_output
        fprintf('--- Ground Truth ---\n');
        fprintf('  Angles:     %s deg\n', sprintf('%.1f  ', angles));
        fprintf('  Delays:     %s taps\n', sprintf('%d  ', delays));
        fprintf('  Dopplers:   %s taps\n', sprintf('%.2f  ', dopplers));
        fprintf('  Ranges:     %s m\n', sprintf('%.1f  ', chanParams.ranges));
        fprintf('  Velocities: %s m/s\n', sprintf('%.1f  ', chanParams.velocities));
        g = gains;
        fprintf('  PathGains:  %s\n', sprintf('%.1f  ', abs(g)));
        fprintf('--------------------\n');
    end
end
