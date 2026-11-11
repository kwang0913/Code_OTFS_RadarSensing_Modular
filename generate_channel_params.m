function chanParams = generate_channel_params(params)
%GENERATE_CHANNEL_PARAMS Generate monostatic target parameters.
%   Handles three target selection modes:
%     Mode 3 (deterministic): params.targets defines exact targets
%     Mode 2 (random):        params.n_targets + aoa/delay/doppler ranges
%     Mode 1 (pool):          params.selected_targets from params.pool
%   epsilon, kappa, and gain are applied as post-processing on any mode.
%   Angles, delays and Dopplers are rounded to params.angle/delay/doppler_resolution
%   (before the epsilon/kappa overrides, which stay exact).
%   Monostatic: AoDs = AoAs. The guard length is fixed by params.padLen; the
%   largest delay must fit in it.
%
%   RNG draw order: targets (selection, angles, delays) -> gains -> Dopplers
%   (magnitude, sign); pilot positions, data, guards, pilot values and noise follow.

    %% Target selection (angles, delays)
    draw_doppler = false;
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
            delays = sample_tap_range(params.delay_range, K, 1, 'delay_range', params.delay_resolution);
        end
        draw_doppler = ~isempty(params.doppler_range);  % drawn after the gains

    else
        % --- Mode 1: pool selection ---
        target_list = params.pool(params.selected_targets);
        angles   = [target_list.angle];
        delays   = [target_list.delay];
        dopplers = [target_list.doppler];
        K = numel(params.selected_targets);
    end
    angles = snap(angles, params.angle_resolution);
    delays = snap(reshape(delays, 1, []), params.delay_resolution);

    %% Post-processing: override fractional delay part (for controlled experiments)
    if ~isempty(params.epsilon)
        epsilon = reshape(params.epsilon, 1, []);
        if (~isscalar(epsilon) && numel(epsilon) ~= K) || ...
                ~isreal(epsilon) || any(~isfinite(epsilon) | epsilon < 0 | epsilon >= 1)
            error('epsilon must be a scalar or one value per target in [0,1).');
        end
        delays = floor(delays) + epsilon;
    end
    if ~isreal(delays) || any(~isfinite(delays) | delays < 0 | delays >= params.M)
        error('Target delays must be finite and in [0,M).');
    end
    if ceil(max(delays)) > params.padLen
        error('generate_channel_params:padLen', ...
            'ceil(max delay) = %d exceeds params.padLen = %d.', ceil(max(delays)), params.padLen);
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
        gains = sqrt(1/2) * (randn(1, K) + 1j * randn(1, K));  % complex Rayleigh
        % gains = exp(1j*2*pi*rand(1, K));                     % unit modulus, random phase
    end

    %% Dopplers (random range: magnitude, then sign)
    if draw_doppler
        dop_mag = sample_tap_range(params.doppler_range, K, 0, 'doppler_range', params.doppler_resolution);
        dop_sign = 2*randi([0,1], 1, K) - 1;
        dopplers = dop_mag .* dop_sign;
    end
    dopplers = snap(dopplers, params.doppler_resolution);

    %% Post-processing: override fractional Doppler part (for controlled experiments)
    if ~isempty(params.kappa)
        dopplers = floor(dopplers) + params.kappa;
    end

    %% Assemble chanParams
    chanParams.AoAs         = angles;
    chanParams.AoDs         = angles;  % monostatic
    chanParams.pathDelays   = delays;
    chanParams.pathDopplers = dopplers;
    chanParams.pathGains    = gains;

    %% Verbose: print ground truth
    if isfield(params, 'verbose_output') && params.verbose_output
        derivedParams = compute_derived_params(params);
        fprintf('--- Ground Truth ---\n');
        fprintf('  Angles:     %s deg\n', sprintf('%.1f  ', angles));
        fprintf('  Delays:     %s taps\n', sprintf('%.6g  ', delays));
        fprintf('  Dopplers:   %s taps\n', sprintf('%.2f  ', dopplers));
        fprintf('  Ranges:     %s m\n', sprintf('%.1f  ', delays * params.c / (2 * derivedParams.fsamp)));
        fprintf('  Velocities: %s m/s\n', sprintf('%.1f  ', dopplers / derivedParams.Tsamp * derivedParams.lambda / 2));
        fprintf('  PathGains:  %s\n', sprintf('%.1f  ', abs(gains)));
        fprintf('--------------------\n');
    end
end

function x = snap(x, res)
% Round to the resolution grid (12-digit rounding removes multiplication noise).
    x = round(round(x / res) * res, 12);
end

function taps = sample_tap_range(bounds, count, scalar_min, field, res)
% Draw uniformly from the res-spaced grid inside the requested bounds.
    bounds = reshape(bounds, 1, []);
    if isscalar(bounds)
        bounds = [scalar_min, bounds];
    end
    if numel(bounds) ~= 2 || ~isreal(bounds) || ...
            any(~isfinite(bounds) | bounds < 0) || bounds(1) > bounds(2)
        error('%s must specify a nonnegative [min,max] range.', field);
    end
    ticks = round(double(bounds) / res, 12);
    first = ceil(ticks(1));
    last = floor(ticks(2));
    if first > last
        error('generate_channel_params:emptyTapGrid', ...
            '%s contains no point on the %g-tap grid.', field, res);
    end
    taps = round(randi([first, last], 1, count) * res, 12);
end
