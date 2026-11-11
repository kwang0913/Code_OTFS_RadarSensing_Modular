function params = config_params_HF(overrides)
%CONFIG_PARAMS_HF High-frequency (24.25 GHz) monostatic DFRC simulation parameters.
%   Monostatic: AoA = AoD always; round-trip delay taps = 2*range/c * fsamp.

    params.config_name = 'hf';

    % --- Control ---
    params.show_figures = true;
    params.verbose_output = true;
    params.coarse_method = 'xcorr';     % 'xcorr' (estimate_coarse_params) | 'fused' (estimate_coarse_params_fused, needs delay and doppler blocks in pilot_blocks_def)

    % --- Pilots ---
    params.pilot_power = 2;
    params.pilot_type = 'random';  % pilot values: 'random' (QPSK) | 'deterministic' (all ones)
    params.pilot_blocks_def = [];  % structured pilot blocks (empty = none)
    params.n_pilotspertx = 16;     % random pilots per Tx without a block in pilot_blocks_def: scalar or one value per Tx

    % --- Guard interval ---
    params.padType = 'RCP';         % 'CP' | 'RCP' | 'ZP' | 'RZP' | 'NONE'
    params.padLen = 20;             % guard length (samples), chosen by design; must cover the max delay

    % --- System ---
    params.Nt = 4;
    params.Nr = 16;
    params.fc = 24.25e9;  % Hz
    params.M  = 128;      % delay bins (subcarriers)
    params.N  = 64;       % Doppler bins (OTFS symbols)
    params.df = 120e3;    % subcarrier spacing (Hz)
    params.c  = 299792458; % exact speed of light in vacuum (m/s)
    params.range_factor = 2;  % range/velocity resolution factor: 1 bistatic, 2 monostatic round trip
    params.ddlt = 0.5;        % Tx antenna spacing (wavelengths)
    params.ddlr = 0.5;        % Rx antenna spacing (wavelengths)
    params.ModOrder  = 4;
    params.modulation = 'QAM';  % 'QAM' | 'PSK'

    % --- Simulation loop ---
    params.SNR_db_list = 20;          % SNR points [dB]
    params.trials = 1;
    params.rng_seed_start = 1;        % trial i uses rng(rng_seed_start + i - 1, 'threefry')
    params.n_workers = 6;             % parfor workers when trials > 1; 0 = serial

    % --- Targets: targets (explicit) > n_targets (random) > selected_targets (pool) ---
    % kappa, epsilon and gain apply to every mode; aoa/delay/doppler_range only to random mode.
    % Explicit targets, e.g.:
    % params.targets = [struct('angle',20,'delay',5,'doppler',-4.2), ...
    %                   struct('angle',35,'delay',8,'doppler',5.4), ...
    %                   struct('angle',-15,'delay',12,'doppler',3.1), ...
    %                   struct('angle',45,'delay',3,'doppler',-6.0)];
    params.targets = [];
    params.n_targets = [];    % random: n_targets drawn from the pool; aoa/delay/doppler ranges override
    params.selected_targets = [5, 10, 15, 20];  % pool indices (targets and n_targets empty)
    params.targets_file = 'targets/pool_5.json';
    params.doppler_range = [3 6];     % magnitude on the doppler_resolution grid, random sign; [] uses pool
    params.kappa = [];                 % fixed fractional Doppler
    params.aoa_range = [];            % discrete candidate angles [deg] (list, not an interval), sampled without replacement: n_targets <= numel; [] uses pool
    params.delay_range = [0 12];      % sampled on the delay_resolution grid; scalar max means [1,max], or [lo,hi]; [] uses pool
    params.epsilon = [];               % fractional delay override; scalar or one value per target in [0,1)
    params.gain = 1;                   % fixed path gain; scalar or one value per target

    % --- Resolutions: generation rounds angles/delays/Dopplers to these grids; SSR stops at them by default ---
    params.angle_resolution   = 0.1;   % deg
    params.delay_resolution   = 0.1;   % taps
    params.doppler_resolution = 0.1;   % taps

    % --- SSR: search half-width (+/- from coarse estimate), stage-1 angle step, min steps ---
    params.ssr.range.a = 5;   % angle
    params.ssr.range.r = 1;   % delay
    params.ssr.range.v = 1;   % Doppler
    % Delay/Doppler steps and the stage-2 angle step are the min steps.
    params.ssr.length.a = 1;
    params.ssr.min_step.a = [];   % [] = angle_resolution
    params.ssr.min_step.r = [];   % [] = delay_resolution
    params.ssr.min_step.v = [];   % [] = doppler_resolution
    % max_J: [] uses targets_num; a positive integer overrides that limit.
    % Residual-based stopping can stop earlier.
    params.ssr.max_J = [];
    params.skip_ssr = false;  % true: skip SSR (err_ssr stays NaN), saves time and memory

    % Apply overrides before derived computations
    if nargin > 0 && ~isempty(overrides)
        params = apply_overrides(params, overrides);
    end

    % === Derived ===
    % n_meas is deterministic from the layout (positions vary, count doesn't)
    rng_state = rng;
    pb = config_pilots(params);
    rng(rng_state);
    params.n_meas = compute_n_meas(pb, params.Nr);

    % Load the target pool only for modes that use it. Explicit targets are
    % self-contained and must not depend on an unrelated pool file.
    if isempty(params.targets)
        targets_file = params.targets_file;
        if ~isfile(targets_file)
            targets_file = fullfile(fileparts(mfilename('fullpath')), targets_file);
        end
        if ~isfile(targets_file)
            error('config_params_HF:targetsFileNotFound', ...
                'Target pool not found: %s', targets_file);
        end
        params.pool = jsondecode(fileread(targets_file));
    else
        params.pool = [];
    end

    % targets_num for pre-allocation
    if ~isempty(params.targets)
        params.targets_num = numel(params.targets);
    elseif ~isempty(params.n_targets)
        params.targets_num = params.n_targets;
    else
        params.targets_num = length(params.selected_targets);
    end

    % Skip AoA estimation when single Rx antenna (no spatial dimension)
    params.skip_aoa = (params.Nr == 1);
end
