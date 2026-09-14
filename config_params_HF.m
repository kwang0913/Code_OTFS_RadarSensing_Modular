function params = config_params_HF(overrides)
%CONFIG_PARAMS_HF High-frequency (24.25 GHz) monostatic DFRC simulation parameters.
%   Monostatic: AoA = AoD always; round-trip delay taps = 2*range/c * fsamp.

    params.config_name = 'hf';

    % --- Control ---
    params.show_figures = true;
    params.verbose_output = true;

    params.pilot_power = 2;
    params.pilot_type = 'random';  % 'random' | 'deterministic'
    params.padType = 'RCP';         % 'CP' | 'RCP' | 'ZP' | 'RZP' | 'NONE'

    % --- System ---
    params.Nt = 2;
    params.Nr = 2;
    params.fc = 24.25e9;  % Hz
    params.M  = 128;      % delay bins (subcarriers)
    params.N  = 64;       % Doppler bins (OTFS symbols)
    params.df = 120e3;    % subcarrier spacing (Hz)
    params.c  = 299792458; % exact speed of light in vacuum (m/s)

    % --- Antenna Spacing (in wavelengths) ---
    params.ddlt = 0.5;
    params.ddlr = 0.5;

    params.ModOrder  = 4;
    params.modulation = 'QAM';  % 'QAM' | 'PSK'

    % --- Simulation Loop ---
    params.SNR_db_list = 20;          % single SNR for debugging; use e.g. -10:5:30 for sweeps
    params.trials = 1;

    % --- Pilot layout ---
    params.n_private = 16;
    % params.n_private = [10 10 10 10];    % pilots per Tx antenna (vector, one per Tx)
    params.pilot_positions = {};     % deterministic mode: cell array of Nx2 [m,n] per Tx

    % --- Target Selection (3 modes, checked in order) ---
    % Mode 3: deterministic — explicit target definitions (highest priority)
    % params.targets = [struct('angle',20,'delay',5,'doppler',-4.2), ...
    %                   struct('angle',35,'delay',8,'doppler',5.4), ...
    %                   struct('angle',-15,'delay',12,'doppler',3.1), ...
    %                   struct('angle',45,'delay',3,'doppler',-6.0)];
    params.targets = [];

    % Mode 2: random — pick n_targets from pool, optionally override with ranges
    params.n_targets = [1];
    params.aoa_range = [-60 60];             % angle list [deg], e.g. -60:5:60 (empty = use pool values)
    params.delay_range = [0 12];      % 0.1-tap grid; scalar max means [1,max], or [lo,hi]; [] uses pool
    params.doppler_range = [3 6];     % magnitude on 0.1-tap grid, random sign; [] uses pool

    % Mode 1: pool selection (default if modes 3 & 2 inactive)
    params.selected_targets = [5];
    params.targets_file = 'targets/pool_5.json';

    % Post-processing (applied on top of any mode)
    params.kappa = [];                 % fixed fractional Doppler
    params.epsilon = [];               % fractional delay override; scalar or one value per target in [0,1)
    params.gain = 1;                   % fixed path gain; scalar or one value per target

    % --- SSR stopping ---
    % [] uses targets_num; a positive integer overrides that limit.
    % Residual-based stopping can stop earlier.
    params.ssr.max_J = [];

    % --- SSR Grid: search range (+/- from coarse estimate) ---
    params.ssr.range.a = 5;   % angle
    params.ssr.range.r = 1;   % delay
    params.ssr.range.v = 1;   % Doppler

    % --- SSR Grid: step length ---
    params.ssr.length.a = 1;
    params.ssr.length.r = 0.1;
    params.ssr.length.v = 0.1;

    % Apply overrides before derived computations
    if nargin > 0 && ~isempty(overrides)
        params = apply_overrides(params, overrides);
    end

    % === Derived quantities ===

    % Sample rate and resolution (padding-independent)
    params.fsamp       = params.M * params.df;                  % sampling rate (Hz)
    params.lambda      = params.c / params.fc;                  % wavelength (m)
    params.range_res   = params.c / (2 * params.M * params.df); % monostatic range resolution (m)
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

    % n_meas is deterministic from n_private (positions vary, count doesn't)
    n_priv = params.n_private;
    if isscalar(n_priv), n_priv = repmat(n_priv, 1, params.Nt); end
    params.n_meas = sum(n_priv) * params.Nr;

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
