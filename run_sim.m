function results = run_sim(params)
%RUN_SIM Monte Carlo simulation: SNR x trials.
%   results = RUN_SIM(params) runs the SNR x trials simulation and returns
%   a struct with fields: pct_err_coarse, pct_err_ssr, n_detected, failed,
%   SNR_db_list, params, and chanParams (serial runs: last trial's ground truth).
%   pct_err dimensions: (trials, nSNR, K, 3) — columns: [angle, delay, doppler].
%   n_detected and failed have dimensions (trials, nSNR); failed marks trials
%   that errored in the parallel path (their errors stay NaN).
%
%   Trial t is seeded with rng(params.seed + t), making all randomness
%   (target selection, pilot placement, path gains, QAM data, noise)
%   deterministic and reproducible (explicit 'threefry' generator, as on parpool workers). Seeds are identical across SNR points
%   (common random numbers). Trials run in parallel (params.n_workers
%   workers) only when trials > 1 and n_workers > 0; otherwise serially,
%   where errors propagate and figures/verbose output are honored.

    use_parallel = params.trials > 1 && params.n_workers > 0;

    nSNR   = length(params.SNR_db_list);
    K      = params.targets_num;

    % Pre-allocate
    pct_err_coarse = NaN(params.trials, nSNR, K, 3);
    pct_err_ssr    = NaN(params.trials, nSNR, K, 3);
    n_detected     = zeros(params.trials, nSNR);
    failed         = false(params.trials, nSNR);

    for snr_idx = 1:nSNR
        SNR_db = params.SNR_db_list(snr_idx);
        fprintf('--- SNR = %d dB ---\n', SNR_db);

        temp_pct_coarse  = NaN(params.trials, K, 3);
        temp_pct_ssr     = NaN(params.trials, K, 3);
        temp_n_detected  = zeros(params.trials, 1);
        temp_failed      = false(params.trials, 1);

        if use_parallel
            update_progress(params.trials);  % reset counter, set total
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) update_progress());
            params_par = params;
            params_par.show_figures   = false;
            params_par.verbose_output = false;
            parfor (trial_idx = 1:params.trials, params.n_workers)
                try
                    [pce, pse, n_det] = run_trial(params_par, SNR_db, params.seed + trial_idx);
                    temp_pct_coarse(trial_idx, :, :) = pce;
                    temp_pct_ssr(trial_idx, :, :)    = pse;
                    temp_n_detected(trial_idx) = n_det;
                catch ME
                    temp_failed(trial_idx) = true;
                    warning('Trial %d failed at %s:%d: %s', trial_idx, ME.stack(1).name, ME.stack(1).line, ME.message);
                end
                send(dq, 0);
            end
            fprintf('\n');
            if any(temp_failed)
                fprintf('  %d/%d trials failed\n', sum(temp_failed), params.trials);
            end
        else
            for trial_idx = 1:params.trials
                [pce, pse, n_det, last_chanParams] = run_trial(params, SNR_db, params.seed + trial_idx);
                temp_pct_coarse(trial_idx, :, :) = pce;
                temp_pct_ssr(trial_idx, :, :)    = pse;
                temp_n_detected(trial_idx) = n_det;

                if params.verbose_output
                    fprintf('  Trial %d done\n', trial_idx);
                end
            end
        end
        fprintf('\n');

        n_detected(:, snr_idx)           = temp_n_detected;
        failed(:, snr_idx)               = temp_failed;
        pct_err_coarse(:, snr_idx, :, :) = temp_pct_coarse;
        pct_err_ssr(:, snr_idx, :, :)    = temp_pct_ssr;
    end

    % Package results
    results.pct_err_coarse = pct_err_coarse;  % (trials, nSNR, K, 3): [angle, delay, doppler]
    results.pct_err_ssr    = pct_err_ssr;
    results.n_detected     = n_detected;      % (trials, nSNR)
    results.failed         = failed;          % (trials, nSNR)
    results.SNR_db_list    = params.SNR_db_list;
    results.params         = params;
    if exist('last_chanParams', 'var')
        results.chanParams = last_chanParams;
    end
end

function [pce, pse, n_det, chanParams] = run_trial(params, SNR_db, seed)
%RUN_TRIAL One seeded trial: channel -> OTFS -> coarse -> SSR -> % errors (K x 3 each).
    rng(seed, 'threefry');  % generator of parpool workers: serial and parfor trials give identical streams
    chanParams = generate_channel_params(params);
    [XTF, XDD, params_t] = generate_TF(params);
    [YTF, YDD, rx_OTFS] = channel_and_noise(XTF, chanParams, params_t, SNR_db);
    n0 = mean(abs(rx_OTFS(:)).^2) / (1 + 10^(SNR_db / 10));
    coarse = estimate_coarse_params(YTF, YDD, XTF, XDD, params_t, chanParams, n0);
    ssr    = estimate_ssr_params(YTF, XTF, coarse, params_t, chanParams, n0);
    [pce, ~] = compute_pct_errors(chanParams, coarse);
    [pse, ~] = compute_pct_errors(chanParams, ssr);
    pce = pce(:, 1:3);
    pse = pse(:, 1:3);
    n_det = coarse.n_detected;
end

function update_progress(total)
%UPDATE_PROGRESS update_progress(total) resets the counter; update_progress() increments.
    persistent count n_total t0
    if nargin
        count = 0;  n_total = total;  t0 = tic;
        return
    end
    count = count + 1;
    fprintf('\rProgress: %d/%d (%.1f s)', count, n_total, toc(t0));
end
