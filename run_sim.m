function results = run_sim(params)
%RUN_SIM Monte Carlo simulation: SNR x trials.
%   results = RUN_SIM(params) runs the SNR x trials simulation and returns
%   a struct with fields: pct_err_coarse, pct_err_ssr, n_detected,
%   SNR_db_list, params.
%   pct_err dimensions: (trials, nSNR, K, 3) — columns: [angle, delay, doppler].
%
%   Each trial is seeded with rng(trial_idx), making all randomness
%   (target selection, pilot placement, path gains, QAM data, noise)
%   deterministic and reproducible.

    use_parallel = params.trials > 1;

    nSNR   = length(params.SNR_db_list);
    K      = params.targets_num;

    % Pre-allocate
    pct_err_coarse = NaN(params.trials, nSNR, K, 3);
    pct_err_ssr    = NaN(params.trials, nSNR, K, 3);
    n_detected     = zeros(params.trials, nSNR);

    for snr_idx = 1:nSNR
        SNR_db = params.SNR_db_list(snr_idx);
        SNR_lin = 10^(SNR_db / 10);
        fprintf('--- SNR = %d dB ---\n', SNR_db);

        temp_pct_coarse  = NaN(params.trials, K, 3);
        temp_pct_ssr     = NaN(params.trials, K, 3);
        temp_n_detected  = zeros(params.trials, 1);


        if use_parallel
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) update_progress());
            params_par = params;
            params_par.show_figures   = false;
            params_par.verbose_output = false;
            parfor (trial_idx = 1:params.trials, 6)
                try
                    rng(trial_idx);
                    chanParams_t = generate_channel_params(params_par);
                    [XTF_t, XDD_t, params_t] = generate_TF(params_par);
                    [YTF_t, YDD_t, rx_OTFS_t] = channel_and_noise(XTF_t, chanParams_t, params_t, SNR_db);
                    n0 = mean(abs(rx_OTFS_t(:)).^2) / (1 + SNR_lin);
                    coarse_t = estimate_coarse_params(YTF_t, YDD_t, XTF_t, XDD_t, params_t, chanParams_t, n0);
                    ssr_t    = estimate_ssr_params(YTF_t, XTF_t, coarse_t, params_t, n0);
                    [pce, ~] = compute_pct_errors(chanParams_t, coarse_t);
                    [pse, ~] = compute_pct_errors(chanParams_t, ssr_t);
                    temp_pct_coarse(trial_idx, :, :) = pce(:, 1:3);
                    temp_pct_ssr(trial_idx, :, :)    = pse(:, 1:3);
                    temp_n_detected(trial_idx) = coarse_t.n_detected;
                catch ME
                    warning('Trial %d failed at %s:%d: %s', trial_idx, ME.stack(1).name, ME.stack(1).line, ME.message);
                end
                send(dq, 0);
            end
        else
            for trial_idx = 1:params.trials
                rng(1);  % fixed seed for single-trial debugging
                chanParams_t = generate_channel_params(params);
                [XTF_t, XDD_t, params_t] = generate_TF(params);
                [YTF_t, YDD_t, rx_OTFS_t] = channel_and_noise(XTF_t, chanParams_t, params_t, SNR_db);
                n0 = mean(abs(rx_OTFS_t(:)).^2) / (1 + SNR_lin);
                coarse_t = estimate_coarse_params(YTF_t, YDD_t, XTF_t, XDD_t, params_t, chanParams_t, n0);
                params_t.chanParams_for_plot = chanParams_t;
                ssr_t    = estimate_ssr_params(YTF_t, XTF_t, coarse_t, params_t, n0);
                [pce, ~] = compute_pct_errors(chanParams_t, coarse_t);
                [pse, ~] = compute_pct_errors(chanParams_t, ssr_t);
                temp_pct_coarse(trial_idx, :, :) = pce(:, 1:3);
                temp_pct_ssr(trial_idx, :, :)    = pse(:, 1:3);
                temp_n_detected(trial_idx) = coarse_t.n_detected;

                last_chanParams = chanParams_t;

                if params.verbose_output
                    fprintf('  Trial %d done\n', trial_idx);
                end
            end
        end
        fprintf('\n');

        n_detected(:, snr_idx)           = temp_n_detected;
        pct_err_coarse(:, snr_idx, :, :) = temp_pct_coarse;
        pct_err_ssr(:, snr_idx, :, :)    = temp_pct_ssr;
        fprintf('\n');
    end

    % Package results
    results.pct_err_coarse = pct_err_coarse;  % (trials, nSNR, K, 3): [angle, delay, doppler]
    results.pct_err_ssr    = pct_err_ssr;
    results.n_detected     = n_detected;      % (trials, nSNR)
    results.SNR_db_list    = params.SNR_db_list;
    results.params         = params;
    if exist('last_chanParams', 'var')
        results.chanParams = last_chanParams;
    end
end

function update_progress()
    persistent count t0
    if isempty(count)
        count = 0;
        t0 = tic;
    end
    count = count + 1;
    fprintf('\rProgress: %d (%.1f s)', count, toc(t0));
end
