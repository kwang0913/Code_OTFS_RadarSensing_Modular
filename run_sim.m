function results = run_sim(params)
%RUN_SIM Monte Carlo simulation: SNR x trials.
%   results = RUN_SIM(params) runs the SNR x trials simulation and returns
%   a struct with fields:
%     err_coarse, err_ssr  (trials, nSNR, K, 3) absolute errors of matched targets,
%                          columns [angle, delay, doppler]; NaN = unmatched/failed.
%     rng_seeds            (trials, 1) seed of each trial.
%     chanParams           {trials, nSNR} cells of structs with the true channel; [] if failed.
%     ests_coarse, ests_ssr  {trials, nSNR} cells of structs with the estimates
%                          (AoAs, AoDs, pathDelays, pathDopplers, pathGains; SSR also
%                          residual_nmse; one entry per estimated path, any count);
%                          [] if the trial failed or SSR was skipped.
%     n_detected_coarse, n_detected_ssr  (trials, nSNR) estimated path counts.
%     n_matched_coarse, n_matched_ssr    (trials, nSNR) estimates matched to a
%                          true target (non-NaN rows of err_*).
%     BER                  (trials, nSNR) LMMSE data-detection BER from the SSR estimates
%                          (coarse if SSR is skipped); computed only when trials == 1,
%                          else NaN (run_ber_batch replays saved results).
%     failed               (trials, nSNR) trials that errored in the parallel
%                          path (their errors stay NaN, counts 0).
%     SNR_db_list, params.
%
%   Trial i is seeded with rng(params.rng_seed_start + i - 1, 'threefry'), making
%   all randomness (target selection, pilot placement, path gains, QAM data,
%   noise) deterministic and reproducible (explicit 'threefry' generator, as on
%   parpool workers). Seeds are identical across SNR points (common random
%   numbers). Trials run in parallel (params.n_workers workers) only when
%   trials > 1 and n_workers > 0; otherwise serially, where errors propagate
%   and figures/verbose output are honored. The caller's RNG state is restored
%   on exit.

    % Per-trial seeding switches the generator; restore the caller's state.
    caller_rng = rng;
    restore_rng = onCleanup(@() rng(caller_rng)); %#ok<NASGU>
    use_parallel = params.trials > 1 && params.n_workers > 0;

    nSNR   = length(params.SNR_db_list);
    K      = params.targets_num;
    rng_seeds = (params.rng_seed_start:params.rng_seed_start + params.trials - 1).';
    derivedParams = compute_derived_params(params);  % from params only, once per run

    % Pre-allocate
    err_coarse        = NaN(params.trials, nSNR, K, 3);
    err_ssr           = NaN(params.trials, nSNR, K, 3);
    chanParams        = cell(params.trials, nSNR);
    ests_coarse       = cell(params.trials, nSNR);
    ests_ssr          = cell(params.trials, nSNR);
    n_detected_coarse = zeros(params.trials, nSNR);
    n_detected_ssr    = zeros(params.trials, nSNR);
    n_matched_coarse  = zeros(params.trials, nSNR);
    n_matched_ssr     = zeros(params.trials, nSNR);
    BER               = NaN(params.trials, nSNR);
    failed            = false(params.trials, nSNR);

    for snr_idx = 1:nSNR
        SNR_db = params.SNR_db_list(snr_idx);
        fprintf('--- SNR = %d dB ---\n', SNR_db);

        tr = cell(params.trials, 1);  % one result struct per trial; [] if failed

        if use_parallel
            update_progress(params.trials);  % reset counter, set total
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) update_progress());
            params_par = params;
            params_par.show_figures   = false;
            params_par.verbose_output = false;
            parfor (trial_idx = 1:params.trials, params.n_workers)
                try
                    tr{trial_idx} = run_trial(params_par, derivedParams, SNR_db, rng_seeds(trial_idx));
                catch ME
                    warning('Trial %d failed at %s:%d: %s', trial_idx, ME.stack(1).name, ME.stack(1).line, ME.message);
                end
                send(dq, 0);
            end
            fprintf('\n');
            n_failed = sum(cellfun(@isempty, tr));
            if n_failed > 0
                fprintf('  %d/%d trials failed\n', n_failed, params.trials);
            end
        else
            for trial_idx = 1:params.trials
                tr{trial_idx} = run_trial(params, derivedParams, SNR_db, rng_seeds(trial_idx));

                if params.verbose_output
                    fprintf('  Trial %d done\n', trial_idx);
                end
            end
        end
        fprintf('\n');

        for trial_idx = 1:params.trials
            if isempty(tr{trial_idx})
                failed(trial_idx, snr_idx) = true;
                continue
            end
            r = tr{trial_idx};
            err_coarse(trial_idx, snr_idx, :, :) = r.err_coarse;
            err_ssr(trial_idx, snr_idx, :, :)    = r.err_ssr;
            chanParams{trial_idx, snr_idx}       = r.chanParams;
            ests_coarse{trial_idx, snr_idx}      = r.ests_coarse;
            ests_ssr{trial_idx, snr_idx}         = r.ests_ssr;
            n_detected_coarse(trial_idx, snr_idx) = r.n_detected_coarse;
            n_detected_ssr(trial_idx, snr_idx)    = r.n_detected_ssr;
            n_matched_coarse(trial_idx, snr_idx)  = r.n_matched_coarse;
            n_matched_ssr(trial_idx, snr_idx)     = r.n_matched_ssr;
            BER(trial_idx, snr_idx)               = r.ber;
        end
        if params.trials == 1
            fprintf('BER: %.4e\n', BER(1, snr_idx));
        end
    end

    % Package results
    results.err_coarse        = err_coarse;  % (trials, nSNR, K, 3): [angle, delay, doppler]
    results.err_ssr           = err_ssr;
    results.rng_seeds         = rng_seeds;   % (trials, 1)
    results.chanParams        = chanParams;  % (trials, nSNR) cells of structs
    results.ests_coarse       = ests_coarse;
    results.ests_ssr          = ests_ssr;
    results.n_detected_coarse = n_detected_coarse;  % (trials, nSNR)
    results.n_detected_ssr    = n_detected_ssr;
    results.n_matched_coarse  = n_matched_coarse;
    results.n_matched_ssr     = n_matched_ssr;
    results.BER               = BER;
    results.failed            = failed;      % (trials, nSNR)
    results.SNR_db_list       = params.SNR_db_list;
    results.params            = params;
end

function r = run_trial(params, derivedParams, SNR_db, seed)
%RUN_TRIAL One seeded trial: simulate -> coarse -> SSR -> errors (K x 3 each), estimates, BER (trials == 1).
    sim = simulate_trial(params, derivedParams, SNR_db, seed);
    chanParams = sim.chanParams;  params_t = sim.params_t;  n0 = sim.n0;
    if strcmp(params_t.coarse_method, 'fused')
        coarse = estimate_coarse_params_fused(sim.YTF, sim.YDD, sim.XTF, sim.XDD, params_t, chanParams, n0, derivedParams);
    else
        coarse = estimate_coarse_params(sim.YTF, sim.YDD, sim.XTF, sim.XDD, params_t, chanParams, n0, derivedParams);
    end
    [pce, ~] = compute_errors(chanParams, coarse, params_t);
    r.chanParams        = chanParams;
    r.ests_coarse       = path_ests(coarse);
    r.n_detected_coarse = numel(coarse.AoAs);
    r.n_matched_coarse  = sum(~isnan(pce(:, 2)));
    if params_t.skip_ssr
        pse = NaN(size(pce, 1), 3);
        r.ests_ssr       = [];
        r.n_detected_ssr = 0;
        r.n_matched_ssr  = 0;
    else
        ssr = estimate_ssr_params(sim.YTF, sim.XTF, coarse, params_t, chanParams, n0, derivedParams);
        [pse, ~] = compute_errors(chanParams, ssr, params_t);
        r.ests_ssr       = path_ests(ssr);
        r.n_detected_ssr = numel(ssr.AoAs);
        r.n_matched_ssr  = sum(~isnan(pse(:, 2)));
    end
    r.err_coarse = pce(:, 1:3);
    r.err_ssr    = pse(:, 1:3);
    r.ber = NaN;
    if params.trials == 1
        ests = r.ests_ssr;
        if isempty(ests), ests = r.ests_coarse; end
        r.ber = detect_lmmse(sim.rx_OTFS, sim.XTF, params_t, ests, derivedParams, ...
            sim.databits, sim.guard_mask, 10^(-SNR_db / 10));
    end
end

function ests = path_ests(e)
%PATH_ESTS Keep only the per-path estimate fields (and the SSR residual used by detect_lmmse).
    keep = {'AoAs', 'AoDs', 'pathDelays', 'pathDopplers', 'pathGains', 'residual_nmse'};
    ests = rmfield(e, setdiff(fieldnames(e), keep));
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
