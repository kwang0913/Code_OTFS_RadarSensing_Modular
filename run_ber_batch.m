function BER = run_ber_batch(result_file, estimate)
%RUN_BER_BATCH BER of every recorded trial of a saved run_sim result, without re-estimation.
%   BER = RUN_BER_BATCH(result_file) regenerates each trial and SNR with simulate_trial
%   (stored params and seed), equalizes with detect_lmmse using the STORED estimates
%   and saves BER (trials x nSNR) to results/ber_<result file name> (ber_coarse_<name> for 'coarse').
%   RUN_BER_BATCH(result_file, estimate): estimate = 'ssr' (default) | 'coarse'.
%   Failed trials, and SSR-skipped runs with 'ssr', give NaN.
%
%   Bash: matlab -batch "run_ber_batch('results/<name>.mat')"

    if nargin < 2, estimate = 'ssr'; end
    if ~any(strcmp(estimate, {'ssr', 'coarse'}))
        error('run_ber_batch:estimate', 'estimate must be ''ssr'' or ''coarse''.');
    end
    results = load(result_file).results;
    params = results.params;
    derivedParams = compute_derived_params(params);
    ests = results.(['ests_' estimate]);
    BER = NaN(size(ests));

    for snr_idx = 1:numel(results.SNR_db_list)
        SNR_db = results.SNR_db_list(snr_idx);
        for trial_idx = 1:params.trials
            if isempty(ests{trial_idx, snr_idx}), continue; end
            sim = simulate_trial(params, derivedParams, SNR_db, results.rng_seeds(trial_idx));
            BER(trial_idx, snr_idx) = detect_lmmse(sim.rx_OTFS, sim.XTF, sim.params_t, ...
                ests{trial_idx, snr_idx}, derivedParams, sim.databits, sim.guard_mask, 10^(-SNR_db / 10));
        end
        fprintf('SNR = %g dB: mean BER %.4e\n', SNR_db, mean(BER(:, snr_idx), 'omitnan'));
    end

    out_dir = 'results';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    [~, name] = fileparts(result_file);
    prefix = 'ber_';  if strcmp(estimate, 'coarse'), prefix = 'ber_coarse_'; end
    out_file = fullfile(out_dir, [prefix name '.mat']);
    SNR_db_list = results.SNR_db_list;  rng_seeds = results.rng_seeds;
    save(out_file, 'BER', 'estimate', 'SNR_db_list', 'rng_seeds');
    fprintf('Saved: %s\n', out_file);
end
