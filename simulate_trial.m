function sim = simulate_trial(params, derivedParams, SNR_db, seed)
%SIMULATE_TRIAL Seeded transmit -> channel -> receive part of one trial (no estimation).
%   sim = SIMULATE_TRIAL(params, derivedParams, SNR_db, seed) seeds rng(seed, 'threefry')
%   and returns a struct with fields chanParams, XTF, XDD (transmitted DD grid incl.
%   pilots), XDD_data, databits, guard_mask, params_t (params with the pilot layout),
%   YTF, YDD, rx_OTFS and n0 (noise variance estimated from rx_OTFS and the SNR).
%   run_sim and run_ber_batch share it, so a stored seed regenerates the same trial.
    rng(seed, 'threefry');  % generator of parpool workers: serial and parfor trials give identical streams
    sim.chanParams = generate_channel_params(params);
    [sim.XTF, sim.XDD, sim.params_t, sim.databits, sim.guard_mask, sim.XDD_data] = generate_TF(params);
    [sim.YTF, sim.YDD, sim.rx_OTFS] = channel_and_noise_rect(sim.XTF, sim.chanParams, sim.params_t, derivedParams, SNR_db);
    sim.n0 = estimate_noise_variance(sim.rx_OTFS, 10^(-SNR_db / 10), struct());  % P_rx/(1+SNR)
end
