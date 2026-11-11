function [YTF, YDD, rx_OTFS, noise_var] = channel_and_noise_rect(XTF, chanParams, params, derivedParams, SNR_db)
%CHANNEL_AND_NOISE_RECT Full PHY pipeline: Heisenberg (TF→time) → time-domain channel → AWGN → Wigner (time→TF/DD).
%   noise_var is the actual receive-noise variance added by awgn (0 for SNR_db = Inf).

    %% 1. Heisenberg Transform (Mod to time domain)
    tx_useful = reshape(ifft(XTF, [], 1) * sqrt(params.M), params.M*params.N, params.Nt);

    %% 2. Time-Domain Channel
    [~, rx_OTFS] = apply_time_channel(tx_useful, chanParams, params, derivedParams, false);

    % 'measured' mode: awgn scales noise relative to signal power (needed
    % because pilot_power modifies the actual signal level)
    noise_var = 0;
    if SNR_db ~= Inf
        [rx_OTFS, noise_var] = awgn(rx_OTFS, SNR_db, 'measured');
    end

    %% 3. Wigner Transform (Demod to TF/DD domain)
    YTF = zeros(params.M, params.N, params.Nr);
    YDD = zeros(params.M, params.N, params.Nr);
    for rx_idx = 1:params.Nr
        temp = strip_padding(rx_OTFS(1:derivedParams.numSamps, rx_idx), params.M, params.N, ...
            derivedParams.padLen, derivedParams.offset, params.padType);
        YTF(:,:,rx_idx) = fft(temp, [], 1) / sqrt(params.M);
        YDD(:,:,rx_idx) = otfs_tf_to_dd(YTF(:,:,rx_idx), params.M, params.N);
    end
end
