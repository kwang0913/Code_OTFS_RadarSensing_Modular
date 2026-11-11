function n0_actual = estimate_noise_variance(rx_OTFS, inv_snr, ests)
%ESTIMATE_NOISE_VARIANCE Effective noise variance from thermal + CSI error.
%   n0_actual = estimate_noise_variance(rx_OTFS, inv_snr, ests)
%   inv_snr = 10^(-SNR_db/10) (linear inverse SNR, not a noise variance). The
%   thermal part is P_rx/(1+SNR); ests.residual_nmse, if present, adds a CSI-error floor.

    P_rx = mean(abs(rx_OTFS(:)).^2);
    n0_thermal = P_rx * inv_snr / (1 + inv_snr);  % P_rx = P_sig + P_noise => P_noise = P_rx/(1+SNR)
    if isfield(ests, 'residual_nmse')
        n0_csi = P_rx * ests.residual_nmse;
    else
        n0_csi = 0;
    end
    n0_actual = max(n0_thermal, n0_csi);  % worst-case: use whichever dominates (thermal or CSI error)
end
