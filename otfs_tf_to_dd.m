function YDD = otfs_tf_to_dd(YTF, M, N)
%OTFS_TF_TO_DD Time-Frequency to Delay-Doppler domain transform.
%   YDD = otfs_tf_to_dd(YTF, M, N) applies the OTFS Wigner transform to
%   convert from the Time-Frequency (TF) domain to the Delay-Doppler (DD)
%   domain.
%
%   Input:
%     YTF - M x N time-frequency grid
%     M   - number of delay taps (subcarriers)
%     N   - number of Doppler taps (symbols)
%
%   Output:
%     YDD - M x N delay-Doppler grid (per Rx antenna)
%
%   Formula: YDD = fft(ifft(YTF, [], 1), [], 2) / sqrt(N/M)

    YDD = fft(ifft(YTF, [], 1), [], 2) / sqrt(N/M);
end
