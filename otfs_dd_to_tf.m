function XTF = otfs_dd_to_tf(XDD, M, N)
%OTFS_DD_TO_TF Delay-Doppler to Time-Frequency domain transform.
%   XTF = otfs_dd_to_tf(XDD, M, N) applies the OTFS ISFFT to convert
%   from the Delay-Doppler (DD) domain to the Time-Frequency (TF) domain.
%
%   Input:
%     XDD - M x N delay-Doppler grid (per Tx antenna)
%     M   - number of delay taps (subcarriers)
%     N   - number of Doppler taps (symbols)
%
%   Output:
%     XTF - M x N time-frequency grid
%
%   Formula: XTF = ifft(fft(XDD, [], 1), [], 2) / sqrt(M/N)

    XTF = ifft(fft(XDD, [], 1), [], 2) / sqrt(M/N);
end
