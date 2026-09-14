function y = dopplerChannel(x, fs, chanParams, params)
%DOPPLERCHANNEL Rectangular OTFS waveform with real-valued delay taps.
%   Fractional delays evaluate each OFDM symbol at its fractional sample
%   times using subcarriers 0:M-1. Padding preserves the symbol boundaries;
%   interpolation never mixes neighboring symbols across a discontinuity.
    numPaths = length(chanParams.pathDelays);
    maxPathDelay = ceil(max(chanParams.pathDelays));
    txOutSize = length(x);

    y = zeros(txOutSize+maxPathDelay,1);

    if any(chanParams.pathDelays ~= floor(chanParams.pathDelays))
        [~, ~, offset] = compute_padding_params(params.M, params.N, ...
            chanParams.padLen, params.padType);
        symbols = strip_padding(x, params.M, params.N, ...
            chanParams.padLen, offset, params.padType);
        spectrum = fft(symbols, [], 1);
    end

    for k = 1:numPaths
        pathOut = zeros(txOutSize+maxPathDelay,1);
        delay = chanParams.pathDelays(k);
        first = ceil(delay);
        if delay == first
            delayed_samples = x;
        else
            % At output sample t, the source time is t-delay.
            % Its containing source interval begins at t-ceil(delay).
            phase = exp(1j*2*pi*(0:params.M-1).'*(first-delay)/params.M);
            shifted_symbols = ifft(spectrum .* phase, [], 1);
            delayed_samples = apply_padding(shifted_symbols, ...
                chanParams.padLen, params.padType);
        end
        pathOut(first+(1:txOutSize)) = delayed_samples;
        pathOut = frequencyOffset(pathOut,fs,chanParams.pathDopplerFreqs(k));
        pathOut = pathOut * chanParams.pathGains(k);

        y = y + pathOut;
    end
end