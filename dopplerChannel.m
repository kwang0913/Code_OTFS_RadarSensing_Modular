function y = dopplerChannel(x, fs, chanParams)
%DOPPLERCHANNEL Time-domain multipath channel: delay + Doppler shift + gain per path.
%   y = dopplerChannel(x, fs, chanParams)
%   Applies delay first, then frequency offset, then complex gain for each path.
    numPaths = length(chanParams.pathDelays);
    maxPathDelay = max(chanParams.pathDelays);
    txOutSize = length(x);

    y = zeros(txOutSize+maxPathDelay,1);

    for k = 1:numPaths
        pathOut = zeros(txOutSize+maxPathDelay,1);
        pathOut(1+chanParams.pathDelays(k):chanParams.pathDelays(k)+txOutSize) = x;
        pathOut = frequencyOffset(pathOut,fs,chanParams.pathDopplerFreqs(k));
        pathOut = pathOut * chanParams.pathGains(k);

        y = y + pathOut;
    end
end