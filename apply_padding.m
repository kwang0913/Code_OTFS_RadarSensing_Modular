function y = apply_padding(x, padLen, padType)
%APPLY_PADDING Add CP/ZP padding to time-domain OFDM block.
%   y = apply_padding(x, padLen, padType)
%
%   Input:  x       - M x N time-domain matrix (M samples per symbol, N symbols)
%           padLen  - padding length (samples)
%           padType - 'CP' | 'ZP' | 'RCP' | 'RZP' | 'NONE'
%
%   Output: y       - numSamps x 1 time-domain vector with padding applied
%
%   Padding modes:
%     CP   - per-symbol cyclic prefix:  prepend last padLen samples of each symbol
%     ZP   - per-symbol zero padding:   append padLen zeros to each symbol
%     RCP  - frame-level cyclic prefix: prepend last padLen samples of entire frame
%     RZP  - frame-level zero padding:  append padLen zeros to entire frame
%     NONE - no padding

    switch upper(padType)
        case 'CP'
            x = [x(end-padLen+1:end, :); x];   % (M+padLen) x N
            y = x(:);
        case 'ZP'
            x = [x; zeros(padLen, size(x, 2))]; % (M+padLen) x N
            y = x(:);
        case 'RCP'
            y = x(:);                           % M*N x 1
            y = [y(end-padLen+1:end); y];       % M*N+padLen x 1
        case 'RZP'
            y = x(:);                           % M*N x 1
            y = [y; zeros(padLen, 1)];          % M*N+padLen x 1
        case 'NONE'
            y = x(:);                           % M*N x 1
        otherwise
            error('apply_padding: unknown padType ''%s''', padType);
    end
end
