function x = strip_padding(y, M, N, padLen, offset, padType)
%STRIP_PADDING Remove CP/ZP padding from time-domain signal.
%   x = strip_padding(y, M, N, padLen, offset, padType)
%
%   Input:  y       - numSamps x 1 time-domain vector (with padding)
%           M       - subcarriers (samples per symbol, without padding)
%           N       - number of symbols
%           padLen  - padding length
%           offset  - sample offset (from derivedParams)
%           padType - 'CP' | 'ZP' | 'RCP' | 'RZP' | 'NONE'
%
%   Output: x       - M x N time-domain matrix (padding stripped)
%
%   Padding modes:
%     CP   - per-symbol: each block is (M+padLen), extract M samples at offset
%     ZP   - per-symbol: each block is (M+padLen), extract first M samples
%     RCP  - frame-level: skip offset samples, take M*N contiguous
%     RZP  - frame-level: take first M*N samples (ZP at end)
%     NONE - no padding: reshape directly

    switch upper(padType)
        case 'CP'
            blocks = reshape(y, M + padLen, N);       % (M+padLen) x N
            x = blocks(offset+1 : offset+M, :);      % M x N
        case 'ZP'
            blocks = reshape(y, M + padLen, N);       % (M+padLen) x N
            x = blocks(1:M, :);                       % M x N
        case 'RCP'
            x = y(offset+1 : offset+M*N);            % M*N x 1
            x = reshape(x, M, N);                    % M x N
        case 'RZP'
            x = y(1 : M*N);                          % M*N x 1
            x = reshape(x, M, N);                    % M x N
        case 'NONE'
            x = reshape(y, M, N);                    % M x N
        otherwise
            error('strip_padding: unknown padType ''%s''', padType);
    end
end
