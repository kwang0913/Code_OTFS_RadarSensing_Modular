function y = apply_padding(x, padLen, padType)
%APPLY_PADDING Add CP/ZP padding to time-domain OTFS signal.
%   x: M x N matrix (per-symbol columns)
%   y: numSamps x 1 column vector with padding applied
    switch upper(padType)
        case 'CP'
            x = [x(end-padLen+1:end, :); x];   % (M+padLen) x N
            y = x(:);
        case 'ZP'
            x = [x; zeros(padLen, size(x, 2))]; % (M+padLen) x N
            y = x(:);
        case 'RCP'
            y = x(:);
            y = [y(end-padLen+1:end); y];
        case 'RZP'
            y = [x(:); zeros(padLen, 1)];
        case 'NONE'
            y = x(:);
        otherwise
            error('Unknown padType: %s', padType);
    end
end
