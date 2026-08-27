function [Meff, numSamps, offset] = compute_padding_params(M, N, padLen, padType)
%COMPUTE_PADDING_PARAMS Derived padding quantities.
    switch upper(padType)
        case {'CP', 'ZP'}
            Meff = M + padLen;
            numSamps = Meff * N;
        case {'RCP', 'RZP'}
            Meff = M;
            numSamps = M * N + padLen;
        case 'NONE'
            Meff = M;
            numSamps = M * N;
        otherwise
            error('Unknown padType: %s', padType);
    end
    if strcmpi(padType, 'CP') || strcmpi(padType, 'RCP')
        offset = padLen;
    else
        offset = 0;
    end
end
