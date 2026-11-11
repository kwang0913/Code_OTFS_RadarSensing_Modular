function [Meff, numSamps, offset] = compute_padding_params(M, N, padLen, padType)
%COMPUTE_PADDING_PARAMS Compute Meff, numSamps, and offset from padding config.
%   [Meff, numSamps, offset] = compute_padding_params(M, N, padLen, padType)

    if strcmpi(padType, 'ZP') || strcmpi(padType, 'CP')
        Meff = M + padLen;  % ZP/CP extend each symbol by padLen samples
    else
        Meff = M;           % RCP/RZP/NONE: no per-symbol padding
    end

    if strcmpi(padType, 'RCP') || strcmpi(padType, 'CP')
        offset = padLen;    % receiver skips the prefix (per symbol for CP, once per frame for RCP)
    else
        offset = 0;
    end

    if strcmpi(padType, 'CP') || strcmpi(padType, 'ZP') || strcmpi(padType, 'NONE')
        numSamps = Meff * N;          % CP/ZP: padding inside each symbol; NONE: none
    else  % RCP, RZP
        numSamps = Meff * N + padLen;  % one frame-level prefix (RCP) or suffix (RZP)
    end
end
