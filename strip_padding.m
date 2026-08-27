function x = strip_padding(y, M, N, padLen, offset, padType)
%STRIP_PADDING Remove CP/ZP padding from time-domain signal.
    switch upper(padType)
        case {'CP', 'ZP'}
            Meff = M + padLen;
            x = zeros(M, N);
            for n = 1:N
                seg = y((n-1)*Meff + offset + (1:M));
                x(:, n) = seg;
            end
        case {'RCP', 'RZP'}
            x = reshape(y(offset + (1:M*N)), M, N);
        case 'NONE'
            x = reshape(y(1:M*N), M, N);
        otherwise
            error('Unknown padType: %s', padType);
    end
end
