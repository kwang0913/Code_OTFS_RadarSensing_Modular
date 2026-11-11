function xi_H = compute_xi_H(nu, delay, n0, m0, N, M, derivedParams, padType)
%COMPUTE_XI_H Same-TF-bin coefficient of a point-sampled rectangular symbol.
%   nu and delay are path arrays that broadcast elementwise.
%   The symbol boundary uses ceil(delay), while the carrier phase retains the
%   physical fractional delay. Doppler uses receive time, including padding.

    first = ceil(delay);
    if strcmpi(padType, 'CP')
        first = max(0, first - derivedParams.padLen);
    elseif strcmpi(padType, 'RCP') && N == 1
        first = max(0, first - derivedParams.padLen);
    end
    count = max(0, M - first);
    omega = 2*pi*nu / (derivedParams.Meff*N);
    denominator = expm1(1j*omega);
    singular = abs(denominator) < 1e-12;
    denominator(singular) = 1;
    overlap = exp(1j*omega.*first) .* expm1(1j*omega.*count) ./ denominator / M;
    if any(singular(:))
        mask = singular & true(size(overlap));
        limit = count / M + zeros(size(overlap));
        overlap(mask) = limit(mask);
    end

    time = derivedParams.offset + n0*derivedParams.Meff;
    xi_H = overlap .* exp(1j*omega*time) .* exp(-1j*2*pi*m0*delay/M);
end
