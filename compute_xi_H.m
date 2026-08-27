function xi_H = compute_xi_H(nu, k, kappa, l, n0, m0, N, M, mode)
%COMPUTE_XI_H Closed-form pilot response xi_j[n]*H^j[n,m] for fractional Doppler.
%
%   xi_H = compute_xi_H(nu, k, kappa, l, n0, m0, N, M)
%   xi_H = compute_xi_H(nu, k, kappa, l, n0, m0, N, M, 'grid')
%
%   mode = 'paired' (default): nu/k/kappa/l are same shape, output matches.
%   mode = 'grid':   nu/k/kappa are Doppler candidates (col), l is delay candidates (row).
%                     Output is n_v x n_r via implicit expansion.

    if nargin < 9, mode = 'paired'; end

    if strcmp(mode, 'grid')
        nu = nu(:);  k = k(:);  kappa = kappa(:);
        l = l(:)';
    end

    % H^j[n0,m0]
    H = exp(-1j*2*pi * nu .* l / (N*M)) ...
      .* exp( 1j*2*pi * k * n0 / N) ...
      .* exp(-1j*2*pi * m0 * l / M);

    % xi_j[n0]
    xi_num = exp(1j*2*pi * nu .* l / (N*M)) - exp(1j*2*pi * nu / N);
    xi_den = 1 - exp(1j*2*pi * nu / (N*M));
    kappa_phase = exp(1j*2*pi * kappa * n0 / N);

    % zero_mask: when kappa=0 (integer Doppler), xi_den→0; guard against division by zero
    zero_mask = abs(xi_den) < 1e-12;
    xi_den(zero_mask) = 1;                   % placeholder to avoid NaN; overwritten below
    xi = (1/M) * (xi_num ./ xi_den) .* kappa_phase;

    % Integer Doppler fallback: xi simplifies to (M-l)/M when kappa=0
    if any(zero_mask(:))
        fallback = ones(size(xi)) .* ((M - l) / M);
        zm_full = logical(zero_mask .* ones(size(xi)));
        xi(zm_full) = fallback(zm_full);
    end

    xi_H = xi .* H;
end
