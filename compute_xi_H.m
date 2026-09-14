function xi_H = compute_xi_H(nu, k, kappa, l, n0, m0, N, M, mode)
%COMPUTE_XI_H Closed-form pilot response xi_j*H^j[n,m] for real delay and Doppler.
%
%   xi_H = compute_xi_H(nu, k, kappa, l, n0, m0, N, M)
%   xi_H = compute_xi_H(nu, k, kappa, l, n0, m0, N, M, 'grid')
%
%   mode = 'paired' (default): nu/k/kappa/l are same shape, output matches.
%   mode = 'grid':   nu/k/kappa are Doppler candidates (col), l is delay candidates (row).
%                     Output is n_v x n_r via implicit expansion.
%   k and kappa are retained for callers; the response uses nu directly.

    if nargin < 9, mode = 'paired'; end

    if strcmp(mode, 'grid')
        nu = nu(:);
        l = l(:)';
    end

    % H^j[n0,m0]: the actual delay remains in both continuous phases.
    H = exp(-1j*2*pi * nu .* l / (N*M)) ...
      .* exp( 1j*2*pi * nu * n0 / N) ...
      .* exp(-1j*2*pi * m0 * l / M);

    % xi_j = (1/M) sum_{q=ceil(l)}^{M-1} exp(1j*2*pi*nu*q/(N*M)).
    first_q = ceil(l);
    n_terms = M - first_q;

    % Only the integer-sample sum is periodic in nu with period N*M.
    % Center its phase to keep expm1 accurate near every zero denominator.
    nu_sum = rem(nu, N*M);
    nu_sum = nu_sum - round(nu_sum / (N*M)) * (N*M);
    omega = 2*pi * nu_sum / (N*M);
    xi_den = expm1(1j * omega);
    zero_mask = (omega == 0);
    xi_den(zero_mask) = 1;
    xi = exp(1j * omega .* first_q) ...
       .* expm1(1j * omega .* n_terms) ./ xi_den;

    % Zero Doppler limit: every summand equals one.
    xi = xi + zero_mask .* n_terms;

    xi_H = (xi / M) .* H;
end
