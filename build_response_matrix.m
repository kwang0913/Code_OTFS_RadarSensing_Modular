function Phi = build_response_matrix(ests, params)
%BUILD_RESPONSE_MATRIX Channel dictionary using full xi_j[n]*H^j[n,m] with fractional Doppler.
%   Monostatic DFRC: Tx steering uses ests.AoAs.

    pilot_blocks = params.pilot_blocks;
    Nr = params.Nr;
    M = params.M;
    N = params.N;
    K = length(ests.AoAs);

    n_meas = compute_n_meas(pilot_blocks, Nr);

    Phi = zeros(n_meas, K);

    % AoA steering: Nr x K
    a_steer = exp(1j*2*pi*params.ddlr*(0:Nr-1).'*sind(ests.AoAs(:)'));

    nu_arr    = ests.pathDopplers(:)';
    k_arr     = floor(nu_arr);       % integer Doppler index
    kappa_arr = nu_arr - k_arr;      % fractional Doppler remainder in [0,1)
    l_arr     = ests.pathDelays(:)';

    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);

        d_pilot = exp(-1j*2*pi*params.ddlt*(blk.tx - 1)*sind(ests.AoAs(:)'));

        for p = 1:blk.n_pilots
            m0 = blk.m(p) - 1;
            n0 = blk.n(p) - 1;
            xi_H_p = compute_xi_H(nu_arr, k_arr, kappa_arr, l_arr, n0, m0, N, M);
            Phi(idx+(1:Nr), :) = a_steer .* d_pilot .* xi_H_p;
            idx = idx + Nr;
        end
    end
end
