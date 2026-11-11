function Phi_est = build_response_matrix(final_ests, params, derivedParams)
%BUILD_RESPONSE_MATRIX Same-bin pilot dictionary with physical receive-time gains.
%   Rows follow build_measurement_vector: Nr rows per unique pilot/reference
%   TF site (pilot_sites order).

    sites = pilot_sites(params.pilot_blocks);
    Nr = params.Nr;
    M = params.M;
    N = params.N;
    K = length(final_ests.AoAs);

    Phi_est = zeros(sites.n_sites*Nr, K);

    % AoA steering: Nr x K
    a_steer = exp(1j*2*pi*params.ddlr*(0:Nr-1).'*sind(final_ests.AoAs(:)'));

    nu_arr    = final_ests.pathDopplers(:)';
    % Fractional delay and Doppler remain unquantized.
    l_arr     = final_ests.pathDelays(:)';

    for s = 1:sites.n_sites
        % AoD steering for the transmit antenna carrying this site
        d_site = exp(-1j*2*pi*params.ddlt*(sites.tx(s) - 1)*sind(final_ests.AoDs(:)'));
        xi_H = compute_xi_H(nu_arr, l_arr, sites.n(s)-1, sites.m(s)-1, N, M, derivedParams, params.padType);
        Phi_est((s-1)*Nr+(1:Nr), :) = a_steer .* d_site .* xi_H;
    end
end
