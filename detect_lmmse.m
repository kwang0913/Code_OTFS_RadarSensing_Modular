function [ber, XDD_hat, info] = detect_lmmse(rx_OTFS, XTF, params, final_ests, derivedParams, databits, guard_mask, n0, noise_var)
%DETECT_LMMSE Matrix-free LMMSE on the transmitted affine DD-data/TF-pilot frame.
%   Fractional rectangular-symbol delays and their exact adjoint are shared
%   with the PHY. Unknowns are only unguarded DD symbols. TF pilot locations
%   are overwritten, not superimposed, so the data operator masks them out.
%   XDD_hat contains soft DD estimates; BER uses constellation hard decisions.
%   info records LSQR flag (0 = converged), relative residual, and iterations.
%   N0 is the nominal inverse SNR; unless NOISE_VAR (actual receive-noise
%   variance, e.g. from channel_and_noise_rect) is given, the noise level is
%   estimate_noise_variance(rx_OTFS, n0, final_ests), i.e. from the full received
%   frame, the same power reference as simulate_trial's n0.
%   final_ests may carry per-Tx gains (gain_per_tx, see apply_time_channel).
    M = params.M; N = params.N; Nt = params.Nt; Nr = params.Nr;
    bps = log2(params.ModOrder);
    sites = pilot_sites(params.pilot_blocks);
    pilot_mask = false(M,N);
    pilot_mask(sub2ind([M,N],sites.m,sites.n)) = true;
    tf_mask = repmat(pilot_mask,1,1,Nt);
    pilots = zeros(M,N,Nt);
    pilots(tf_mask) = XTF(tf_mask);
    pilot_time = reshape(ifft(pilots,[],1)*sqrt(M),M*N,Nt);
    pilot_rx = apply_time_channel(pilot_time,final_ests,params,derivedParams,false);
    useful_rx = zeros(M*N,Nr);
    for r = 1:Nr
        block = strip_padding(rx_OTFS(1:derivedParams.numSamps,r),M,N, ...
            derivedParams.padLen,derivedParams.offset,params.padType);
        useful_rx(:,r) = block(:);
    end
    y = useful_rx-pilot_rx;
    data_indices = find(~guard_mask);
    n_data = numel(data_indices);
    if nargin >= 9 && ~isempty(noise_var)
        noise = noise_var;
    else
        noise = estimate_noise_variance(rx_OTFS,n0,final_ests);
    end
    root_noise = sqrt(noise);
    operator = @apply_augmented;
    [symbols,flag,relres,iterations] = lsqr(operator,[y(:);zeros(n_data,1)],1e-6,500);
    info = struct('flag',flag,'relres',relres,'iterations',iterations);
    if flag ~= 0
        warning('detect_lmmse:Convergence', ...
            ['LMMSE LSQR did not converge: flag=%d, iterations=%d, relative residual=%.3g. ' ...
             'The returned BER uses the available iterate, not a converged solution.'], ...
            flag,iterations,relres);
    end
    XDD_hat = zeros(M,N,Nt);
    XDD_hat(data_indices) = symbols;
    total_errors = 0; total_bits = 0;
    for t = 1:Nt
        if isfield(params,'modulation') && strcmpi(params.modulation,'PSK')
            detected = pskdemod(XDD_hat(:,:,t),params.ModOrder);
        else
            detected = qamdemod(XDD_hat(:,:,t),params.ModOrder,UnitAveragePower=true);
        end
        detected_bits = reshape(int2bit(detected(:),bps),bps,M*N);
        original_bits = reshape(databits(:,:,t),bps,M*N);
        mask = reshape(~guard_mask(:,:,t),1,M*N);
        total_errors = total_errors + sum(detected_bits(:,mask) ~= original_bits(:,mask),'all');
        total_bits = total_bits + bps*sum(mask);
    end
    ber = total_errors/total_bits;

    function out = apply_augmented(in,transp)
        if strcmp(transp,'notransp')
            dd = zeros(M,N,Nt);
            dd(data_indices) = in;
            tf = otfs_dd_to_tf(dd,M,N);
            tf(tf_mask) = 0;
            tx = reshape(ifft(tf,[],1)*sqrt(M),M*N,Nt);
            rx = apply_time_channel(tx,final_ests,params,derivedParams,false);
            out = [rx(:);root_noise*in];
        else
            rx = reshape(in(1:M*N*Nr),M*N,Nr);
            tx = apply_time_channel(rx,final_ests,params,derivedParams,true);
            tf = fft(reshape(tx,M,N,Nt),[],1)/sqrt(M);
            tf(tf_mask) = 0;
            dd = otfs_tf_to_dd(tf,M,N);
            out = dd(data_indices)+root_noise*in(M*N*Nr+1:end);
        end
    end
end
