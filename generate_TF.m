function [XTF, XDD, params] = generate_TF(params)
%GENERATE_TF Build MIMO TF grid with pilot placement and power boosting.
%   [XTF, XDD, params] = generate_TF(params)
%
%   DD guard bins are random positions (independent of TF pilots),
%   n_private(tx) per Tx antenna.
%   TF pilots are placed by config_pilots and boosted by pilot_power.

    % --- Pilot placement (consumes RNG for random mode) ---
    [params.pilot_blocks, params.n_meas, params.n_pilots_per_tx] = config_pilots(params);
    pilot_blocks = params.pilot_blocks;

    % --- Random QAM data grid ---
    databits = randi([0 1], params.M*log2(params.ModOrder), params.N, params.Nt);
    int_syms = bit2int(reshape(databits, params.M*log2(params.ModOrder), []), log2(params.ModOrder));
    if isfield(params, 'modulation') && strcmpi(params.modulation, 'PSK')
        XDD = reshape(pskmod(int_syms, params.ModOrder), params.M, params.N, params.Nt);
    else
        XDD = reshape(qammod(int_syms, params.ModOrder, UnitAveragePower=true), params.M, params.N, params.Nt);
    end

    % --- DD guard: zero random positions per Tx antenna ---
    n_priv = params.n_private;
    if isscalar(n_priv), n_priv = repmat(n_priv, 1, params.Nt); end
    for tx = 1:params.Nt
        np = n_priv(tx);
        if np == 0, continue; end
        idx = randperm(params.M * params.N, np);
        XDD_tx = XDD(:, :, tx);
        XDD_tx(idx) = 0;
        XDD(:, :, tx) = XDD_tx;
    end

    % ISFFT: DD→TF — fft along dim 1 (delay→freq), ifft along dim 2 (Doppler→time)
    XTF = ifft(fft(XDD, [], 1), [], 2) / sqrt(params.M/params.N);

    % TF pilots: zero non-designated Tx, boost designated Tx by pilot_power
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        other_tx = setdiff(1:params.Nt, blk.tx);
        for p = 1:blk.n_pilots
            XTF(blk.m(p), blk.n(p), other_tx) = 0;
            XTF(blk.m(p), blk.n(p), blk.tx) = params.pilot_power * XTF(blk.m(p), blk.n(p), blk.tx);
        end
    end
end
