function [XTF, XDD, params, databits, guard_mask, XDD_data] = generate_TF(params)
%GENERATE_TF Build MIMO TF grid with pilot placement and power boosting.
%   [XTF, XDD, params, databits, guard_mask, XDD_data] = generate_TF(params)
%
%   Pilots (config_pilots): the DD guard zeroes one random bin per unique
%   pilot/reference TF bin on every Tx antenna; all pilot/reference TF sites
%   are cleared on all Tx and filled with constant-modulus pilots of
%   amplitude pilot_power. params.n_meas is the number of measurement rows.
%   XDD is the DD grid actually transmitted (pilot bins included); XDD_data is
%   the data-only DD grid (guard bins zeroed, no pilots).
%   databits:   M*log2(ModOrder) x N x Nt — original data bits per Tx
%   guard_mask: M x N x Nt logical — true at DD positions zeroed as guards

    if ~any(strcmpi(params.pilot_type, {'random', 'deterministic'}))
        error('generate_TF:UnsupportedPilotType', 'pilot_type must be ''random'' or ''deterministic''.');
    end

    % --- Pilot placement (consumes RNG) ---
    params.pilot_blocks = config_pilots(params);
    pilot_blocks = params.pilot_blocks;
    % Guarded DD symbols per Tx = TF bins overwritten by references/pilots.
    n_bins = pilot_sites(pilot_blocks).n_bins;
    params.n_meas = compute_n_meas(pilot_blocks, params.Nr);

    % --- Random QAM data grid ---
    databits = randi([0 1], params.M*log2(params.ModOrder), params.N, params.Nt);
    int_syms = bit2int(reshape(databits, params.M*log2(params.ModOrder), []), log2(params.ModOrder));
    if isfield(params, 'modulation') && strcmpi(params.modulation, 'PSK')
        XDD = reshape(pskmod(int_syms, params.ModOrder), params.M, params.N, params.Nt);
    else
        XDD = reshape(qammod(int_syms, params.ModOrder, UnitAveragePower=true), params.M, params.N, params.Nt);
    end

    % --- DD guard: zero random positions per Tx antenna ---
    guard_mask = false(params.M, params.N, params.Nt);
    for tx_idx = 1:params.Nt
        temp = XDD(:,:,tx_idx);
        temp = temp(:);
        idx = randperm(params.M*params.N, n_bins);
        temp(idx) = 0;
        XDD(:,:,tx_idx) = reshape(temp, params.M, params.N);

        mask_lin = false(params.M*params.N, 1);
        mask_lin(idx) = true;
        guard_mask(:,:,tx_idx) = reshape(mask_lin, params.M, params.N);
    end
    XDD_data = XDD;

    % ISFFT: DD→TF — fft along dim 1 (delay→freq), ifft along dim 2 (Doppler→time)
    XTF = ifft(fft(XDD, [], 1), [], 2) / sqrt(params.M/params.N);

    % Deduplicate shared references: multiple pilot blocks may share one ref symbol
    ref_keys = zeros(length(pilot_blocks), 3);
    for b = 1:length(pilot_blocks)
        ref_keys(b,:) = [pilot_blocks(b).ref_m, pilot_blocks(b).ref_n, pilot_blocks(b).ref_tx];
    end
    [unique_refs, ~, ref_map] = unique(ref_keys, 'rows');

    ref_vals = zeros(size(unique_refs, 1), 1);
    for r = 1:size(unique_refs, 1)
        if strcmpi(params.pilot_type, 'random')
            ref_vals(r) = exp(1j * (pi/4 + randi([0 3]) * pi/2));
        else
            ref_vals(r) = 1;
        end
    end

    % Pass 1: Zero guard zone — clear all pilot/reference TF positions before inserting pilots
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        if strcmp(blk.dim, 'random')
            for p = 1:length(blk.m)
                XTF(blk.m(p), blk.n(p), 1:params.Nt) = 0;
            end
        else
            XTF(blk.m, blk.n, 1:params.Nt) = 0;
        end
        XTF(blk.ref_m, blk.ref_n, 1:params.Nt) = 0;
    end

    % Pass 2: Place pilots and references
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);

        if strcmpi(params.pilot_type, 'random')
            pilot_seq = exp(1j * (pi/4 + randi([0 3], 1, blk.n_pilots) * pi/2));
        else
            pilot_seq = ones(1, blk.n_pilots);
        end

        if strcmp(blk.dim, 'random')
            for p = 1:blk.n_pilots
                XTF(blk.m(p), blk.n(p), blk.tx) = params.pilot_power * pilot_seq(p);
            end
        elseif strcmp(blk.dim, 'delay')
            XTF(blk.m, blk.n, blk.tx) = params.pilot_power * reshape(pilot_seq, length(blk.m), []);
        else
            XTF(blk.m, blk.n, blk.tx) = params.pilot_power * pilot_seq(:).';
        end

        XTF(blk.ref_m, blk.ref_n, blk.ref_tx) = params.pilot_power * ref_vals(ref_map(b));
    end

    % DD grid actually transmitted (pilot bins included); reference for coarse xcorr
    XDD = fft(ifft(XTF, [], 1), [], 2) / sqrt(params.N/params.M);
end
