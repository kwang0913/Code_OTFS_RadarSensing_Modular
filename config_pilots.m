function pilot_blocks = config_pilots(params)
%CONFIG_PILOTS Pilot layout: structured blocks (pilot_blocks_def) + n_pilotspertx random pilots per Tx without a block.
%
%   Fields per block:
%     dim    - 'delay', 'doppler' or 'random' (random = explicit (m,n) pairs)
%     m, n   - subcarrier/time-slot indices (array in the spanning dim)
%     tx     - Tx antenna for pilots (1-based)
%     ref_m, ref_n, ref_tx - reference position and antenna
%
%   Derived automatically: offsets, n_pilots.
%
%   Blocks come from params.pilot_blocks_def (e.g. from JSON; none if empty).
%   The m and n fields accept MATLAB range strings such as "2:128" as well as
%   numeric scalars or arrays. Every Tx without a block then gets
%   params.n_pilotspertx random pilots (scalar, or one entry per Tx; 0 = none).

    % ================================================================
    %  PILOT BLOCK DEFINITIONS
    % ================================================================
    pb = struct('dim', {}, 'm', {}, 'n', {}, 'tx', {}, 'ref_m', {}, 'ref_n', {}, 'ref_tx', {});
    if isfield(params, 'pilot_blocks_def') && ~isempty(params.pilot_blocks_def)
        % --- From JSON / override ---
        defs = params.pilot_blocks_def;
        if iscell(defs)
            % jsondecode returns a cell array when struct fields differ
            for i = numel(defs):-1:1
                pb(i) = parse_pilot_block_def(defs{i});
            end
        else
            % struct array (all blocks share the same field set)
            for i = numel(defs):-1:1
                pb(i) = parse_pilot_block_def(defs(i));
            end
        end
    end

    % Random pilots on Tx without a block
    n_pilotspertx = params.n_pilotspertx;
    if isscalar(n_pilotspertx)
        n_pilotspertx = repmat(n_pilotspertx, 1, params.Nt);
    elseif numel(n_pilotspertx) ~= params.Nt
        error('config_pilots:badPilotCount', ...
            'n_pilotspertx must be a scalar or have one entry per Tx (Nt = %d).', params.Nt);
    end
    M = params.M;
    N = params.N;

    % Collect all [m,n] positions used by blocks above
    used = false(M, N);
    b = numel(pb);
    for i = 1:b
        blk_i = pb(i);
        if strcmp(blk_i.dim, 'delay')
            for mi = 1:length(blk_i.m)
                used(blk_i.m(mi), blk_i.n) = true;
            end
        elseif strcmp(blk_i.dim, 'random')
            for p = 1:length(blk_i.m)
                used(blk_i.m(p), blk_i.n(p)) = true;
            end
        else  % doppler
            for ni = 1:length(blk_i.n)
                used(blk_i.m, blk_i.n(ni)) = true;
            end
        end
        used(blk_i.ref_m, blk_i.ref_n) = true;
    end

    [avail_m, avail_n] = find(~used);  % all TF bins (m,n) not occupied by the blocks above
    perm = randperm(length(avail_m));
    avail_m = avail_m(perm);
    avail_n = avail_n(perm);

    remaining = setdiff(1:params.Nt, [pb.tx]);
    if sum(n_pilotspertx(remaining)) > length(avail_m)
        error('config_pilots:notEnoughBins', ...
            'Not enough free TF bins for %d random pilots.', sum(n_pilotspertx(remaining)));
    end
    offset_pos = 0;
    for tx = remaining  % assign random pilots to Tx without a block
        if n_pilotspertx(tx) == 0, continue; end
        sel = offset_pos + (1:n_pilotspertx(tx));
        b = b + 1;
        pb(b).dim    = 'random';
        pb(b).m      = avail_m(sel)';
        pb(b).n      = avail_n(sel)';
        pb(b).tx     = tx;
        pb(b).ref_m  = avail_m(sel(1));
        pb(b).ref_n  = avail_n(sel(1));
        pb(b).ref_tx = tx;
        offset_pos = offset_pos + n_pilotspertx(tx);
    end

    % ================================================================
    %  Derived quantities
    % ================================================================
    for i = 1:length(pb)
        if strcmp(pb(i).dim, 'delay')
            pb(i).offsets = pb(i).m - pb(i).ref_m;
        elseif strcmp(pb(i).dim, 'random')
            pb(i).offsets = 0:(length(pb(i).m) - 1);
        else
            pb(i).offsets = pb(i).n - pb(i).ref_n;
        end
        pb(i).n_pilots = length(pb(i).offsets);
    end

    sites = pilot_sites(pb);
    if sites.n_sites > sites.n_bins  % some TF bin carries two Tx antennas
        [bins, ~, ic] = unique([sites.m sites.n], 'rows', 'stable');
        k = find(accumarray(ic, 1) > 1, 1);
        error('config_pilots:sharedBin', ...
            'TF bin (m,n)=(%d,%d) is used by Tx %s; one Tx per bin.', ...
            bins(k,1), bins(k,2), mat2str(sites.tx(ic == k).'));
    end

    pilot_blocks = pb;
end

% -----------------------------------------------------------------
function blk = parse_pilot_block_def(def)
%PARSE_PILOT_BLOCK_DEF Convert a JSON-decoded struct to a pilot block.
    blk.dim    = char(def.dim);
    blk.m      = parse_range(def.m);
    blk.n      = parse_range(def.n);
    blk.tx     = def.tx;
    blk.ref_m  = parse_range(def.ref_m);
    blk.ref_n  = parse_range(def.ref_n);
    blk.ref_tx = def.ref_tx;
end
