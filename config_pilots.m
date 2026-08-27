function [pilot_blocks, n_meas_total, n_pilots_per_tx] = config_pilots(params)
%CONFIG_PILOTS Pilot layout for monostatic DFRC.
%   Each Tx antenna gets its own pilot block. Positions are either random
%   or deterministic (explicit [m, n] pairs).
%
%   Fields per block:
%     dim      - always 'random' (position-explicit indexing)
%     m        - subcarrier indices (1-based)
%     n        - symbol indices (1-based)
%     tx       - Tx antenna (1-based)
%     n_pilots - number of pilot measurements in this block
%
%   Config fields:
%     params.n_private    - vector [n1, n2, ..., nNt]: pilots per Tx antenna
%     params.pilot_type   - 'random' | 'deterministic'
%     params.pilot_positions - cell array of Nx2 [m, n] matrices (deterministic only)
%     params.pilot_blocks_def - if non-empty, used directly (JSON override)

    % ================================================================
    %  PILOT BLOCK DEFINITIONS
    % ================================================================
    if isfield(params, 'pilot_blocks_def') && ~isempty(params.pilot_blocks_def)
        % --- From JSON / override ---
        defs = params.pilot_blocks_def;
        if iscell(defs)
            for i = numel(defs):-1:1
                pb(i) = parse_pilot_block_def(defs{i});
            end
        else
            for i = numel(defs):-1:1
                pb(i) = parse_pilot_block_def(defs(i));
            end
        end
    else
        % --- Build from n_private + pilot_type ---
        M  = params.M;
        N  = params.N;
        Nt = params.Nt;
        n_private = params.n_private;

        % Allow scalar shorthand: expand to vector
        if isscalar(n_private)
            n_private = repmat(n_private, 1, Nt);
        end

        used = false(M, N);

        b = 0;
        for tx = 1:Nt
            np = n_private(tx);
            if np == 0, continue; end

            if strcmp(params.pilot_type, 'random')
                [avail_m, avail_n] = find(~used);
                perm = randperm(length(avail_m));
                if np > length(avail_m)
                    error('config_pilots: not enough free DD grid positions for %d random pilots on Tx %d', np, tx);
                end
                sel = perm(1:np);
                pm = avail_m(sel)';
                pn = avail_n(sel)';
            else
                % Deterministic: read explicit positions
                pos = params.pilot_positions{tx};  % np x 2 matrix [m, n]
                pm = pos(:, 1)';
                pn = pos(:, 2)';
            end

            % Mark positions as used
            for i = 1:np
                used(pm(i), pn(i)) = true;
            end

            b = b + 1;
            pb(b).dim    = 'random';
            pb(b).m      = pm;
            pb(b).n      = pn;
            pb(b).tx     = tx;
        end
    end

    % ================================================================
    %  Derived quantities
    % ================================================================
    for i = 1:length(pb)
        pb(i).n_pilots = length(pb(i).m);
    end

    n_meas_total = sum([pb.n_pilots]) * params.Nr;

    % Pilots per Tx antenna
    n_pilots_per_tx = zeros(1, params.Nt);
    for i = 1:length(pb)
        n_pilots_per_tx(pb(i).tx) = n_pilots_per_tx(pb(i).tx) + pb(i).n_pilots;
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
end
