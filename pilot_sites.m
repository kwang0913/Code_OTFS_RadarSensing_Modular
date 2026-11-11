function sites = pilot_sites(pilot_blocks)
%PILOT_SITES Unique pilot/reference TF sites of a pilot layout, in measurement order.
%   sites = PILOT_SITES(pilot_blocks) lists every (tx, m, n) TF site that
%   carries a reference or a pilot ONCE, even when several blocks, or a
%   block's own reference and first pilot, share it. Single source of truth
%   for measurement rows, dictionary rows and guard counts.
%
%   Order: block by block, each block's reference first, then its pilots in
%   pilot order; sites already seen are skipped.
%   Fields (all indices 1-based):
%     tx, m, n      - S-by-1 site coordinates
%     n_sites       - S, number of unique (tx,m,n) sites
%     n_bins        - number of unique TF bins (m,n), independent of tx; this
%                     is the number of data symbols guarded per Tx antenna
%     ref_site(b)   - site index of block b's reference
%     pilot_site{b} - n_pilots(b)-by-1 site indices of block b's pilots
%
%   Pilot p of a block: 'random' -> (m(p),n(p)); 'delay' -> column-major
%   over (m,n), m fastest; otherwise ('doppler') -> (m, n(p)).
    num_blocks = numel(pilot_blocks);
    raw = cell(num_blocks, 1);
    first_pilot = zeros(num_blocks, 1);
    count = 0;
    for b = 1:num_blocks
        blk = pilot_blocks(b);
        p = (1:blk.n_pilots).';
        switch blk.dim
            case 'random'
                pm = blk.m(p); pn = blk.n(p);
            case 'delay'
                num_m = numel(blk.m);
                pm = blk.m(mod(p-1, num_m) + 1);
                pn = blk.n(floor((p-1) / num_m) + 1);
            otherwise
                pm = repmat(blk.m(1), numel(p), 1);
                pn = blk.n(p);
        end
        raw{b} = [blk.ref_tx, blk.ref_m, blk.ref_n; ...
                  repmat(blk.tx, numel(p), 1), pm(:), pn(:)];
        first_pilot(b) = count + 2;
        count = count + 1 + numel(p);
    end
    all_rows = vertcat(raw{:});
    [~, first_idx, site_of_row] = unique(all_rows, 'rows', 'stable');

    sites.tx = all_rows(first_idx, 1);
    sites.m = all_rows(first_idx, 2);
    sites.n = all_rows(first_idx, 3);
    sites.n_sites = numel(first_idx);
    sites.n_bins = size(unique(all_rows(:, 2:3), 'rows'), 1);
    sites.ref_site = zeros(num_blocks, 1);
    sites.pilot_site = cell(num_blocks, 1);
    for b = 1:num_blocks
        n_rows = size(raw{b}, 1);
        rows = first_pilot(b) - 1 + (0:n_rows-1);
        sites.ref_site(b) = site_of_row(rows(1));
        sites.pilot_site{b} = site_of_row(rows(2:end));
    end
end
