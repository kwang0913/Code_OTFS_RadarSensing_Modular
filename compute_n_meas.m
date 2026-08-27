function n_meas = compute_n_meas(pilot_blocks, Nr)
%COMPUTE_N_MEAS Total measurement rows: n_pilots * Nr per block.
    n_meas = sum([pilot_blocks.n_pilots]) * Nr;
end
