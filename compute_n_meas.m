function n_meas = compute_n_meas(pilot_blocks, Nr)
%COMPUTE_N_MEAS Total measurement rows: Nr per UNIQUE pilot/reference TF site.
%   A reference that coincides with a pilot (or with another block's site)
%   is one observation; see pilot_sites.
    n_meas = pilot_sites(pilot_blocks).n_sites * Nr;
end
