function y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr)
%BUILD_MEASUREMENT_VECTOR Construct Y/X ratio vector from pilot blocks.
%   y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr)
%
%   Layout: [pilot_1(Nr) | pilot_2(Nr) | ...] per block, concatenated.

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = zeros(n_meas, 1);
    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        for p = 1:blk.n_pilots
            y(idx + (1:Nr)) = squeeze(YTF(blk.m(p), blk.n(p), 1:Nr)) ./ XTF(blk.m(p), blk.n(p), blk.tx);
            idx = idx + Nr;
        end
    end
end
