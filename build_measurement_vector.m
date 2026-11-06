function y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr, n0)
%BUILD_MEASUREMENT_VECTOR Construct regularized Y/X ratio vector from pilot blocks.
%   y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr, n0)
%
%   Layout: [pilot_1(Nr) | pilot_2(Nr) | ...] per block, concatenated.
%   Division is Wiener-regularized, Y*conj(X)/(|X|^2 + n0), so weak pilots
%   do not amplify data interference and noise.

    n_meas = compute_n_meas(pilot_blocks, Nr);
    y = zeros(n_meas, 1);
    idx = 0;
    for b = 1:length(pilot_blocks)
        blk = pilot_blocks(b);
        for p = 1:blk.n_pilots
            X = XTF(blk.m(p), blk.n(p), blk.tx);
            Y = squeeze(YTF(blk.m(p), blk.n(p), 1:Nr));
            y(idx + (1:Nr)) = Y * conj(X) / (abs(X)^2 + n0);
            idx = idx + Nr;
        end
    end
end
