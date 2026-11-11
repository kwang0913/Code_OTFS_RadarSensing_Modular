function y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr, n0)
%BUILD_MEASUREMENT_VECTOR Construct regularized Y/X ratio vector from pilot blocks.
%   y = build_measurement_vector(YTF, XTF, pilot_blocks, Nr, n0)
%
%   Layout: Nr rows per unique pilot/reference TF site, in pilot_sites order
%   (block by block: reference, then new pilot sites). Row (s-1)*Nr+r is
%   receive antenna r at site s divided by the symbol transmitted there.
%   Division is Wiener-regularized, Y*conj(X)/(|X|^2 + n0), so weak pilots
%   do not amplify data interference and noise.

    sites = pilot_sites(pilot_blocks);
    M = size(YTF, 1);
    received = reshape(YTF, M*size(YTF, 2), []);
    received = received(sites.m + (sites.n-1)*M, 1:Nr);                       % S-by-Nr
    transmitted = XTF(sub2ind(size(XTF), sites.m, sites.n, sites.tx));       % S-by-1
    y = reshape((received .* conj(transmitted) ./ (abs(transmitted).^2 + n0)).', [], 1);
end
