function [c, s2] = pilot_noise_model(XTF, params, n0)
%PILOT_NOISE_MODEL Scale and noise variance of the measurement vector y = Y*conj(X)/(|X|^2+n0).
%   [c, s2] = PILOT_NOISE_MODEL(XTF, params, n0) for constant-modulus pilots
%   (|X| = params.pilot_power at every pilot/reference site):
%     y = c*Phi*alpha + w,  c = |X|^2/(|X|^2+n0),  var(w) = s2 = n0*|X|^2/(|X|^2+n0)^2.
%   LS/ridge/OMP gain fits use the model c*Phi and the ridge s2 (unit-variance gain prior).
    sites = pilot_sites(params.pilot_blocks);
    amp = abs(XTF(sub2ind(size(XTF), sites.m, sites.n, sites.tx)));
    if any(abs(amp - params.pilot_power) > 1e-9 * params.pilot_power)
        error('pilot_noise_model:notConstantModulus', ...
            'Pilots and references must have constant modulus params.pilot_power.');
    end
    p2 = params.pilot_power^2;
    c  = p2 / (p2 + n0);
    s2 = n0 * p2 / (p2 + n0)^2;
end
