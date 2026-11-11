function derivedParams = compute_derived_params(params)
%COMPUTE_DERIVED_PARAMS Sampling, padding and resolution quantities from params alone.
%   params.range_factor: 2 monostatic (round trip), 1 bistatic; scales range and velocity resolution.
    derivedParams.fsamp = params.M * params.df;
    derivedParams.padLen = params.padLen;
    [derivedParams.Meff, derivedParams.numSamps, derivedParams.offset] = ...
        compute_padding_params(params.M, params.N, params.padLen, params.padType);
    derivedParams.dt = derivedParams.Meff / derivedParams.fsamp;
    derivedParams.Tsamp = params.N * derivedParams.dt;
    derivedParams.T_s = 1 / derivedParams.fsamp;
    derivedParams.lambda = params.c / params.fc;
    derivedParams.range_res = params.c / (params.range_factor * params.M * params.df);
    derivedParams.velocity_res = derivedParams.lambda / (params.range_factor * derivedParams.Tsamp);
end
