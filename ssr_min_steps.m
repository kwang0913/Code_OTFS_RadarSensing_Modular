function min_step = ssr_min_steps(params)
%SSR_MIN_STEPS Smallest SSR search steps (estimator resolution) and their validation.
%   min_step = SSR_MIN_STEPS(params) returns a struct with fields a, r, v.
%   Defaults are a = params.angle_resolution, r = params.delay_resolution and
%   v = params.doppler_resolution; params.ssr.min_step.<field> overrides them
%   (empty/missing field = default).
    min_step = struct('a', params.angle_resolution, 'r', params.delay_resolution, ...
        'v', params.doppler_resolution);
    if isfield(params, 'ssr') && isfield(params.ssr, 'min_step')
        given = params.ssr.min_step;
        for name = fieldnames(min_step).'
            key = name{1};
            if isfield(given, key) && ~isempty(given.(key))
                min_step.(key) = given.(key);
            end
        end
    end
    validateattributes(min_step.a, {'numeric'}, {'real','finite','scalar','positive'}, ...
        mfilename, 'ssr.min_step.a');
    validateattributes(min_step.r, {'numeric'}, {'real','finite','scalar','positive','<=',1}, ...
        mfilename, 'ssr.min_step.r');
    validateattributes(min_step.v, {'numeric'}, {'real','finite','scalar','positive'}, ...
        mfilename, 'ssr.min_step.v');
end
