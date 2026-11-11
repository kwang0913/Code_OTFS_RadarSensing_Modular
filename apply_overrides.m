function params = apply_overrides(params, overrides)
%APPLY_OVERRIDES Recursively merge override fields into a params struct.
%   Nested structs are merged field-by-field (not replaced wholesale).
%
%   Example:
%       params = apply_overrides(params, struct('ssr', struct('range', struct('a', 0.3))))
%       % Only params.ssr.range.a is changed; other ssr fields are preserved.

    fields = fieldnames(overrides);
    for i = 1:numel(fields)
        key = fields{i};
        val = overrides.(key);
        if isstruct(val) && isfield(params, key) && isstruct(params.(key))
            % Recurse into nested structs
            params.(key) = apply_overrides(params.(key), val);
        else
            params.(key) = val;
        end
    end
end
