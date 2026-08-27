function base = apply_overrides(base, overrides)
%APPLY_OVERRIDES Recursively merge override fields into base struct.
    if isempty(overrides), return; end
    fields = fieldnames(overrides);
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(overrides.(f)) && isfield(base, f) && isstruct(base.(f))
            base.(f) = apply_overrides(base.(f), overrides.(f));
        else
            base.(f) = overrides.(f);
        end
    end
end
