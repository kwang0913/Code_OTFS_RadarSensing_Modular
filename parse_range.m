function v = parse_range(val)
%PARSE_RANGE Convert numeric or MATLAB range string to array.
%   Examples: parse_range(5) -> 5, parse_range("2:128") -> [2 3 ... 128]
    if isnumeric(val)
        v = val;
    elseif ischar(val) || isstring(val)
        v = str2num(char(val)); %#ok<ST2NM>
    else
        error('parse_range: unsupported type');
    end
end
