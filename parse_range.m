function arr = parse_range(val)
%PARSE_RANGE Convert a numeric value or MATLAB-style range string to an array.
%
%   arr = parse_range(5)         % -> 5
%   arr = parse_range([1 3 5])   % -> [1 3 5]
%   arr = parse_range('2:128')   % -> 2:128
%   arr = parse_range('1:2:11')  % -> [1 3 5 7 9 11]

    if isnumeric(val)
        arr = val(:).';
        return
    end

    try
        arr = eval(char(val));
    catch ex
        error('parse_range:invalid', 'Cannot parse range string "%s": %s', val, ex.message);
    end
end
