function out = apply_rect_delay(in, delay, params, derived, adjoint, full_output)
%APPLY_RECT_DELAY Point-sampled rectangular-symbol delay and its exact adjoint.
%   Useful-time input/output has M*N rows, with one column per stream.
%   For d samples, first advance EACH source symbol by ceil(d)-d using its
%   subcarrier expansion, then delay the padded serial waveform by ceil(d).
%   This is not a whole-frame fractional FFT shift. Integer delays skip FFTs.
%   FULL_OUTPUT=true retains padding and the channel tail (forward only).
%   The fixed waveform must cover the delay; no prefix is inferred from CSI.
    if nargin < 5, adjoint = false; end
    if nargin < 6, full_output = false; end
    M = params.M; N = params.N; Q = M*N; L = derived.padLen;
    if ~isscalar(delay) || ~isfinite(delay) || ~isreal(delay) || delay < 0 || ceil(delay) > L
        error('apply_rect_delay:InvalidDelay', ...
            'Delay must be finite, nonnegative, and within the fixed padding length.');
    end
    if size(in,1) ~= Q
        error('apply_rect_delay:InvalidShape', 'Input must have M*N useful-sample rows.');
    end
    if full_output && adjoint
        error('apply_rect_delay:FullAdjoint', 'The adjoint accepts useful samples, not a padded receive frame.');
    end
    shift = ceil(delay);
    advance = shift-delay;
    streams = size(in,2);
    if adjoint
        in = integer_delay(in, shift, M, N, params.padType, true);
    end
    if adjoint, advance = -advance; end
    in = advance_rect_symbols(in,advance,M);
    if adjoint
        out = in;
    elseif full_output
        out = complex(zeros(derived.numSamps+L,streams,'like',in));
        for stream = 1:streams
            padded = apply_padding(reshape(in(:,stream),M,N),L,params.padType);
            out(shift+(1:derived.numSamps),stream) = padded;
        end
    else
        out = integer_delay(in, shift, M, N, params.padType, false);
    end
end

function out = integer_delay(in, shift, M, N, padType, adjoint)
%   Useful-window restriction of the padded linear delay. Delays fit padding,
%   so CP is circular per symbol, RCP circular per frame, and ZP has no ISI.
    if adjoint, shift = -shift; end
    switch upper(padType)
        case 'RCP'
            out = circshift(in,shift,1);
        case 'CP'
            out = reshape(circshift(reshape(in,M,N,[]),shift,1),M*N,[]);
        case 'ZP'
            out = reshape(zero_shift(reshape(in,M,N,[]),shift),M*N,[]);
        case {'RZP','NONE'}
            out = zero_shift(in,shift);
        otherwise
            error('apply_rect_delay:UnknownPadding', 'Unknown padding mode: %s.', padType);
    end
end

function out = zero_shift(in, shift)
    out = zeros(size(in),'like',in);
    count = size(in,1);
    if abs(shift) >= count, return; end
    if shift >= 0
        out(shift+1:end,:,:) = in(1:count-shift,:,:);
    else
        out(1:count+shift,:,:) = in(1-shift:end,:,:);
    end
end
