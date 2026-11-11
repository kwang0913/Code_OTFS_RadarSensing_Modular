function [out, raw] = apply_time_channel(in, chan, params, derived, adjoint)
%APPLY_TIME_CHANNEL Matrix-free useful-time MIMO channel or Euclidean adjoint.
%   Forward: M*N-by-Nt -> M*N-by-Nr. Adjoint reverses these shapes.
%   Optional forward second output retains padding and the channel tail,
%   (derived.numSamps+derived.padLen)-by-Nr, for physical AWGN simulation.
%   Delays remain fractional. Padding and Doppler timing are fixed waveform
%   metadata, shared by truth and estimates, never recomputed from chan.
%   Doppler multiplies AFTER delay, at the original receive sample index.
%   Gains: chan.pathGains is one gain per path and Tx steering uses chan.AoDs.
%   With chan.gain_per_tx = true, chan.pathGains is paths-by-Nt (one gain per
%   path and Tx antenna) and replaces the Tx steering; AoDs is then unused.
    if nargin < 5, adjoint = false; end
    M = params.M; N = params.N; Q = M*N; L = derived.padLen;
    [Meff, numSamps, offset] = compute_padding_params(M,N,L,params.padType);
    if ~isscalar(L) || ~isfinite(L) || L < 0 || L ~= round(L) || ...
            derived.numSamps ~= numSamps || derived.offset ~= offset || ...
            derived.Meff ~= Meff || ...
            (strcmpi(params.padType,'CP') && L > M) || ...
            (strcmpi(params.padType,'RCP') && L > Q)
        error('apply_time_channel:InvalidWaveform', ...
            'Waveform metadata must describe the fixed transmitted padding and useful window.');
    end
    expected_antennas = params.Nt;
    if adjoint, expected_antennas = params.Nr; end
    if ~isequal(size(in),[Q,expected_antennas])
        error('apply_time_channel:InvalidShape', 'Input has incorrect useful-sample or antenna dimensions.');
    end
    full_output = nargout > 1;
    if full_output && adjoint
        error('apply_time_channel:FullAdjoint', 'The adjoint returns only useful transmit samples.');
    end
    per_tx_gains = isfield(chan,'gain_per_tx') && chan.gain_per_tx;
    if per_tx_gains
        if ~isequal(size(chan.pathGains),[numel(chan.pathDelays),params.Nt])
            error('apply_time_channel:InvalidShape', ...
                'Per-Tx gains must be paths-by-Nt.');
        end
        at = chan.pathGains.';
    else
        at = exp(-1j*2*pi*(0:params.Nt-1).'*params.ddlt*sind(chan.AoDs(:).'));
    end
    ar = exp(1j*2*pi*(0:params.Nr-1).'*params.ddlr*sind(chan.AoAs(:).'));
    time = reshape((0:M-1).'+(0:N-1)*Meff+offset,Q,1);
    if full_output, time = (0:numSamps+L-1).'; end
    denominator = derived.fsamp*derived.Tsamp;
    if adjoint
        signals = in*conj(ar);
    else
        signals = in*at;
    end
    if full_output
        padded_signals = complex(zeros(numSamps+L,numel(chan.pathDelays),'like',in));
    end
    for j = 1:numel(chan.pathDelays)
        phase = exp(1j*2*pi*chan.pathDopplers(j)*time/denominator);
        if ~per_tx_gains, phase = chan.pathGains(j)*phase; end
        if adjoint
            signals(:,j) = apply_rect_delay(conj(phase).*signals(:,j), ...
                chan.pathDelays(j),params,derived,true);
        elseif full_output
            padded_signals(:,j) = phase.*apply_rect_delay(signals(:,j), ...
                chan.pathDelays(j),params,derived,false,true);
        else
            signals(:,j) = phase.*apply_rect_delay(signals(:,j), ...
                chan.pathDelays(j),params,derived,false);
        end
    end
    if adjoint
        out = signals*at';
    elseif full_output
        raw = padded_signals*ar.';
        out = complex(zeros(Q,params.Nr,'like',in));
        for r = 1:params.Nr
            useful = strip_padding(raw(1:numSamps,r),M,N,L,offset,params.padType);
            out(:,r) = useful(:);
        end
    else
        out = signals*ar.';
    end
end
