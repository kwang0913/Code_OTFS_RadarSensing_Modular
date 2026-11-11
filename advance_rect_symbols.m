function out = advance_rect_symbols(in, advance, M)
%ADVANCE_RECT_SYMBOLS Evaluate each rectangular source symbol at q+advance.
%   IN is M*N-by-stream; frequency bins are m=0:M-1, not centered bins.
%   Negative advance is the Euclidean adjoint. Zero bypasses all transforms.
    if advance == 0
        out = in;
        return;
    end
    phase = exp(1j*2*pi*(0:M-1).'*advance/M);
    out = reshape(ifft(fft(reshape(in,M,[],size(in,2)),[],1).*phase,[],1),size(in));
end
