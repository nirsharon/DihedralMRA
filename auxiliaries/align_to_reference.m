function x_aligned = align_to_reference(x, xref)

assert(all(size(x) == size(xref)), 'x and xref must have identical size');
x = x(:);
L = length(x);
xref = xref(:);

x_fft = fft(x);
xref_fft = fft(xref);

correlation_x_xref = real(ifft(conj(x_fft) .* xref_fft));
correlation_rev_x_xref = real(ifft(conj(x_fft) .* conj(xref_fft)));
[max1, ind1] = max(correlation_x_xref);
[max2, ind2] = max(correlation_rev_x_xref);
if max1>max2
    x_aligned = circshift(x, ind1-1);
else
    x_aligned = circshift(reverse(x), -(ind2-1));
end
end