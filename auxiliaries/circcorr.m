function[cc_xy] =  circcorr(x,y)
cc_xy = ifft(fft(x).*conj(fft(y))); %fftshift
%c = cconv(x,conj(fliplr(y)));
end

