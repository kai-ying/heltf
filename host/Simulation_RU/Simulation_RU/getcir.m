function cir2=getcir(data,framelen,rawframe)
    data=reshape(data,framelen,length(data)/framelen);
    data=fftshift(fft(data(1:end,:),[],1));
    cir2=ifft(fftshift(data.*(rawframe).'));
end
