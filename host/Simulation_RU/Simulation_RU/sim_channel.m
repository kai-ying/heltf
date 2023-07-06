function sig=sim_channel(cfg,bb,doppler,delay,gain)
bblen = length(bb);
taps=length(delay);
temp=bb.*exp(-1j * (0 : (bblen-1)).*2*pi.*doppler./cfg.bw);   %simulate the doppler shift, per sample point phase change
sig=zeros(size(bb));
%plot(abs(temp))
%hold on
for i=1:taps
    d=zeros(size(bb));
    dd=round(delay(i)*1e-9*cfg.bw); %here only use the simple integer delay, will consider real non-integer delay later
    d(1+dd:end)=temp(1:end-dd);
    sig=sig+d*gain(i)*exp(-1j*2*pi*cfg.fc*delay(i)*1e-9);
end
end

