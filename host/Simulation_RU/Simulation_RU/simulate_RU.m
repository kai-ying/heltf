%%  ------ set simulation parameters
 cfg = getconfig();

%%  ------ Baseband part
bb=zeros(cfg.users,cfg.symbollen);

%get the sending he-ltf baseband for each user
for u=1:cfg.users
    [bb(u,:),cfg.userltf(u,:)]=getheltf(cfg,u);
%     bb(u,1:10) = 0;
%     bb(u,end-10:end) = 0;
%     plot(real(bb(u,:)));
%     hold on
end
% 全部子载波
transmit_bb = bb(1,:)./ max(imag(bb(1,:)));
figure()
plot(real(transmit_bb))
hold on
plot(imag(transmit_bb))
hold on
% 分布式子载波
% transmit_bb = sum(bb(:,:));

figure
cfr_bb=fftshift(fft(transmit_bb(cfg.gilen+1:end)));
plot(abs(cfr_bb))
hold on


imag_bb = imag(transmit_bb);
real_bb = real(transmit_bb);
writematrix(real_bb, 'real_bb.txt')
writematrix(imag_bb, 'imag_bb.txt')


%%  ------ Simple baseband channel model, will consider real up/down conversion later

%channel parameters

doppler=zeros(cfg.users,1);
doppler(1)=3000; % doppler
doppler(2)=-2000;
doppler(3)=1000;
doppler(4)=2000;

delay=zeros(cfg.users,3); %simple delay profile, each user 3 taps
delay(1,:)=[130 200 330]; %delay in ns 
delay(2,:)=[250 350 400];
delay(3,:)=[200 250 300];
delay(4,:)=[250 290 350];

gain=zeros(cfg.users,3); %simple attenuation profile
gain(1,:)=[ 1*exp(-1j*pi/4) 0.3*exp(-1j*5*pi/4) 0.4*exp(-1j*pi/2)];
gain(2,:)=[ 3*exp(-1j*pi/3) 0.7*exp(-1j*13*pi/4) 0.3*exp(-1j*pi/4)];
gain(3,:)=[ 1*exp(-1j*2*pi/3) 0.3*exp(-1j*3*pi/4) 0.4*exp(-1j*7*pi/4)];
gain(4,:)=[ 0.8*exp(-1j*pi) 2*exp(-1j*5*pi/3) 0.5*exp(-1j*4*pi/3)];

recv_bb=zeros(1,cfg.symbollen); %received baseband

for u=1:cfg.users
    recv_bb=recv_bb+sim_channel(cfg,bb(u,:),doppler(u),delay(u,:),gain(u,:));
%     plot(real(recv_bb));
%     hold on
end

%%  --decode part
figure

plot(abs(cfr_all))
legend('receive')
% hold on
% plot(real(cfr_all))
% hold on 
% plot(imag(cfr_all))
% hold on

cir=zeros(cfg.users,cfg.fftlen);
figure
% get CIR for each user
for u=1:cfg.users
    cir(u,:)=ifft( fftshift(cfr_all.*conj(cfg.userltf(u,:)))); %convolution
    plot(abs(cir(u,:)));
    hold on
%     plot(abs(fft(cir(2,:))))
%     hold on
end