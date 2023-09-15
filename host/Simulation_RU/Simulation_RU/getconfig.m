function cfg=getconfig()

cfg.fc=5.825e9; %carrier frequency
cfg.users=2;    %number of users
cfg.ltfmode=2;  %HE-LTF mode : 1- x1; 2- x2; 3- x4;
cfg.user_offset=[0 -1 -2 1]; %user subcarrier offset, for 20M x4 should be -2,-1,0,1 other x4 should be 0,1,2,3
cfg.user_phase=[0 pi/2 pi -pi/2];
%cfg.user_phase=[0 0 0 0];
subinterval=[4 2 1]; %interval between subcarriers
cfg.user_phase=cfg.user_phase/subinterval(cfg.ltfmode);
fftlength=[256 256 256];   %FFT length for each mode in 20MHz band, always use x4 for ax
cfg.bandwidth=4; %bandwidth: 1- 20M; 2- 40M; 3: 80M; 4: 160M;
realbw=[20e6 40e6 80e6 160e6]; %real bandwidth
cfg.bw=realbw(cfg.bandwidth);
fftlenmul=[1 2 4 8];
cfg.fftlen=fftlength(cfg.ltfmode)*fftlenmul(cfg.bandwidth);
cfg.gitype=2;    %GI length: 1-0.8us; 2 1.6us; 3-3.2us
gilength=[16 32 64];
cfg.gilen=gilength(cfg.gitype)*fftlenmul(cfg.bandwidth);
cfg.symbollen=cfg.gilen+cfg.fftlen;
cfg.userltf=zeros(cfg.users,cfg.fftlen);
end