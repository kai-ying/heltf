%%//2022.5.26更新后代码
%%//2022.6.12更新后代码
clear all;
close all;
Td=0.0001; %可容忍时延
F=1e6; %需要实现的速率
bit_per_zidai=4;%每个子带一个符号承载的bit数目

T_cp=4*Td;%循环前缀的时延
T_zidai=6*T_cp;%子带的周期
F_zidai=1/(5*T_cp);%子带带宽间隔
bit_per_fuhao=F/(1/T_zidai);%划分子带前一个符号要传输的bit
carrier_count = round(bit_per_fuhao/bit_per_zidai); % 子载波数
symbol_count = round(1/T_zidai); %总符号数
ifft_length = 2^nextpow2(carrier_count); % IFFT长度
CP_length = carrier_count/5; % 循环前缀

fc=7e6; %载波频率
Fs=6e6;%采样频率
rate = [];%误码率统计
SNR=30;%典型信噪比
SNR_range =10:2:20;%信噪比范围
bit_per_symbol = 4;%每个符号的bit数目
alpha = 0.5; % 根升余弦窗系数
h_mut=F_zidai*(carrier_count+1)*Td;%多径数目
Am = [0.1,0.2];%多径幅度
Td_bit_max=round(h_mut);%最大bit延迟

% ================1.产生随机序列=======================

bit_length = carrier_count*symbol_count*bit_per_symbol;
bit_sequence = round(rand(1,bit_length))'; % 列向量

% =================2.串并转换==========================
% ==================16QAM调制=========================

% carrier_position = 1:200;%把ifft后直流分量的地方空出来
bit_moded = qammod(bit_sequence,16,'InputType','bit');

figure(1);
scatter(real(bit_moded),imag(bit_moded));
title('调制后的星座图');
grid on;

bit_moded = reshape(bit_moded,carrier_count,symbol_count);

%===================加入导频==========================
Nps=16;             %导频间隔
Np=ceil(ifft_length/Nps)  %导频个数
Xp = round(rand(1,symbol_count*Np*bit_per_symbol))';  %随机产生导频
Xp_moded = qammod(Xp,16,'InputType','bit');
Xp_moded = reshape(Xp_moded,Np,symbol_count);
bit_add_pilot=zeros(ifft_length,symbol_count);
bit_prep=zeros(ifft_length,symbol_count);
bit_prep(1:carrier_count,:)=bit_moded(1:carrier_count,:);
ip=0;
pilot_loc=[];   %导频位置
for k=1:ifft_length
    if mod(k,Nps)==1
        bit_add_pilot(k,:) = Xp_moded(floor(k/Nps)+1,:); 
        pilot_loc = [pilot_loc k]; 
        ip = ip+1;
    else        
        bit_add_pilot(k,:) = bit_prep(k-ip,:);
    end
end

bit_moded = bit_add_pilot;

% ===================3.IFFT===========================

signal_time = ifft(bit_moded,ifft_length);

% ==================4.加循环前缀==================

signal_time_C = [signal_time(end-CP_length+1:end,:);signal_time]; % 单个完整符号为1024+128=1152
signal_window=signal_time_C;

% ================5.根升成型滤波器设计========== 

delay=10;                  %当前码元影响周围码元的个数
Nsam=Fs/F;                %每个符号抽样的数
hsqrt = rcosdesign(alpha,delay,Fs/F,'sqrt');  % 产生平方根升余弦滤波器
Nh=0.5*(length(hsqrt)-1);                   %滤波器长度的一半 

% ================6.信号成型滤波=================

for i=1:size(signal_window,1)
signal_time_C_up(i,:)= upsample(signal_window(i,:),Fs/F);      %内插Nsam-1个0,一码元多采样，产生基带脉冲信号
signal_time_C_f(i,:)=conv(signal_time_C_up(i,:),hsqrt);        %滤波并截断
signal_rcos(i,:)=signal_time_C_f(i,1+Nh:end-Nh);
end

figure;
ly = length(signal_rcos(1,:));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];
Syy = 10*log10(fftshift(abs(fft(signal_rcos(1,:))/Fs)));
plot(freq,Syy);
title('成形后基带子带信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

% ===================7.发送信号，多径信道====================

signal_Tx = reshape(signal_rcos,1,[]); % 并串转换，变成时域一个完整信号，待传输

figure;
ly = length(signal_Tx(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(abs(fftshift(fft(signal_Tx(1:carrier_count))/Fs))); 
plot(freq,Syy);
title('成形后基带信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

% 上变频
loc_carr = exp(1i*2*pi*(0:length(signal_Tx)-1)*fc/Fs);
signal_Tx = real(signal_Tx.*loc_carr).*2;

figure;
ly = length(signal_Tx(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(signal_Tx(1:carrier_count))/Fs))); 
plot(freq,Syy);
title('成形后频带单径信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

h=[(randn+j*randn),(randn+j*randn)/2];
% path2 = Am(1)*[zeros(1,round(Td_bit_max/2)) signal_Tx(1:end-round(Td_bit_max/2)) ];
% path3 = Am(2)*[zeros(1,Td_bit_max) signal_Tx(1:end-Td_bit_max) ];
% signal_Tx_mult =signal_Tx + path2 + path3 ; % 多径信号
signal_Tx_mult =conv(signal_Tx,h)  ; % 多径信号（另一种方式）
signal_Tx_mult = signal_Tx_mult(1:end-1);

figure;
ly = length(signal_Tx_mult(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(signal_Tx_mult(1:carrier_count))/Fs))); 
plot(freq,Syy);
title('成形后频带多径信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

% ========================8.加AWGN==========================

signal_power_sig = var(signal_Tx); % 单径发送信号功率
signal_power_mut = var(signal_Tx_mult); % 多径发送信号功率
SNR_linear = 10^(SNR/10);
noise_power_mut = signal_power_mut/SNR_linear;
noise_power_sig = signal_power_sig/SNR_linear;
noise_sig = randn(size(signal_Tx))*sqrt(noise_power_sig);
noise_mut = randn(size(signal_Tx_mult))*sqrt(noise_power_mut);
Rx_data_sig = signal_Tx+noise_sig;
Rx_data_mut = signal_Tx_mult+noise_mut;


figure;
ly = length(Rx_data_sig(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(Rx_data_sig(1:carrier_count))/Fs)));   
plot(freq,Syy);
title('成形后频带单径加噪信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

figure;
ly = length(Rx_data_mut(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(Rx_data_mut(1:carrier_count))/Fs)));   
plot(freq,Syy);
title('成形后频带多径加噪信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

% 下变频
loc_carr = exp(-1i*2*pi*(0:length(Rx_data_sig)-1)*fc/Fs);
Rx_data_sig =Rx_data_sig.*loc_carr;
Rx_data_mut =Rx_data_mut.*loc_carr;

figure;
ly = length(Rx_data_sig(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(Rx_data_sig(1:carrier_count))/Fs)));   
plot(freq,Syy);
title('成形后基带单径加噪信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

figure;
ly = length(Rx_data_mut(1:carrier_count));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(Rx_data_mut(1:carrier_count))/Fs)));   
plot(freq,Syy);
title('成形后基带多径加噪信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

% =======================接收=============================
% =======================1.串并转换==========================

Rx_data_mut = reshape(Rx_data_mut,ifft_length+CP_length,[]);
Rx_data_sig = reshape(Rx_data_sig,ifft_length+CP_length,[]);
Rx_data_sig_final=Rx_data_sig;
Rx_data_mut_final=Rx_data_mut;

% =========================2.恢复采样===========================

%接收滤波
for i=1:size(Rx_data_sig_final,1)     
Rx_data_sig_f(i,:)=conv(Rx_data_sig_final(i,:),hsqrt);  
Rx_data_mut_f(i,:)=conv(Rx_data_mut_final(i,:),hsqrt);  %滤波并截断
Rx_data_sig_cut(i,:)=Rx_data_sig_f(i,1+Nh:end-Nh);
Rx_data_mut_cut(i,:)=Rx_data_mut_f(i,1+Nh:end-Nh);
end

figure;
ly = length(Rx_data_sig_cut(1,:));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(Rx_data_sig_cut(1,:))/Fs)));   
plot(freq,Syy);
title('匹配滤波后基带单径加噪信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

figure;
ly = length(Rx_data_mut_cut(1,:));                            
freq = [-Fs/2 : Fs/ly : Fs/2 - Fs/ly];  
Syy = 10*log10(fftshift(abs(fft(Rx_data_mut_cut(1,:))/Fs)));   
plot(freq,Syy);
title('匹配滤波后基带多径加噪信号的功率谱密度');xlabel('频率'),ylabel('功率谱(dB)');

for i=1:size(Rx_data_mut_cut,1)
Rx_data_sig_last(i,:)=downsample(Rx_data_sig_cut(i,:),Nsam);
Rx_data_mut_last(i,:)=downsample(Rx_data_mut_cut(i,:),Nsam);   %最佳采样点采样
end

%========================3.加CFO==============================

CFO=0.15;
nn_mut_CFO=0:length(Rx_data_mut_last)-1; 
nn_mut_CFO=repmat(nn_mut_CFO.',1,symbol_count);
Rx_data_mut_CFO = Rx_data_mut_last.*exp(1i*2*pi*CFO*nn_mut_CFO/ifft_length);


nn_sig_CFO=0:length(Rx_data_sig_last)-1; 
nn_sig_CFO=repmat(nn_sig_CFO.',1,symbol_count);
Rx_data_sig_CFO = Rx_data_sig_last.*exp(1i*2*pi*CFO*nn_sig_CFO/ifft_length);


%=========================4.估计CFO===========================
%=========================使用cp估计========================
MSE_CFO_CP_mut=0;
MSE_CFO_CP_sig=0;

for iter=1:100
    nn_est_CFO_sig=1:CP_length; 
    CFO_est_sig = angle(Rx_data_sig_CFO(nn_est_CFO_sig,1)'*Rx_data_sig_CFO(nn_est_CFO_sig+ifft_length,1))/(2*pi);
    MSE_CFO_CP_sig=MSE_CFO_CP_sig+(CFO_est_sig-CFO)^2;
    
    nn_est_CFO_mut=1:CP_length; 
    CFO_est_mut = angle(Rx_data_mut_CFO(nn_est_CFO_mut,1)'*Rx_data_mut_CFO(nn_est_CFO_mut+ifft_length,1))/(2*pi);
    MSE_CFO_CP_mut=MSE_CFO_CP_mut+(CFO_est_mut-CFO)^2;
end
MSE_CFO_CP_mut=MSE_CFO_CP_mut/100
MSE_CFO_CP_sig=MSE_CFO_CP_sig/100


%=========================5.去除CFO===========================

Rx_data_mut_deCFO = Rx_data_mut_CFO.*exp(-1i*2*pi*CFO_est_mut*nn_mut_CFO/ifft_length);

Rx_data_sig_deCFO = Rx_data_sig_CFO.*exp(-1i*2*pi*CFO_est_sig*nn_sig_CFO/ifft_length);

% =======================6.串并转换==========================

Rx_data_sig_final=Rx_data_sig_deCFO;
Rx_data_mut_final=Rx_data_mut_deCFO;

% ====================7.去循环前缀和后缀======================

Rx_data_sig_final(1:CP_length,:) = [];
Rx_data_mut_final(1:CP_length,:) = [];

% =========================8.FFT=============================

fft_sig = fft(Rx_data_sig_final);
fft_mut = fft(Rx_data_mut_final);

%========================信道估计============================

for i=1:symbol_count
    H_mut(:,i)=LS_CE(fft_mut(:,i).',Xp_moded(:,i).',pilot_loc,ifft_length,Nps,'spline');
end
H = fft(h,ifft_length);
H_realtime = fft_mut./bit_moded;
fft_mut_ref = fft_mut;
fft_mut = fft_mut./H_mut;


MSE_realtime=[];
MSEs_realtime=[];
MSE_est=[];
MSEs_est=[];
for i=1:symbol_count
    MSE_est(i)=(H-H_mut(:,i).')*(H-H_mut(:,i).')';
    MSEs_est(i)=MSE_est(i)/symbol_count;
    MSE_realtime(i)=(H(1:(carrier_count+floor(carrier_count/Nps)))-H_realtime(1:(carrier_count+floor(carrier_count/Nps)),i).')*(H(1:(carrier_count+floor(carrier_count/Nps)))-H_realtime(1:(carrier_count+floor(carrier_count/Nps)),i).')';
    MSEs_realtime(i)=MSE_realtime(i)/symbol_count;
end
MSE_H_realtime=sum(MSE_realtime)/(carrier_count+floor(carrier_count/Nps));
MSE_H_realtime_point=sum(MSEs_realtime)/(carrier_count+floor(carrier_count/Nps));
MSE_H_est=sum(MSE_est)/(carrier_count+floor(carrier_count/Nps));
MSE_H_est_point=sum(MSEs_est)/(carrier_count+floor(carrier_count/Nps));

%=========================去除导频===========================

fft_sig(pilot_loc,:) = [];
fft_mut(pilot_loc,:) = [];
fft_mut_ref(pilot_loc,:) = [];

% =========================9.采样（去除ifft中补的位数）===========

data_sig = fft_sig(1:carrier_count,:);
data_mut = fft_mut(1:carrier_count,:);
data_mut_ref = fft_mut_ref(1:carrier_count,:);

figure(6)
scatter(real(reshape(data_sig,1,[])),imag(reshape(data_sig,1,[])),'.')
grid on;
title('单径接收信号星座图')

figure(7)
scatter(real(reshape(data_mut_ref,1,[])),imag(reshape(data_mut_ref,1,[])),'.')
grid on;
title('多径未矫正信道接收信号星座图')

figure(8)
scatter(real(reshape(data_mut,1,[])),imag(reshape(data_mut,1,[])),'.')
grid on;
title('多径信道矫正后接收信号星座图')

% =========================10.16QAM逆映射===========================

bit_demod_sig = reshape(qamdemod(data_sig,16,'OutputType','bit'),[],1);
bit_demod_mut = reshape(qamdemod(data_mut,16,'OutputType','bit'),[],1);
bit_demod_mut_ref = reshape(qamdemod(data_mut_ref,16,'OutputType','bit'),[],1);
axis([-5 5 -5 5])

% =========================11.误码率===========================
error_bit_sig = sum(bit_demod_sig~=bit_sequence);
error_bit_mut = sum(bit_demod_mut~=bit_sequence);
error_bit_mut_ref = sum(bit_demod_mut_ref~=bit_sequence);
error_rate_sig = error_bit_sig/bit_length;
error_rate_mut = error_bit_mut/bit_length;
error_rate_mut_ref = error_bit_mut_ref/bit_length;
rate = [error_rate_sig error_rate_mut error_rate_mut_ref];
disp(['单径误码率为：',num2str(rate(1))]);
disp(['多径信道未校正时误码率为：',num2str(rate(3))])
disp(['多径信道LS(spline插值)矫正后的误码率为：',num2str(rate(2))])
rate_comp=[];
for SNR=SNR_range %不同信噪比下的误码率对比

% ========================8.加AWGN==========================

signal_power_sig = var(signal_Tx); % 单径发送信号功率
signal_power_mut = var(signal_Tx_mult); % 多径发送信号功率
SNR_linear = 10^(SNR/10);
noise_power_mut = signal_power_mut/SNR_linear;
noise_power_sig = signal_power_sig/SNR_linear;
noise_sig = randn(size(signal_Tx))*sqrt(noise_power_sig);
noise_mut = randn(size(signal_Tx_mult))*sqrt(noise_power_mut);
Rx_data_sig = signal_Tx+noise_sig;
Rx_data_mut = signal_Tx_mult+noise_mut;

% 下变频
loc_carr = exp(-1i*2*pi*(0:length(Rx_data_sig)-1)*fc/Fs);
Rx_data_sig =Rx_data_sig.*loc_carr;
Rx_data_mut =Rx_data_mut.*loc_carr;

% =======================接收=============================
% =======================1.串并转换==========================

Rx_data_mut = reshape(Rx_data_mut,ifft_length+CP_length,[]);
Rx_data_sig = reshape(Rx_data_sig,ifft_length+CP_length,[]);
Rx_data_sig_final=Rx_data_sig;
Rx_data_mut_final=Rx_data_mut;

% =========================2.恢复采样===========================

%接收滤波
for i=1:size(Rx_data_sig_final,1)     
Rx_data_sig_f(i,:)=conv(Rx_data_sig_final(i,:),hsqrt);  
Rx_data_mut_f(i,:)=conv(Rx_data_mut_final(i,:),hsqrt);  %滤波并截断
Rx_data_sig_cut(i,:)=Rx_data_sig_f(i,1+Nh:end-Nh);
Rx_data_mut_cut(i,:)=Rx_data_mut_f(i,1+Nh:end-Nh);
end

for i=1:size(Rx_data_mut_cut,1)
Rx_data_sig_last(i,:)=downsample(Rx_data_sig_cut(i,:),Nsam);
Rx_data_mut_last(i,:)=downsample(Rx_data_mut_cut(i,:),Nsam);   %最佳采样点采样
end

%========================3.加CFO==============================

CFO=0.15;
nn_mut_CFO=0:length(Rx_data_mut_last)-1; 
nn_mut_CFO=repmat(nn_mut_CFO.',1,symbol_count);
Rx_data_mut_CFO = Rx_data_mut_last.*exp(1i*2*pi*CFO*nn_mut_CFO/ifft_length);

nn_sig_CFO=0:length(Rx_data_sig_last)-1; 
nn_sig_CFO=repmat(nn_sig_CFO.',1,symbol_count);
Rx_data_sig_CFO = Rx_data_sig_last.*exp(1i*2*pi*CFO*nn_sig_CFO/ifft_length);

%=========================4.估计CFO===========================
%=========================使用cp估计========================

MSE_CFO_CP_mut=0;
MSE_CFO_CP_sig=0;

for iter=1:100
    nn_est_CFO_sig=1:CP_length; 
    CFO_est_sig = angle(Rx_data_sig_CFO(nn_est_CFO_sig,1)'*Rx_data_sig_CFO(nn_est_CFO_sig+ifft_length,1))/(2*pi);
    MSE_CFO_CP_sig=MSE_CFO_CP_sig+(CFO_est_sig-CFO)^2;
    
    nn_est_CFO_mut=1:CP_length; 
    CFO_est_mut = angle(Rx_data_mut_CFO(nn_est_CFO_mut,1)'*Rx_data_mut_CFO(nn_est_CFO_mut+ifft_length,1))/(2*pi);
    MSE_CFO_CP_mut=MSE_CFO_CP_mut+(CFO_est_mut-CFO)^2;
end
MSE_CFO_CP_mut=MSE_CFO_CP_mut/100;
MSE_CFO_CP_sig=MSE_CFO_CP_sig/100;


%=========================5.去除CFO===========================

Rx_data_mut_deCFO = Rx_data_mut_CFO.*exp(-1i*2*pi*CFO_est_mut*nn_mut_CFO/ifft_length);
Rx_data_sig_deCFO = Rx_data_sig_CFO.*exp(-1i*2*pi*CFO_est_sig*nn_sig_CFO/ifft_length);

% =======================6.串并转换==========================

Rx_data_sig_final=Rx_data_sig_deCFO;
Rx_data_mut_final=Rx_data_mut_deCFO;

% ====================7.去循环前缀和后缀======================

Rx_data_sig_final(1:CP_length,:) = [];
Rx_data_mut_final(1:CP_length,:) = [];

% =========================8.FFT=============================

fft_sig = fft(Rx_data_sig_final);
fft_mut = fft(Rx_data_mut_final);

%========================信道估计============================

for i=1:symbol_count
    H_mut_s(:,i)=LS_CE(fft_mut(:,i).',Xp_moded(:,i).',pilot_loc,ifft_length,Nps,'spline');
end
fft_mut_s = fft_mut./H_mut_s;
for i=1:symbol_count
    H_mut_l(:,i)=LS_CE(fft_mut(:,i).',Xp_moded(:,i).',pilot_loc,ifft_length,Nps,'linear');
end
fft_mut_l = fft_mut./H_mut_l;

%=========================去除导频===========================

fft_sig(pilot_loc,:) = [];
fft_mut_s(pilot_loc,:) = [];
fft_mut_l(pilot_loc,:) = [];

% =========================9.采样（去除ifft中补的位数）===========

data_sig = fft_sig(1:carrier_count,:);
data_mut_s = fft_mut_s(1:carrier_count,:);
data_mut_l = fft_mut_l(1:carrier_count,:);

% =========================10.16QAM逆映射===========================

bit_demod_sig = reshape(qamdemod(data_sig,16,'OutputType','bit'),[],1);
bit_demod_mut_s = reshape(qamdemod(data_mut_s,16,'OutputType','bit'),[],1);
bit_demod_mut_l = reshape(qamdemod(data_mut_l,16,'OutputType','bit'),[],1);

% =========================11.误码率===========================
error_bit_sig = sum(bit_demod_sig~=bit_sequence);
error_bit_mut_s = sum(bit_demod_mut_s~=bit_sequence);
error_bit_mut_l = sum(bit_demod_mut_l~=bit_sequence);
error_rate_sig = error_bit_sig/bit_length;
error_rate_mut_s = error_bit_mut_s/bit_length;
error_rate_mut_l = error_bit_mut_l/bit_length;
rate = [error_rate_sig error_rate_mut_s error_rate_mut_l];

rate_comp = [rate_comp; error_rate_sig error_rate_mut_s error_rate_mut_l];

end

figure;
semilogy(SNR_range,rate_comp(:,2),'-sk');
hold on
semilogy(SNR_range,rate_comp(:,3),'-sb');
legend('多径误码率-spline插值','多径误码率-linear插值');
title('LS信道估计不同插值方法下的误码率比较');
xlabel('SNR(dB)');ylabel('Ratio');grid on;


function H_LS = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,int_opt)
% LS channel estimation function
% Inputs:
%       Y         = Frequency-domain received signal
%       Xp        = Pilot signal
%       pilot_loc = Pilot location
%       N         = FFT size
%       Nps       = Pilot spacing
%       int_opt   = 'linear' or 'spline'
% output:
%       H_LS      = LS channel etimate

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

Np=ceil(Nfft/Nps); k=1:Np; LS_est(k) = Y(pilot_loc(k))./Xp(k);  % LS channel estimation
if  lower(int_opt(1))=='l', method='linear'; 
else method='spline';  
end
H_LS = interpolate(LS_est,pilot_loc,Nfft,method); % Linear/Spline interpolation
end

function H_interpolated = interpolate(H_est,pilot_loc,Nfft,method)
% Input:        H_est    = Channel estimate using pilot sequence
%           pilot_loc    = location of pilot sequence
%                Nfft    = FFT size
%              method    = 'linear'/'spline'
% Output: H_interpolated = interpolated channel

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

if pilot_loc(1)>1
    slope = (H_est(2)-H_est(1))/(pilot_loc(2)-pilot_loc(1));
    H_est = [H_est(1)-slope*(pilot_loc(1)-1)  H_est]; pilot_loc = [1 pilot_loc];
end
if pilot_loc(end)<Nfft
    slope = (H_est(end)-H_est(end-1))/(pilot_loc(end)-pilot_loc(end-1));
    H_est = [H_est  H_est(end)+slope*(Nfft-pilot_loc(end))]; pilot_loc = [pilot_loc Nfft];
end
if lower(method(1))=='l', H_interpolated = interp1(pilot_loc,H_est,[1:Nfft]);
else      H_interpolated = interp1(pilot_loc,H_est,[1:Nfft],'spline');
end
end

