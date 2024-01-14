function x1 = MovingAverageFilter(x,win_sz)
% x:待滑动平均的数据
% win_sz：窗宽

if nargin < 2  %默认窗宽等于7
    win_sz=7;
end

L = length(x); %数据长度
x1 = zeros(L,1); %平均之后的数据

half_win = ceil(win_sz/2);
half_win_ = floor(win_sz/2);
if half_win==half_win_
    half_win = half_win+1;
end

x1(1:half_win) = x(1:half_win);
x1(L-half_win:L) = x(L-half_win:L);

for i = half_win:L-half_win
    k=0;
    for j = i-half_win_:i+half_win_  %对第i个窗里面的数求平均
        k = k+1;
        temp(k) = x(j) ; %临时存储第i个窗的数据
    end
    x1(i) = mean(temp); %第i个窗里面的平均值给第i个数
end
end
