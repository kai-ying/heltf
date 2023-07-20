function [fitresult, gof] = createFit(x, abs_I1)
%CREATEFIT(X,ABS_I1)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      X 输入: x
%      Y 输出: abs_I1
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 14-Jul-2023 14:07:18 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( x, abs_I1 );

% 设置 fittype 和选项。
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [3.92520948840252 7.23835918526748e-05 -0.0140286703985415];

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', '无标题拟合 1' );
h = plot( fitresult, xData, yData );
legend( h, 'abs_I1 vs. x', '无标题拟合 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'abs_I1', 'Interpreter', 'none' );
grid on


