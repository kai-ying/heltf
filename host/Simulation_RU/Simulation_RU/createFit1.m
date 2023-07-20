function [fitresult, gof] = createFit1(x, abs_I1)
%CREATEFIT1(X,ABS_I1)
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

%  由 MATLAB 于 14-Jul-2023 14:12:25 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( x, abs_I1 );

% 设置 fittype 和选项。
ft = fittype( 'smoothingspline' );

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft );

% 为绘图创建一个图窗。
figure( 'Name', '无标题拟合 1' );

% 绘制数据拟合图。
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData );
legend( h, 'abs_I1 vs. x', '无标题拟合 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'abs_I1', 'Interpreter', 'none' );
grid on

% 绘制残差图。
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, '无标题拟合 1 - 残差', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'abs_I1', 'Interpreter', 'none' );
grid on


