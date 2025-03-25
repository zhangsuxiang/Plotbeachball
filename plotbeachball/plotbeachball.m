function plotbeachball(str, dip, rake)
% PLOTBEACHBALL 绘制震源机制解（海滩球图）
%
% 输入参数（单位均为度）:
%   str   - 断层面走向
%   dip   - 断层面倾角
%   rake  - 断层面滑动角 
%
% 说明:
%   该函数设置了默认参数用于绘图：
%       siz      = 1;     % 绘图大小
%       n        = 20;    % 绘图分辨率（自动限制在5-50之间）
%       clr      = 'r';   % 填充颜色
%       (lat,lon)= (0,0); % 绘图中心位置
%       fillcode = 1;     % 采用填充方式
%
%   注意：若滑动角大于180°，请先将其转换为负角。

% 默认参数设置
siz = 1;      % 绘图大小
n = 20;       % 分辨率（建议值在5-50之间）
clr = 'r';    % 填充颜色
lat = 0; lon = 0;  % 事件位置
fillcode = 1; % 是否采用填充

if n < 5, n = 5; elseif n > 50, n = 50; end
d2r = pi/180; 
pi2 = 2*pi;

% 将输入角度转换为弧度
str_rad = str * d2r; 
dip_rad = dip * d2r; 
rake_rad = rake * d2r;

% 计算滑动矢量 A 和断层面法向量 N（单位方向矢量）
A = [cos(rake_rad)*cos(str_rad) + sin(rake_rad)*cos(dip_rad)*sin(str_rad), ...
     cos(rake_rad)*sin(str_rad) - sin(rake_rad)*cos(dip_rad)*cos(str_rad), ...
     -sin(rake_rad)*sin(dip_rad)];
N = [-sin(str_rad)*sin(dip_rad), cos(str_rad)*sin(dip_rad), -cos(dip_rad)];

% 利用外部函数 an2dsr_wan 得到第二个节面的参数（单位：度）
% 请确保你有该函数，或自行实现相应转换
[str2, dip2, rake2] = an2dsr_wan(N, A);

% --------------------------
% 计算第一个断层面（基于输入的 str, dip, rake）
% --------------------------
rak = 0:-pi/n:-pi;
cosih = -sin(dip_rad)*sin(rak);
ih = acos(cosih);
cosdet = sqrt(1 - cosih.^2);
fai = acos(cos(rak)./cosdet);
str1 = str_rad + fai;  % 在走向基础上增加夹角
xs1 = siz * sqrt(2) * sin(ih/2) .* sin(str1);
ys1 = siz * sqrt(2) * sin(ih/2) .* cos(str1);

% --------------------------
% 计算第二个断层面（基于计算得到的 str2, dip2）
% --------------------------
str2_rad = str2 * d2r; 
dip2_rad = dip2 * d2r;
cosih = -sin(dip2_rad)*sin(rak);
ih = acos(cosih);
cosdet = sqrt(1 - cosih.^2);
fai = acos(cos(rak)./cosdet);
str21 = str2_rad + fai;
xs2 = siz * sqrt(2) * sin(ih/2) .* sin(str21);
ys2 = siz * sqrt(2) * sin(ih/2) .* cos(str21);

% --------------------------
% 调整弧线起止角度
% --------------------------
str1_adj = str_rad + pi;
if str1_adj > pi2
    str1_adj = str1_adj - pi2;
end
d = str2_rad;
d1 = d + pi;
if (str1_adj - d) > pi
    d = d - pi2;
elseif (str1_adj - d) >= pi
    str1_adj = str1_adj - pi2;
end
if (d1 - str_rad) > pi
    d1 = d1 - pi2;
elseif (d1 - str_rad) >= pi
    str_rad = str_rad - pi2;
end

st1 = linspace(str1_adj, d, n);
st2 = linspace(d1, str_rad, n);

% 连接两个断层面小弧和计算完整海滩球轮廓坐标
p = [xs1, siz*sin(st1), xs2, siz*sin(st2)] + lon;
q = [ys1, siz*cos(st1), ys2, siz*cos(st2)] + lat;

% --------------------------
% 绘制海滩球图
% --------------------------
n_circle = n * 2;
theta_circle = linspace(0, 2*pi, n_circle);
x_circle = siz * cos(theta_circle);
y_circle = siz * sin(theta_circle);

figure;
if fillcode
    % 根据滑动角判断填充哪一部分
    if rake > 0
        % 若滑动角为正（逆冲分量），先填充白色，再用指定颜色填充断层面区域
        fill(x_circle, y_circle, 'w');
        hold on;
        fill(p, q, clr);
    else
        % 若滑动角为负（正断层分量），相反填充
        fill(x_circle, y_circle, clr);
        hold on;
        fill(p, q, 'w');
    end
else
    plot(x_circle, y_circle, 'k');
    hold on;
    plot(xs1, ys1, clr, xs2, ys2, clr);
end
axis equal;
axis off;
title('震源机制解');
return