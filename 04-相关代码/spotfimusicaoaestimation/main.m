clear;
% sample_csi_trace 是一个90*1的向量矩阵，前30个元素对应第一个rx天线的子载波，后30个元素对应第二个rx天线的子载波，以此推类。
% 将从5300网卡获取的CSI转化为90*1的向量矩阵，赋值到变量sample_csi_trace
% ------------------------------------------------------------------------------------------------------------------------------------------------
% 加载csi数据
sample_csi_traceTmp = load('sample_csi_trace');% 将sample_csi_trace文件中的变量加载到工作区
sample_csi_trace = sample_csi_traceTmp.sample_csi_trace;% 获取sample_csi_trace变量

% 定义常规变量
fc = 5.63e9; % 中心频率 
M = 3;    % rx天线的数量
fs = 40e6; % 信道带宽
c = 3e8;  % 光速
d = 2.6e-2;  % 线性天线阵中相邻天线之间的距离

SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; %WiFi子载波索引，其中CSI是可用的
N = length(SubCarrInd); % 子载波的数量 30
fgap = 312.5e3; % WiFi中连续子载波之间的频率间隔(Hz)
lambda = c/fc;  % 波长=光速/频率
T = 1; % 发射天线的数量
%从变量的定义来看，是一个天线发射，三个天线接收，会得到一个3*30*1的向量矩阵，但在这里看，接收到的csi数据已经被重构成90*1的向量

% ------------------------------------------------------------------------------------------------------------------------------------------------
% MUSIC算法要求在一个网格中估算MUSIC spectrum。paramRange为这个网格捕获的参数
% 下面的示例, MUSIC spectrum 计算101个间距在-25ns~25ns之间的ToF（飞行时间）值. MUSIC spectrum 计算101个-90度到90度的AoA（到达角）值.

paramRange = struct; % 定义paramRange结构体
paramRange.GridPts = [101 101 1]; % 格式为[ToF，AoA， 1]
paramRange.delayRange = [-50 50]*1e-9; % 要考虑的ToF网格的最大值和最小值。[-25ns,25ns]
paramRange.angleRange = 90*[-1 1]; % 要考虑的AoA网格的最大值和最小值。[-90,90]
do_second_iter = 0; % 第二通路？
paramRange.K = floor(M/2)+1; % 与平滑相关的参数。 
paramRange.L = floor(N/2); % 与平滑相关的参数。 
paramRange.T = 1; % ?
paramRange.deltaRange = [0 0];  %变量程

maxRapIters = Inf; % inf为无穷大的意思
useNoise = 0;
paramRange.generateAtot = 2;% 生成Atot ？

% TOF消除代码（Spotfi论文中的算法1）  注：也就是消去不同包不同TOF。   
csi_plot = reshape(sample_csi_trace, N, M);
[PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
csi_plot = csi_plot.*ToMult;
relChannel_noSlope = reshape(csi_plot, N, M, T);
sample_csi_trace_sanitized = relChannel_noSlope(:);

% 到达角估计的MUSIC算法
% aoaEstimateMatrix是（nComps x 5）矩阵，其中nComps是环境中的路径数。第一列为tof的ns，第二列为spotfi论文中定义的aoa度数。
aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
                    T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(2))   
tofEstimate = aoaEstimateMatrix(:,1); % ToF in nanoseconds
aoaEstomate = aoaEstimateMatrix(:,2); % AoA in degrees

% clear;
% % sample CSI trace is a 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and so on.
% % replace sample_csi_trace with CSI from Intel 5300 converted to 90x1 vector
% sample_csi_traceTmp = load('sample_csi_trace'); 
% sample_csi_trace = sample_csi_traceTmp.sample_csi_trace;

% fc = 5.63e9; % center frequency 
% M = 3;    % number of rx antennas 
% fs = 40e6; % channel bandwidth 
% c = 3e8;  % speed of light 
% d = 2.6e-2;  % distance between adjacent antennas in the linear antenna array 
% % dTx = 2.6e-2; 

% SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available WiFi子载波索引，其中CSI是可用的
% N = length(SubCarrInd); % number of subcarriers 
% % subCarrSize = 128;  % total number fo
% fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi 
% lambda = c/fc;  % wavelength
% T = 1; % number of transmitter antennas

% % MUSIC algorithm requires estimating MUSIC spectrum in a grid. paramRange captures parameters for this grid
% % For the following example, MUSIC spectrum is caluclated for 101 ToF (Time of flight) values spaced equally between -25 ns and 25 ns. MUSIC spectrum is calculated for for 101 AoA (Angle of Arrival) values between -90 and 90 degrees.
% paramRange = struct;
% paramRange.GridPts = [101 101 1]; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
% paramRange.delayRange = [-50 50]*1e-9; % lowest and highest values to consider for ToF grid. 
% paramRange.angleRange = 90*[-1 1]; % lowest and values to consider for AoA grid.
% do_second_iter = 0;
% % paramRange.seconditerGridPts = [1 51 21 21];
% paramRange.K = floor(M/2)+1; % parameter related to smoothing.  
% paramRange.L = floor(N/2); % parameter related to smoothing.  
% paramRange.T = 1;
% paramRange.deltaRange = [0 0]; 

% maxRapIters = Inf;
% useNoise = 0;
% paramRange.generateAtot = 2;

% % ToF sanitization code (Algorithm 1 in SpotFi paper)
% csi_plot = reshape(sample_csi_trace, N, M);
% [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
% ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
% csi_plot = csi_plot.*ToMult;
% relChannel_noSlope = reshape(csi_plot, N, M, T);
% sample_csi_trace_sanitized = relChannel_noSlope(:);

% % MUSIC algorithm for estimating angle of arrival
% % aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment. First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
% aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
                    % T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(2))   
% tofEstimate = aoaEstimateMatrix(:,1); % ToF in nanoseconds
% aoaEstomate = aoaEstimateMatrix(:,2); % AoA in degrees