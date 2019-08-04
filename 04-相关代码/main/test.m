clear;
csi_trace = read_bf_file('data/Test1/64_HT20_test1.dat'); % 读取CSI文件
% --------------------------------------------------------------------------------------------
% 变量定义
fc = 5.63e9; % 中心频率
M = 3;    % rx天线数量
fs = 40e6; % 信道带宽
c = 3e8;  %  光速
d = 2.6e-2;  % 线性天线阵中相邻天线之间的距离

SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi子载波索引，其中CSI是可用的
N = length(SubCarrInd); % 子载波数
fgap = 312.5e3; % WiFi中连续子载波之间的频率间隔(Hz)
lambda = c/fc;  % 波长
T = 1; % 发射天线数量

paramRange = struct; % 创建一个结构体
paramRange.GridPts = [101 101 1]; % 格式为[ToF格点数(飞行时间)，到达角格点数(AoA)， 1]
paramRange.delayRange = [-50 50]*1e-9; % 要考虑的ToF网格的最大值和最小值。
paramRange.angleRange = 90*[-1 1]; % 为AoA网格考虑的最小值和值。
do_second_iter = 0; 
paramRange.K = floor(M/2)+1; % 与平滑相关的参数。
paramRange.L = floor(N/2); % 与平滑相关的参数。
paramRange.T = 1; 
paramRange.deltaRange = [0 0];  

maxRapIters = Inf;
useNoise = 0;
paramRange.generateAtot = 2;
AoA = [];%存放aoa值的数组
% --------------------------------------------------------------------------------------------
% 循环处理csi数据
for i=1:10 %这里是取的数据包的个数
    % --------提取数据包----------------
    csi_entry = csi_trace{i};
    csi = get_scaled_csi(csi_entry); %提取csi矩阵   
    % --------TOF----------------
    csi_plot = reshape(csi, N, M);% 转换为30*3的矩阵 
    [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
    ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
    csi_plot = csi_plot.*ToMult;
    relChannel_noSlope = reshape(csi_plot, N, M, T);
    sample_csi_trace_sanitized = relChannel_noSlope(:);
    % --------AOA----------------
    aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
                    T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(2)) ;  
    tofEstimate = aoaEstimateMatrix(:,1); % ToF in nanoseconds
    aoaEstomate = aoaEstimateMatrix(:,2); % AoA in degrees
    % --------存放进入A0A数组----------------
    aoaEstomate = aoaEstomate';
    AoA = [AoA;aoaEstomate];
end
