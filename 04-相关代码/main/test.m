clear;
csi_trace = read_bf_file('data/Test1/64_HT20_test1.dat');
% --------------------------------------------------------------------------------------------
% 变量定义
fc = 5.63e9; 
M = 3;    
fs = 40e6; 
c = 3e8;  
d = 2.6e-2;  

SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; %WiFi瀛愯浇娉㈢储寮曪紝鍏朵腑CSI鏄彲鐢ㄧ殑
N = length(SubCarrInd); % 瀛愯浇娉㈢殑鏁伴噺 30
fgap = 312.5e3; % WiFi涓繛缁瓙杞芥尝涔嬮棿鐨勯鐜囬棿闅?(Hz)
lambda = c/fc;  % 娉㈤暱=鍏夐??/棰戠巼
T = 1; % 鍙戝皠澶╃嚎鐨勬暟閲?

paramRange = struct; % 瀹氫箟paramRange缁撴瀯浣?
paramRange.GridPts = [101 101 1]; % 鏍煎紡涓篬ToF锛孉oA锛? 1]
paramRange.delayRange = [-50 50]*1e-9; % 瑕佽?冭檻鐨凾oF缃戞牸鐨勬渶澶у?煎拰鏈?灏忓?笺?俒-25ns,25ns]
paramRange.angleRange = 90*[-1 1]; % 瑕佽?冭檻鐨凙oA缃戞牸鐨勬渶澶у?煎拰鏈?灏忓?笺?俒-90,90]
do_second_iter = 0; % 绗簩閫氳矾锛?
paramRange.K = floor(M/2)+1; % 涓庡钩婊戠浉鍏崇殑鍙傛暟銆? 
paramRange.L = floor(N/2); % 涓庡钩婊戠浉鍏崇殑鍙傛暟銆? 
paramRange.T = 1; % ?
paramRange.deltaRange = [0 0];  %鍙橀噺绋?

maxRapIters = Inf; % inf涓烘棤绌峰ぇ鐨勬剰鎬?
useNoise = 0;
paramRange.generateAtot = 2;% 鐢熸垚Atot 锛?
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
