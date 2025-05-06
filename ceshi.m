clear; clc; close all;

%% 系统参数
M = 16;                % 延迟维度
N = 16;                % 多普勒维度
numTrials = 100;      % 蒙特卡洛仿真次数（增加以提高统计准确性）
QPSK = comm.QPSKModulator('PhaseOffset',pi/4, 'BitInput',true);
QPSKDemod = comm.QPSKDemodulator('PhaseOffset',pi/4,'BitOutput',true,'OutputDataType','double');
maxSpeed = 500;        % 最大速度(km/h)
fc = 4e9;              % 载波频率
df = 15e3;             % 子载波间隔
fs = M*df;             % 采样率
padLen = 15;           % CP长度

%% 信道参数（5径信道，与论文表II一致）
tau_max = [0, 15, 30, 45, 60]*1e-9; % 时延（ns）
gain = db2mag([1.0, -1.804, -3.565, -5.376, -8.860]); % 增益（dB转线性）
max_doppler = (maxSpeed*1000/3600)*fc/3e8; % 最大多普勒
doppler = max_doppler*cos(2*pi*rand(1,5)); % Jakes模型

%% 导频参数（按论文优化值）
sigma_p_sq = 0.3;      % 导频功率
sigma_d_sq = 1 - sigma_p_sq; % 数据功率
pilot = sqrt(sigma_p_sq)*exp(1i*pi/4); % QPSK导频
dataConst = qammod(0:3,4,'UnitAveragePower',true); % 4-QAM

%% 初始化BER存储
SNR_dB = 0:2:20;       % SNR范围
BER_SPNI = zeros(size(SNR_dB));
BER_SPI = zeros(size(SNR_dB));
BER_EP = zeros(size(SNR_dB));
BER_OTFS = zeros(size(SNR_dB));

%% 主仿真循环
for snrIdx = 1:length(SNR_dB)
    snr = SNR_dB(snrIdx);
    fprintf('Processing SNR = %d dB...\n', snr);
    
    % 噪声方差（考虑符号功率归一化）
    N0 = 10^(-snr/10);
    berSPNI = 0; berSPI = 0; berEP = 0; berOTFS = 0;
    
    parfor trial = 1:numTrials % 并行加速
        %% 生成数据符号
        bits = randi([0 1], 2*M*N, 1);
        data = step(QPSK, bits);
        dataGrid = reshape(data, M, N);
        
        %% ------------------ SP方案 ------------------
        % 叠加导频
        pilotGrid = zeros(M,N);
        pilotGrid(1, N/2+1) = pilot; % 中心导频
        txGrid_SP = sqrt(sigma_d_sq)*dataGrid + pilotGrid;
        
        % OTFS调制
        txSig_SP = helperOTFSmod(txGrid_SP, padLen, 'CP');
        
        % 多普勒信道
        rxSig_SP = dopplerChannel(txSig_SP, fs, struct(...
            'pathDelays',round(tau_max*fs),...
            'pathGains',gain,...
            'pathDopplerFreqs',doppler));
        rxSig_SP = awgn(rxSig_SP, snr, 'measured');
        
        % SP-NI信道估计
        H_est_NI = SP_ChannelEstimation(rxSig_SP, pilotGrid, sigma_d_sq, sigma_p_sq, N0, M, N, padLen);
        
        % SP-NI数据检测
        dataEst_NI = SP_DataDetection(rxSig_SP, H_est_NI, pilotGrid, sigma_p_sq, sigma_d_sq, N0, 'NI', M, N, padLen);
        berSPNI = berSPNI + biterr(bits, step(QPSKDemod, dataEst_NI(:)));
        
        % SP-I迭代检测
        dataEst_I = SP_DataDetection(rxSig_SP, H_est_NI, pilotGrid, sigma_p_sq, sigma_d_sq, N0, 'I', M, N, padLen);
        berSPI = berSPI + biterr(bits, step(QPSKDemod, dataEst_I(:)));
        
        %% ------------------ EP方案 ------------------
        [txSig_EP, dataGrid_EP] = EP_Transmit(data, M, N, padLen);
        rxSig_EP = awgn(dopplerChannel(txSig_EP, fs, struct(...
            'pathDelays',round(tau_max*fs),...
            'pathGains',gain,...
            'pathDopplerFreqs',doppler)), snr, 'measured');
        H_est_EP = EP_ChannelEstimation(rxSig_EP, M, N, padLen);
        dataEst_EP = EP_DataDetection(rxSig_EP, H_est_EP, dataGrid_EP, N0);
        berEP = berEP + biterr(bits, step(QPSKDemod, dataEst_EP(:)));
        
        %% ------------------ 理想信道估计 ------------------
        rxSig_OTFS = awgn(dopplerChannel(txSig_SP, fs, struct(...
            'pathDelays',round(tau_max*fs),...
            'pathGains',gain,...
            'pathDopplerFreqs',doppler)), snr, 'measured');
        dataEst_OTFS = Ideal_Detection(rxSig_OTFS, M, N, padLen, tau_max, gain, doppler);
        berOTFS = berOTFS + biterr(bits, step(QPSKDemod, dataEst_OTFS(:)));
    end
    
    % 计算平均BER
    BER_SPNI(snrIdx) = berSPNI/(numTrials*length(bits));
    BER_SPI(snrIdx) = berSPI/(numTrials*length(bits));
    BER_EP(snrIdx) = berEP/(numTrials*length(bits));
    BER_OTFS(snrIdx) = berOTFS/(numTrials*length(bits));
end

%% 绘图
figure;
semilogy(SNR_dB, BER_SPNI, 'b-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB, BER_SPI, 'r-s', 'LineWidth', 2);
semilogy(SNR_dB, BER_EP, 'g-d', 'LineWidth', 2);
semilogy(SNR_dB, BER_OTFS, 'k-^', 'LineWidth', 2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
legend('SP-NI', 'SP-I', 'EP', 'Ideal OTFS');
title('OTFS BER性能比较（500 km/h）');
axis([0 20 1e-5 1]);

%% ================== 辅助函数 ==================
%% ================== SP信道估计函数（公式20） ==================
function H_est = SP_ChannelEstimation(rxSig, pilotGrid, sigma_d, sigma_p, N0, M, N, padLen)
    % 论文式(20) MMSE信道估计
    Y_dd = helperOTFSdemod(rxSig, M, padLen, 0, 'CP');
    y = Y_dd(:);
    
    % 构建导频矩阵Omega_p（论文式13）
    Q = 5; % 5径信道
    Omega_p = zeros(M*N, Q);
    pilotVec = pilotGrid(:);
    for i = 1:Q
        % 计算Gamma_i矩阵作用后的导频（简化为循环移位）
        shift_l = i-1; % 延迟移位
        shift_k = i-1; % 多普勒移位
        shifted_pilot = circshift(pilotGrid, [shift_l, shift_k]);
        Omega_p(:,i) = shifted_pilot(:);
    end
    
    % 计算噪声+干扰协方差矩阵（论文式19）
    C_wd = (sum(abs(pilotVec).^2))*sigma_d + N0;
    C_wd = C_wd * eye(M*N);
    
    % 信道协方差矩阵（已知路径增益方差）
    C_h = diag([1.0, 0.7, 0.5, 0.3, 0.2]); % 论文表II
    
    % MMSE估计（论文式20）
    H_est = (Omega_p' * (C_wd \ Omega_p) + inv(C_h)) \ (Omega_p' * (C_wd \ y));
end

%% ================== SP数据检测函数（NI/I模式） ==================
function dataEst = SP_DataDetection(rxSig, H_est, pilotGrid, sigma_p, sigma_d, N0, mode, M, N, padLen)
    % 输入参数处理
    Y_dd = helperOTFSdemod(rxSig, M, padLen, 0, 'CP');
    y = Y_dd(:);
    
    % 去除导频分量（论文式23）
    y_data = y - pilotGrid(:).*H_est;
    
    % 构建有效信道矩阵（论文式9）
    H_eff = construct_Heff(H_est, M, N);
    
    if strcmpi(mode, 'NI')
        % 非迭代MMSE检测
        dataEst = (H_eff'*H_eff + (N0/sigma_d)*eye(M*N)) \ (H_eff'*y_data);
    else
        % 迭代检测（算法2）
        maxIter = 5;   % 论文参数
        dataEst = y_data; % 初始化
        prevEst = dataEst;
        for iter = 1:maxIter
            % 数据辅助信道估计（论文式41）
            H_est = SP_IterativeChannelEstimation(y, dataEst, pilotGrid, sigma_d, sigma_p, N0, M, N);
            
            % 消息传递检测（算法1）
            dataEst = MessagePassingDetection(y_data, H_est, sigma_p, sigma_d, N0, M, N);
            
            % 早停检查
            if norm(dataEst - prevEst) < 1e-6
                break;
            end
            prevEst = dataEst;
        end
    end
    dataEst = reshape(dataEst, M, N);
end

function H_eff = construct_Heff(H_est, M, N)
    % 构建稀疏有效信道矩阵（论文式9）
    Q = length(H_est);
    H_eff = sparse(M*N, M*N);
    for i = 1:Q
        H_eff = H_eff + H_est(i)*circshift(speye(M*N), [i-1, i-1]);
    end
end

%% ================== 消息传递检测（算法1） ==================
function x_hat = MessagePassingDetection(y, H_est, sigma_p, sigma_d, N0, M, N)
    % 参数设置
    maxIter = 10;   % 最大迭代次数
    damping = 0.5;  % 阻尼因子
    Q = length(H_est);
    MN = M*N;
    constellation = qammod(0:3,4,'UnitAveragePower',true); % 4-QAM
    S = length(constellation);
    
    % 初始化变量节点
    [varNodes, obsNodes] = initFactorGraph(H_est, M, N);
    
    % 迭代消息传递
    for iter = 1:maxIter
        % 观测节点到变量节点
        for a = 1:MN
            neighbors = obsNodes(a).connections;
            for b = neighbors
                % 计算干扰项均值和方差（论文式32-33）
                mu = 0;
                sigma2 = 0;
                for c = setdiff(neighbors, b)
                    symProb = varNodes(c).pmf;
                    mu = mu + sum(symProb .* constellation) * H_est(a,c);
                    sigma2 = sigma2 + sum(symProb .* abs(constellation).^2) * abs(H_est(a,c))^2;
                end
                sigma2 = sigma2 - abs(mu)^2 + sigma_p^2*norm(H_est)^2 + N0;
                
                % 更新观测节点消息
                obsNodes(a).mu(b) = mu;
                obsNodes(a).sigma2(b) = sigma2;
            end
        end
        
        % 变量节点到观测节点
        for b = 1:MN
            neighbors = varNodes(b).connections;
            for a = neighbors
                % 计算概率分布（论文式34-35）
                beta = zeros(S,1);
                for s = 1:S
                    beta(s) = exp(-abs(y(a) - obsNodes(a).mu(b) - H_est(a,b)*constellation(s))^2 / obsNodes(a).sigma2(b));
                end
                prob = beta / sum(beta);
                
                % 阻尼更新
                varNodes(b).pmf = damping*prob + (1-damping)*varNodes(b).pmf;
            end
        end
    end
    
    % 硬判决输出
    x_hat = zeros(MN,1);
    for b = 1:MN
        [~, idx] = max(varNodes(b).pmf);
        x_hat(b) = constellation(idx);
    end
end

%% ================== EP方案函数 ==================
function [txSig, dataGrid] = EP_Transmit(data, M, N, padLen)
    % 嵌入式导频传输（论文图1(a)）
    dataGrid = reshape(data, M, N);
    
    % 插入导频和保护间隔（保护带尺寸按论文设置）
    pilotPos = [M/2, N/2]; % 中心位置
    guardLen = 4;          % 保护间隔
    
    % 创建EP帧结构
    dataGrid(pilotPos(1)-guardLen:pilotPos(1)+guardLen, ...
             pilotPos(2)-guardLen:pilotPos(2)+guardLen) = 0; 
    dataGrid(pilotPos(1), pilotPos(2)) = 1; % 导频符号
    
    % OTFS调制
    txSig = helperOTFSmod(dataGrid, padLen, 'CP');
end

function H_est = EP_ChannelEstimation(rxSig, M, N, padLen)
    % EP信道估计（阈值法，论文式6）
    Y_dd = helperOTFSdemod(rxSig, M, padLen, 0, 'CP');
    
    % 阈值检测路径
    threshold = 3*std(Y_dd(:));
    [lIdx, kIdx] = find(abs(Y_dd) > threshold);
    
    % 提取路径增益
    H_est = zeros(M, N);
    for i = 1:length(lIdx)
        H_est(lIdx(i), kIdx(i)) = Y_dd(lIdx(i), kIdx(i));
    end
end

function dataEst = EP_DataDetection(rxSig, H_est, dataGrid, N0)
    % EP数据检测（MMSE）
    Y_dd = helperOTFSdemod(rxSig, size(dataGrid,1), 0, 'CP');
    dataEst = (H_est'*H_est + N0*eye(numel(dataGrid))) \ (H_est'*Y_dd(:));
end

%% ================== 理想检测函数 ==================
function dataEst = Ideal_Detection(rxSig, M, N, padLen, tau, gain, doppler)
    % 已知完美信道信息的MMSE检测
    H_perfect = getPerfectChannel(M, N, tau, gain, doppler);
    Y_dd = helperOTFSdemod(rxSig, M, padLen, 0, 'CP');
    dataEst = (H_perfect'*H_perfect + 1e-6*eye(M*N)) \ (H_perfect'*Y_dd(:));
end

function H = getPerfectChannel(M, N, tau, gain, doppler)
    % 生成完美信道矩阵
    Q = length(tau);
    H = zeros(M,N);
    for i = 1:Q
        l = round(tau(i)*M/(1/(M*15e3))); % 延迟采样点
        k = round(doppler(i)*N/(1/(N/15e3))); % 多普勒采样点
        H(mod(l,M)+1, mod(k,N)+1) = gain(i);
    end
end