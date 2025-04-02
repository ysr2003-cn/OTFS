clear;
clc;
close all;

% Set Monte Carlo simulation iterations
numTrials = 100;  % Increase simulation iterations for more accurate results

% Initialize parameters
M = 128;          % Number of subcarriers
N = 50;          % Number of symbols per frame
df = 15e3;       % Frequency spacing
fc = 4e9;        % Carrier frequency in Hz
padLen = 20;     % Padding length, greater than channel delay spread
padType = 'ZP';  % Use zero padding (ZP)

% Pilot generation and grid filling
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % Fill only one grid

% OTFS modulation
txOut = helperOTFSmod(Pdd, padLen, padType);

% Configure channel parameters, consider fractional Doppler values
chanParams.pathDelays      = [0  5   8];    % Path delays (in samples)
chanParams.pathGains       = [1  0.7 0.5];  % Path gains
chanParams.pathDopplers    = [0 -3.5 5.25]; % Doppler indices, can be fractional

fsamp = M*df;            % Sampling frequency
Meff = M + padLen;       % Number of samples per OTFS sub-symbol
numSamps = Meff * N;     % Number of samples per OTFS symbol
T = ((M+padLen)/(M*df)); % Symbol duration (seconds)

% Calculate actual Doppler frequencies, including fractional parts
chanParams.pathDopplerFreqs = chanParams.pathDopplers / (N*T); % In Hz

% Define range of pilot SNR
pilotSNRdB = [30 35 40];  % Pilot SNR values
dataSNRdB = [10 15 18];   % Data SNR values
Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));

nmse = zeros(length(pilotSNRdB), 1);  % Store NMSE for each pilot SNR value

% Original channel response for NMSE calculation
dopplerOut_pilot = dopplerChannel(txOut, fsamp, chanParams);
chOut_pilot = awgn(dopplerOut_pilot, 0, 'measured');  % Use fixed pilot SNR to generate true channel
rxIn_pilot = chOut_pilot(1:numSamps);
Ydd_pilot = helperOTFSdemod(rxIn_pilot, M, padLen, 0, padType);
H_true = Ydd_pilot * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2);

% for idx = 1:length(pilotSNRdB)
%     pilotSNR = pilotSNRdB(idx);
%     n0_pilot = Es / (10^(pilotSNR / 10));
% 
%     nmse_trials = zeros(1, numTrials);
%     for trial = 1:numTrials
%         dopplerOut_pilot = dopplerChannel(txOut, fsamp, chanParams);
%         chOut_pilot = awgn(dopplerOut_pilot, pilotSNR, 'measured');
%         rxIn_pilot = chOut_pilot(1:numSamps);
%         Ydd_pilot = helperOTFSdemod(rxIn_pilot, M, padLen, 0, padType);
%         H_est = Ydd_pilot * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2 + n0_pilot);
%         nmse_trials(trial) = sum(sum(abs(H_est - H_true).^2)) / sum(sum(abs(H_true).^2));
%     end
%     nmse(idx) = mean(nmse_trials); 
% end
% 
% % Plot NMSE curve
% figure;
% plot(pilotSNRdB, 10*log10(nmse), '-o');
% grid on;
% xlabel('Pilot SNR (dB)');
% ylabel('NMSE (dB)');
% title('NMSE of Channel Estimation vs. Pilot SNR');

% Extract channel path information
[lp, vp] = find(abs(H_true) >= 0.05);
chanEst.pathGains = diag(H_true(lp, vp));   
chanEst.pathDelays = lp - 1;            
chanEst.pathDopplers = vp - pilotBin;   

% Generate data
Xgrid = zeros(M, N);
Xdata = randi([0, 1], 2*M, N);
Xgrid(1:M, :) = pskmod(Xdata, 4, pi/4, InputType="bit");

% OTFS modulation
txOut_data = helperOTFSmod(Xgrid, padLen, padType);

% Store BER results for different pilot SNR values
berOTFS = zeros(length(dataSNRdB), length(pilotSNRdB));  % Store BER for OTFS

for pIdx = 1:length(pilotSNRdB)
    pilotSNR = pilotSNRdB(pIdx);
    n0_pilot = Es / (10^(pilotSNR / 10));

    % Recalculate channel estimation for the given pilot SNR
    dopplerOut_pilot = dopplerChannel(txOut, fsamp, chanParams);
    chOut_pilot = awgn(dopplerOut_pilot, pilotSNR, 'measured');
    rxIn_pilot = chOut_pilot(1:numSamps);
    Ydd_pilot = helperOTFSdemod(rxIn_pilot, M, padLen, 0, padType);
    H_est = Ydd_pilot * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2 + n0_pilot);

    % Update channel path information for data estimation
    [lp, vp] = find(abs(H_est) >= 0.05);
    chanEst.pathGains = diag(H_est(lp, vp));   
    chanEst.pathDelays = lp - 1;            
    chanEst.pathDopplers = vp - pilotBin;

    for dIdx = 1:length(dataSNRdB)
        dataSNR = dataSNRdB(dIdx);
        n0_data = Es / (10^(dataSNR / 10));

        berOTFS_trials = zeros(1, numTrials);
        
        for trial = 1:numTrials
            % OTFS signal through channel and noise
            dopplerOut_otfs = dopplerChannel(txOut_data, fsamp, chanParams);
            chOut_otfs = awgn(dopplerOut_otfs, dataSNR, 'measured');
            
            % OTFS reception and demodulation
            rxWindow_otfs = chOut_otfs(1:numSamps);
            G = getG(M, N, chanEst, padLen, padType);
            y_otfs = ((G'*G) + n0_data * eye(Meff * N)) \ (G' * rxWindow_otfs); % LMMSE
            Xhat_otfs = helperOTFSdemod(y_otfs, M, padLen, 0, padType); % OTFS demodulation
            XhatDataOTFS = pskdemod(Xhat_otfs, 4, pi/4, OutputType="bit", OutputDataType="logical");
            [~, berOTFS_trials(trial)] = biterr(Xdata, XhatDataOTFS);
        end
        
        % Calculate average BER
        berOTFS(dIdx, pIdx) = mean(berOTFS_trials);
    end
end

% Plot BER curves for different pilot SNR values
figure;
hold on;
markers = {'-s', '-o', '-x'};
for pIdx = 1:length(pilotSNRdB)
    semilogy(dataSNRdB, berOTFS(:, pIdx), markers{pIdx}, 'DisplayName', ['Pilot SNR = ', num2str(pilotSNRdB(pIdx)), ' dB']);
end
grid on;
xlabel('Data SNR (dB)');
ylabel('BER');
set(gca,'YScale','log','YLim',[1e-5,1e-1]);
title(['BER Performance of OTFS (', num2str(numTrials), ' Monte Carlo Simulations)']);
legend show;
hold off;