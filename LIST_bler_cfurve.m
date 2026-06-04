clear; clc; close all;
tic;

%% Parameters
N = 256;
K = 128;
R = K/N;
L1 = 1; % the base size list to compare
L2 = 4; % the compared size list
SNR_dB_range = 4.4;
SNR_dB_range = -1.01:0.5:1.02;

num_trials_per_snr = 1000;

IS_SC = false;
L = 1;
GA_CRC_Length = 0;

bler1 = zeros(size(SNR_dB_range)); % AWGN
bler2 = zeros(size(SNR_dB_range)); % 

ber1 = zeros(size(SNR_dB_range));
ber2 = zeros(size(SNR_dB_range));

%% ===== Loop over SNR values =====
for SNR = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(SNR);
    EbN0  = SNR_dB - 10*log10(R);
    EbN0_lin = 10^(EbN0/10);

    % Correct sigma
    Sigma = sqrt(1 / (2 * 10^(SNR_dB/10)));

    % BSC crossover probability
    p = qfunc(sqrt(2*R*EbN0_lin));

    total_errors1 = 0;
    total_errors2 = 0;
    block_errors1 = 0;
    block_errors2 = 0;
    total_bits = 0;

    fprintf('SNR = %.2f dB\n', SNR_dB);
    fprintf('Eb/N0 = %.2f dB\n', EbN0);

    %% ===== Load frozen set per SNR (your logic preserved) =====
    filename = sprintf('SNR_4.4_dB_BSC_HARD_256_1.mat');
    %filename = sprintf('SNR_%0.1f_dB_AWGN_16.mat', SNR_dB);

    data = load(filename);

    % ===== ADJUST THIS LINE BASED ON YOUR FILE =====
    frozen_bits_indicator = data.frozen_indicator;
    %% ===== Monte Carlo =====
    for trial = 1:num_trials_per_snr

        if mod(trial,500)==0
            fprintf('Trial %d/%d\n', trial, num_trials_per_snr);
        end

        %% Generate valid u (respect frozen bits)
        u = zeros(1,N);
        info_bits = randi([0 1],1,K);
        u(frozen_bits_indicator==0) = info_bits;

        %% Encode
        x = polar_encoder(u);
        s = 1 - 2*x;

        %% ===== noise shared channel =====
        noise = Sigma * randn(1,N);
        Y = s + noise;
        Y = Y(:)';


        %% ===== LLRs =====
        Lambda = BSC_BPSK_LLR(Y, Sigma);
        Lambda = Lambda(:)'; % Force 1xN

        %% ===== Decode =====
        if IS_SC
            [Estimated_U,~,~] = SC_Decoder(Lambda, u, N, frozen_bits_indicator, false);
            [Estimated_U,~,~] = SC_Decoder(Lambda, u, N, frozen_bits_indicator, false);
        else
            % AWGN
            [Estimated_U1,~] = gArikan_BPSK_SCL_Decoder(L1,GA_CRC_Length,Y,u,Sigma,SNR_dB,frozen_bits_indicator,'BSC');
            
            % BSC
            [Estimated_U2,~] = gArikan_BPSK_SCL_Decoder(L2,GA_CRC_Length,Y,u,Sigma,SNR_dB,frozen_bits_indicator,'BSC');
        end

        %% Errors
        num_errors1 = sum(Estimated_U1 ~= u);
        num_errors2 = sum(Estimated_U2 ~= u);

        total_errors1 = total_errors1 + num_errors1;
        total_errors2 = total_errors2 + num_errors2;
        total_bits = total_bits + N;

        if any(Estimated_U1 ~= u)
            block_errors1 = block_errors1 + 1;
        end

        if any(Estimated_U2 ~= u)
            block_errors2 = block_errors2 + 1;
        end

    end

    %% Store results
    bler1(SNR) = block_errors1 / num_trials_per_snr;
    bler2(SNR) = block_errors2 / num_trials_per_snr;

    fprintf('List 1 BLER = %.3e | List 2 BLER= %.3e\n\n', bler1(SNR), bler2(SNR));

    ber1(SNR) = total_errors1 / total_bits;
    ber2(SNR) = total_errors2 / total_bits;

end

%% ===== Plot =====
figure;
semilogy(SNR_dB_range, bler1, 'b-o','LineWidth',2); hold on;
semilogy(SNR_dB_range, bler2, 'r-s','LineWidth',2);
grid on;
xlabel('SNR [dB]');
ylabel('BLER');
legend('L = 1','L = 4');
title(sprintf('Polar Code BLER Comparison (N=%d, K=%d)', N, K));

figure;
semilogy(SNR_dB_range, ber1, 'b-o','LineWidth',2); hold on;
semilogy(SNR_dB_range, ber2, 'r-s','LineWidth',2);
grid on;
xlabel('SNR [dB]');
ylabel('BER');
legend('LIST 1 BER','LIST 2 BER');
title(sprintf('Polar Code BER Comparison (N=%d, K=%d)', N, K));

%% Save
savefig('LIST_SIZE_BLER_COMPARISON.fig');

elapsedTime = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);