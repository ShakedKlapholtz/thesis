function SNR_CHECK()
tic;                  % start time
%% Parameters
N = 256;
K = 96;
R = K/N;                          % code rate
%EbN0_range = 1:0.5:2.5;
%SNR_dB_range = EbN0_range + 10*log10(R);
IS_SC = false ;
GA_CRC_Length = 0 ;
seed = 42;
rng(seed);     % set the seed
SNR_dB = 4.4;
EbN0 = SNR_dB - 10*log10(R);
target_BLER = 0.01;
min_errors = 50;
num_trials_per_snr = min_errors/target_BLER;  % Increase for more accuracy
num_simulations = 100;
filename = sprintf('EbN0_-0.85__BSC.mat');
data = load(filename);
frozen_indicator = data.frozen_indicator;
ber_vec = data.ber_vec;
bit_order = data.bit_order;
frozen_mask = (frozen_indicator == 1);
%[~, order] = sort(ber_vec, 'descend');    % order = מיקום הביטים בסדר ההקפאה
ber_sorted = sort(ber_vec);
%ber_sorted = ber_vec(order);              % BER במיון יורד
%plot_key_length_probability_from_sorted_ber(ber_sorted);
pct_diff_list = []
for sim = 1:num_simulations
    if mod(sim, 5) == 0
        fprintf('Trial %d/%d\n', sim, num_simulations);
    end
    frozen_indicator = zeros(1, N);   % vector of length N
    U_Vec = randi([0, 1], 1, N);
    u = U_Vec;
    % --- Encode
    x = polar_encoder(u);
    bpsk_x = 1 - 2 * x;

    % --- Noise
    SNR_linear = 10^(SNR_dB / 10);
    Sigma = sqrt(1 / (2 * SNR_linear));
    noise = Sigma * randn(1, N);
    Y = bpsk_x + noise;

    % --- LLR
    %Lambda = AWGN_BPSK_LLR(Y, Sigma);
    Lambda = BSC_BPSK_LLR_HARD(Y,EbN0);
    z = sign(Lambda);
    ber_vec = (bpsk_x ~= z);      % 1 = error, 0 = correct
    ber_pct = mean(ber_vec) * 100;
    %display(mean(ber_pct));
    temp = mean(ber_pct);
    pct_diff_list = [pct_diff_list, temp];
end 
display(mean(pct_diff_list));

end