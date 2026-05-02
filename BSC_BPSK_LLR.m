function lambda_BSC = BSC_BPSK_LLR(Y,SNR_dB)
            
% Inputs:
%   symbols    - Vector of received BPSK symbols (+1/-1)
%   EbN0_dB    - Energy per bit to noise ratio in dB
%
% Output:
%   LLR        - Vector of log-likelihood ratios

    % Convert Eb/N0 from dB to linear
    SNR_dB_Lin = 10.^(SNR_dB/10);
    
    % Compute error probability
    Pe = qfunc(sqrt(2*SNR_dB_Lin));

    % Compute reliability factor: ln((1-Pe)/Pe)
    reliability = log((1-Pe)./Pe);

    % Compute LLRs
    lambda_BSC = reliability .* Y;
end
 