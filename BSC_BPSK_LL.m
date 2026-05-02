 function [ll_p0,ll_p1] = BSC_BPSK_LL(Y, Sigma, SNR_dB)
 
    % Convert Eb/N0 from dB to linear
    SNR_dB_Lin = 10.^(SNR_dB/10);
    
    % Compute error probability
    Pe = qfunc(sqrt(2*SNR_dB_Lin));

    reliability1 = log(1-Pe);
    reliability2 = log(Pe);

    [ll_p0,ll_p1] = AWGN_BPSK_LL(Y,Sigma);
    
    
    % ---- FORCE COLUMN VECTORS ----
    ll_p0 = ll_p0(:);
    ll_p1 = ll_p1(:);
    
    for index = 1:length(ll_p1)
        if ll_p0(index) > ll_p1(index)
            ll_p0(index) = reliability1;
            ll_p1(index) = reliability2;
        else 
            ll_p1(index) = reliability1;
            ll_p0(index) = reliability2;
        end
    end 
 end