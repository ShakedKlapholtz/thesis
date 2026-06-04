function [Estimated_U,Estimated_L] = gArikan_BPSK_SCL_Decoder(L,GA_CRC_Length,Y,U,Sigma,SNR_dB,Is_Frozen_Bit_Index_Vec,mode)

    N = length(Y);
    
    if nargin < 8
        mode = 'BSC';
    end

    R = sum(Is_Frozen_Bit_Index_Vec==0)/length(Is_Frozen_Bit_Index_Vec);
    EbN0 = SNR_dB - 10*log10(R);
    EbN0_lin = 10^(EbN0/10);

    switch mode
    
        case 'AWGN'
            % Soft LLR
            [LL0,LL1] = AWGN_BPSK_LL(Y, Sigma);
    
        case 'BSC'
            % % Hard decision
            % Y_hard = (Y < 0);
            % 
            %   % Crossover probability
            % p = qfunc(sqrt(2*R*EbN0_lin));
    
            % Your function
            [LL0, LL1] = BSC_BPSK_LL(Y, Sigma, SNR_dB); % Assuming this returns N x 1
            
            % Force to Row Vector (1 x N)
            LL0 = LL0(:)'; 
            LL1 = LL1(:)'; 

    end
    
    
    [Estimated_U_List,Estimated_X_List,Estimated_L0_List,Estimated_L1_List,~] = gArikan_SCL_Decoder(LL0,LL1,L,U,N,Is_Frozen_Bit_Index_Vec);
    
    is_correct = ismember(U,Estimated_U_List,'rows');
    
    if(and(GA_CRC_Length>=0,is_correct))
        
        Index = find(sum(bsxfun(@eq,Estimated_U_List,U),2)==N);
        
        Estimated_U = U;

        Estimated_L = Estimated_L0_List(Index,:)-Estimated_L1_List(Index,:);
        
    else

        p = (~Estimated_X_List)*(LL0.')+(Estimated_X_List)*(LL1.');

        [~,ML_Index] = max(p);
        
        Estimated_U = Estimated_U_List(ML_Index,:);

        Estimated_L = Estimated_L0_List(ML_Index,:)-Estimated_L1_List(ML_Index,:);
        
    end
    
 end