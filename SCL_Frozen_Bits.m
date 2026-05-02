function [SCL_Is_Frozen_Vec] = SCL_Frozen_Bits(BER_Vec, Is_Valid_Statistics, k, L, Flavor, Is_Print)

    BER_Vec = min(BER_Vec,0.5);

    N = length(BER_Vec);

    X_Vec = 1:N;

    L_Stages = log2(L);
    
    if(Is_Print)
        figure; hold on; grid minor; plot(BER_Vec,'.b'); ylim([0,0.5]); ylabel('BER'); xlabel('Index');
    end
    
    SCL_Is_Frozen_Vec = zeros(L_Stages+1,N);
    Best_Score_Vec = zeros(1,L_Stages+1);
    for i_L = 0:1:L_Stages
       
        if(i_L == 0)
            [~,Indexes] = sort(BER_Vec);
            SCL_Is_Frozen_Vec(1,Indexes(k+1:end)) = 1;
            if(SCL_Is_Frozen_Vec(~Is_Valid_Statistics)==1)
                error("Insufficient statistics");
            end
            Best_Score_Vec_Temp = SCL_Frozen_Bits_Cost(BER_Vec, Flavor, 0, SCL_Is_Frozen_Vec(1,:));
            
            [SCL_Is_Frozen_Vec_Temp,Best_Score_Vec(1)] = SCL_Frozen_Bits_Optimization_With_Init_Guess(BER_Vec, Flavor, k, i_L, SCL_Is_Frozen_Vec(1,:));
            
            if(abs(Best_Score_Vec_Temp-max(BER_Vec(~SCL_Is_Frozen_Vec(1,:)),[],'all'))>Best_Score_Vec(1)/100)
                error("Wrong freezing for SC");
            end
            
            if( (~isequal(SCL_Is_Frozen_Vec_Temp,SCL_Is_Frozen_Vec(1,:))) && (Best_Score_Vec(1)~=Best_Score_Vec_Temp) )
                error("Bug: output for L=1 was different than regular freezing for SC");
            end
        else
            [SCL_Is_Frozen_Vec(i_L+1,:), Best_Score_Vec(i_L+1)] = SCL_Frozen_Bits_Optimization_With_Init_Guess(BER_Vec, Flavor, k, i_L, SCL_Is_Frozen_Vec(i_L,:));
            Best_Score_Vec_Temp = SCL_Frozen_Bits_Cost(BER_Vec, Flavor, i_L, SCL_Is_Frozen_Vec(i_L,:));
            if(Best_Score_Vec_Temp<Best_Score_Vec(i_L+1))
                error("Initial guess is better than output");
            end
        end
        
        if(Is_Print)
            plot(X_Vec(~SCL_Is_Frozen_Vec(i_L+1,:)),BER_Vec(~SCL_Is_Frozen_Vec(i_L+1,:)),'o');
        end
        
    end

end

