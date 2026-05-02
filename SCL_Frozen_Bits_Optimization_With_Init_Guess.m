function [Is_Frozen_Vec, Best_Score] = SCL_Frozen_Bits_Optimization_With_Init_Guess(BER_Vec, Flavor, k, L_Stages, Is_Frozen_Vec_Guess)

% maxBatch = L_Stages + 1;  % flip at most windowLen bits simultaneously instead of fixed 10

% maxBatch = 2;
maxBatch = 1;
maxBatch = min(maxBatch,L_Stages + 1);  % flip at most windowLen bits simultaneously instead of fixed 10

% Change_Percentage = 0.01; % Compared with min(k,n-k) (e.g. if n=2048 and k=1024 and Change_Percentage = 0.1, round(0.1*1024)=102 changes from guess)
% Min_Max_Change = 10; % Try to change up to at least Min_Max_Change bits
% TH_Keep = 0.01;  % BER threshold for “always keep” bits
% TH_Keep = 0.00;  % BER threshold for “always keep” bits
% TH_Frozen = 0.1;  % BER threshold for “always freeze” bits
% TH_Frozen = 0.4;  % BER threshold for “always freeze” bits
% TH_Frozen = 0.5;  % BER threshold for “always freeze” bits

TH_Keep = max(BER_Vec(~Is_Frozen_Vec_Guess))/4;  % BER threshold for “always keep” bits
TH_Frozen = max(BER_Vec(~Is_Frozen_Vec_Guess))*4;  % BER threshold for “always freeze” bits

% TH = 4*max(BER_Vec(~Is_Frozen_Vec_Guess));
    
n = numel(BER_Vec);
windowLen = L_Stages + 1;

% ---------- initial guess ----------
if isempty(Is_Frozen_Vec_Guess)
    [~, idxAsc] = sort(BER_Vec,'ascend');
    mask = ones(1,n);
    mask(idxAsc(1:k)) = 0;
else
    mask = double(Is_Frozen_Vec_Guess(:)');
    if numel(mask) ~= n, error('Is_Frozen_Vec_Guess must have length n'); end
    if sum(mask==0) ~= k, error('Is_Frozen_Vec_Guess must contain exactly k zeros'); end
end

% ---------- score of initial guess ----------
Guess_Score = SCL_Frozen_Bits_Cost(BER_Vec, Flavor, L_Stages, mask);

fprintf('\nTrying to improve for L=%d\nGuess score = %.4f\n', 2^(L_Stages), Guess_Score);

% ---------- permanent-kept optimization ----------
% permanentKeep = false(1,n);
% for idx = 1:n
%     if BER_Vec(idx) == 0
%         startIdx = max(1, idx - windowLen + 1);
%         endIdx   = min(n, idx + windowLen - 1);
%         if all(BER_Vec(startIdx:endIdx) == 0)
%             permanentKeep(startIdx:endIdx) = true;
%         end
%     end
% end

permanentKeep = false(1,n);

for i = 1:n
    startIdx = max(1, i-windowLen);
    endIdx   = min(n, i+windowLen);
    
    if all(BER_Vec(startIdx:endIdx) < TH_Keep)
        permanentKeep(startIdx:endIdx) = true;
    end
end

% Only make permanent keep if it was initially keep
permanentKeep = permanentKeep & (Is_Frozen_Vec_Guess == 0);

permanentFrozen = false(1,n);

for i = 1:n
    startIdx = max(1, i-(L_Stages+1));
    endIdx   = min(n, i+(L_Stages+1));
    
    if all(BER_Vec(startIdx:endIdx) > TH_Frozen)
        permanentFrozen(startIdx:endIdx) = true;
    end
end

% Only make permanent frozen if it was initially frozen
permanentFrozen = permanentFrozen & (Is_Frozen_Vec_Guess == 1);

% 0 = kept, 1 = frozen
mask(permanentFrozen) = 1;  
mask(permanentKeep) = 0;

% ---------- batch-flip heuristic near guess ----------
% maxTotalChange = max(1, max(Min_Max_Change,round(Change_Percentage*min(k,n-k)))); % consider up to 10% of bits

improved = true;
while improved
    improved = false;
    
    % candidate bits for swap (excluding permanent-kept)
    keptIdx   = find(mask==0 & ~permanentKeep & ~permanentFrozen);
    frozenIdx = find(mask==1 & ~permanentKeep & ~permanentFrozen);
%     keptIdx   = find(mask==0 & ~permanentKeep & ~permanentFrozen & BER_Vec<TH);
%     frozenIdx = find(mask==1 & ~permanentKeep & ~permanentFrozen & BER_Vec>TH);
    
    % limit total change candidates to 10% of k
    keptIdx   = keptIdx(1:numel(keptIdx));
    frozenIdx = frozenIdx(1:numel(frozenIdx));
    
    % generate batch swaps: up to maxBatch bits at a time
    for batchSize = 1:min(maxBatch, min(numel(keptIdx), numel(frozenIdx)))
               
        keptComb   = nchoosek(keptIdx, batchSize);
        frozenComb = nchoosek(frozenIdx, batchSize);
        
        fprintf('Trying to improve with %d simultaneous flips (will continue up to %d simultaneous flips) - %d options\n', batchSize, min(maxBatch, min(numel(keptIdx), numel(frozenIdx))), size(keptComb,1)*size(frozenComb,1));
        
        for i = 1:size(keptComb,1)
            for j = 1:size(frozenComb,1)
                maskTmp = mask;
                maskTmp(keptComb(i,:))   = 1;
                maskTmp(frozenComb(j,:)) = 0;
                
                if sum(maskTmp==0) ~= k
                    continue;
                end
                
                scoreTmp = SCL_Frozen_Bits_Cost(BER_Vec, Flavor, L_Stages, maskTmp);
                if scoreTmp < Guess_Score - 1e-12
                    mask = maskTmp;
                    Guess_Score = scoreTmp;
                    improved = true;
                    fprintf('Improved solution with score = %.4f\n', scoreTmp);
                    break; % accept first improvement
                end
            end
            if improved
                break;
            end
        end
        if improved
            break;
        end
    end
end

if(~improved)
    fprintf('Did not improve the initial guess\n');
end

% ---------- finalize ----------
Is_Frozen_Vec = mask;
Best_Score = SCL_Frozen_Bits_Cost(BER_Vec, Flavor, L_Stages, Is_Frozen_Vec);

% ensure not worse than initial guess
if Best_Score > Guess_Score + 1e-12
    error('Function failed: Best_Score > Guess_Score');
end

end
