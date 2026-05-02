function Score = SCL_Frozen_Bits_Cost(BER_Vec, Flavor, L_Stages, Is_Frozen_Vec)
% Prototype preserved.
% Returns Score in BER units (0..0.5), normalized by windowLen.

windowLen = L_Stages + 1;

% ---------- compute badness same as optimizer ----------
p = BER_Vec(:)';
switch upper(Flavor)
    case 'MI'
        H = zeros(size(p));
        m = (p>0) & (p<1);
        H(m) = - p(m).*log2(p(m)) - (1-p(m)).*log2(1-p(m));
        H(p==0) = 0;
        MI = 1 - H;
        Badness = -MI;
    case 'LLR'
        LLR = zeros(size(p));
        m = (p>0) & (p<1);
        LLR(m) = log((1-p(m))./p(m));
        LLR(p==0) = 1e6;
        LLR(p==1) = -1e6;
        Badness = -LLR;
    otherwise
        error('Unsupported Flavor "%s"', Flavor);
end

% ---------- kept subsequence ----------
seq = Badness(Is_Frozen_Vec==0);
K = numel(seq);

if K == 0
    raw = 0;
elseif K < windowLen
    raw = sum(seq);
else
    c = cumsum([0 seq]);
    wins = c(windowLen+1:end) - c(1:end-windowLen);
    raw = max(wins);
end

avgMetric = raw / windowLen;
Score = invert_avg_to_BER_cost(avgMetric, Flavor);

end

% helper invert for cost (same logic as optimizer's invert)
function p = invert_avg_to_BER_cost(avgMetric, Flavor)
switch upper(Flavor)
    case 'MI'
        Havg = avgMetric+1;
        if Havg <= 0, p = 0; return; end
        if Havg >= 1, p = 0.5; return; end
        f = @(x) - x.*log2(max(x,eps)) - (1-x).*log2(max(1-x,eps)) - Havg;
        try
            p = fzero(f,[eps 0.5]);
            if ~isfinite(p) || p<0 || p>0.5
                p = fminbnd(@(x) abs(f(x)), eps, 0.5);
            end
        catch
            p = fminbnd(@(x) abs(f(x)), eps, 0.5);
        end
    case 'LLR'
        avgLLR = -avgMetric;
        p = 1 / (1 + exp(avgLLR));
        p = min(max(p,0),0.5);
    otherwise
        error('Unsupported Flavor "%s"', Flavor);
end
end
