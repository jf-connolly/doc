%-- Confidence interval
alpha        = 0.10;
confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
pValue       = tinv(confInterval, nReplications-1) / sqrt(nReplications);

tested1 = 4;
tested2 = 7;

mu1 = graphErr(tested1 *3-1,:);
mu2 = graphErr(tested2 *3-1,:);

std1 = graphErr(tested1 *3,:)/pValue;
std2 = graphErr(tested2 *3,:)/pValue;

pValue = tinv(1-alpha/2, nReplications*2-2);
tested = (mu1-mu2)./(sqrt(std1.^2/50 + std2.^2/50));


nBlock = size(graphErr,2);
result = zeros(1,nBlocks);
for t = 1:nBlock
    if abs(tested(t)) > pValue,  result(t) = 1;
    else                         result(t) = 0;
    end
end

result

