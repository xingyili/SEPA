function res = SEPA(geneExpressionData, adjMatrix)
% SEPA:
[geneNum, geneSample] = size(geneExpressionData);
res = zeros(1, geneSample);
adjMatrix(speye(geneNum,geneNum)==1) = 1;
[row, col] = find(adjMatrix);

% Parameters of solving alpha and beta; 
stopError = 10^-1;
maxIter = 10^4;

for i = 1:geneSample
    singleSample = geneExpressionData(:, i)/sum(geneExpressionData(:, i));
    varB = sparse(row, col, singleSample(row), geneNum, geneNum);
    
    % Iteration parameters
    flagIter = 1;
    curIter = 0;
    alphaPre = singleSample;
    betaPre = ones(geneNum,1);
    varA = adjMatrix;
    while(flagIter)
        alpha = 1./(varA*betaPre);
        beta = singleSample./(varB'*alpha);
        curError = norm([alpha-alphaPre; beta-betaPre], Inf);
        alphaPre = alpha;
        betaPre = beta;
        curIter = curIter + 1;
        flagIter = (curError>stopError && curIter<maxIter);
    end
    % Result of Eq(7)
    res(i) = - dot(singleSample, log(alpha.*beta));
end
eigval=eig(adjMatrix);
lamda=max(eigval);

maxSEPA = log(lamda); 
res = res/maxSEPA;
end

