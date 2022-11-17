

n = 60; % number of observations
[X,Y] = generate_data(n); % observations from a random normal multivariate distribution

addpath(genpath('/rds/user/xl454/hpc-work/Cam-CAN/FC_SC/csa_tbx'))

%%
% ----------------------------------------
% Construct CCA structure for CCA analysis
% ----------------------------------------
CCA        = [];
CCA.X      = X;
CCA.Y      = Y;

CCA.mode.cv.do              = 0; % Cross-validation settings
CCA.mode.cv.numFolds        = 5;
CCA.mode.cv.numPart         = 1; % How many times to repeat CV with different patitions
CCA.mode.cv.permutePW       = 1; % Partition-wise permutations
CCA.mode.cv.numPermPW       = 1000;
CCA.mode.cv.permuteFW       = 0;
CCA.mode.cv.numPermFW       = 100;
CCA.mode.cv.doSplitHalfDW   = 0;
CCA.mode.cv.doSplitHalfFW   = 0;
CCA.mode.cv.numBootDW       = 100;
CCA.mode.cv.doSaveNullsFW   = 0;
CCA.mode.cv.doSaveNullsDW   = 0;

CCA.mode.permClassic.do     = 0;        % Classical Permutations
CCA.mode.permClassic.numPerm= 10000;

CCA.mode.permBootstr.do       = 0;      % Split-half permutations
CCA.mode.permBootstr.numPerm  = 1000;

CCA.mode.standard.do    = 1;

CCA.numComp             = min([5,...
                       size(CCA.X,2)...
                       size(CCA.Y,2)]);
CCA.doSaveNulls         = 0;
CCA.usePresetRandOrder  = 0;
CCA.nameAnalysis        = 'prelim';
% CCA.dirOut              = S.paths.results;
CCA.lambdaX = 1;% 'auto' or [0-1], where 0 is CCA, 1 is PLS, 0-1 is regularized CCA
CCA.lambdaY = 1;%
tic, [cca] = csa_stats_rCVA_wrapper(CCA); toc


%%
% Show correlation between first canonical variates
figure;scatter(cca.XS(:,1),cca.YS(:,1))
figure;scatter(cca.XS(:,2),cca.YS(:,2))


nnodes =246;
cm_mask = triu(ones(nnodes),1);
cm_idx = find(idx);

numcomp =3;
for icomp =1:numcomp
    cmX = zeros(nnodes);
    cmX(cm_idx) = cca.XS(:,icomp);
    figure;imagesc(cmX)

    cm = zeros(nnodes);
    cm(cm_idx) = cca.YS(:,icomp);
    figure;imagesc(cm);
end


figure;scatter(cca.XL(:,1),cca.YL(:,1))