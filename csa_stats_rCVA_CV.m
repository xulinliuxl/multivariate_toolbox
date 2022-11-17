function [CVA,CCApart,CCAperm] = csa_stats_rCVA_CV(CVA)

% Unpack CVA
% varnames = fieldnames(CVA);
% for iVar = 1:numel(varnames)
%     varname = varnames{iVar};
%     eval(sprintf('%s = CVA.%s;',varname,varname));
% end

% ---------------------------
% Unpack CV settings
cfg             = CVA.mode.cv;
CVtype          = cfg.CVtype;
numPart         = cfg.numPart;
numFolds        = cfg.numFolds;
permutePW       = cfg.permutePW;
numPermPW       = cfg.numPermPW;
permuteFW       = cfg.permuteFW';
numPermFW       = cfg.numPermFW;
partSeed        = cfg.partSeed;
doCovFW         = cfg.doCovFW;
labels          = cfg.labels;
numBootDW       = cfg.numBootDW;
doSplitHalfDW   = cfg.doSplitHalfDW;
doSplitHalfFW   = cfg.doSplitHalfFW;
doSaveNullsFW   = cfg.doSaveNullsFW;
doSaveNullsDW   = cfg.doSaveNullsDW;
% ----------------------------------
lambdaX         = CVA.lambdaX;
lambdaY         = CVA.lambdaY;
X               = CVA.X;
Yorig           = CVA.Y;
numVarX         = CVA.numVarX;
numVarY         = CVA.numVarX;
doCovs          = CVA.doCovs;
Ns              = CVA.Ns;
numComp         = CVA.numComp;
usePresetOrder  = CVA.usePresetRandOrder;
presetRandOrder = CVA.randOrder; 
% ----------------------------------
numComp   = min([numComp,size(X,2)...% set nComp, needs number of CVs less than number of variables of the smaller set
                size(X,1)*(1-(1/numFolds))...
                size(Yorig,2)]);

numIter   = numPart;
if strcmp(cfg.CVtype,'LeaveOut')% Reset numFolds and numPart when LeaveOut is called
    numFolds    = 1;
    numPart     = 1;
end


% ------------------------
% Flag Paraller processing
if ~CVA.runInSerial
    parforArg   = Inf;
else
    parforArg   = 0;
end
                    
% --------------------------------------------
% Initialise variable names to collect outputs
% --------------------------------------------
varNames    = {};
varNames    = {'Rfw','Pfw','Cfw'};%,'Rdw','Pdw','Cdw','RAraw','NRAraw','RAte','NRAte','RAtr','NRAtr','RAoos','NRAoos'};     

numOffdiag = [];
if doSplitHalfFW % Additional settings for Split-half Resampling
    varNames = [varNames 'RhsfwXL' 'RhsfwYL' 'RhsfwXW' 'RhsfwYW' 'RhsfwS'];
%     CVtype   = 'kFold'; 

    % -------------------------------------------------------------
    % Build binary matrix to identify folds withing a component
    % as an efficient way to extract data avoiding the use of loops
    % ------------------------------------------------------------
    sizeCM = numFolds*numComp;
    idxOffdiag = zeros(sizeCM);
    for ii = 1:numFolds % create the matrix
        idxTemp = triu(ones(sizeCM),numComp*ii)-triu(ones(sizeCM),(numComp*ii)+1);% Select elements on the off nth offdiagonal
        idxOffdiag = (idxOffdiag + idxTemp);
    end
    numOffdiag = double(idxOffdiag);
    idxOffdiag = logical(idxOffdiag);
    numtemp    = find(idxOffdiag);
    numOffdiag(numtemp) = numtemp;
end
  
 % Fold-wise permutations are implemented by additional (random-order)
 % paritionings of the entire sample
if permutePW
    numIter = numIter + numPermPW;
    varNames = [varNames 'Rdw','Pdw','Cdw'];
end

if ~permuteFW
    numPermFW = 1;
    parforArg = 0;
end

% Labelling samples to ensure partitions with equal numbersdrawn from each group/decile
if isempty(labels) 
    labels    = ones(Ns,1);       
end 


% spmd  % Give each worker the same stream, one that supports substreatms
%     rng(partSeed,'combRecursive');
% end   

%%
% randOrderOut = palm_quickperms(Ns,[],numIter);
%  
% ----------------------
% Spin-test permutations for spatially correspondence of brain maps
coord_l = [];
coord_r = [];
randOrderOut = [];


if usePresetOrder
   randOrderOut = presetRandOrder(:,1:numIter);
end 


% PARFOR-------------------------------------------------------------------------
%parfor iIter = 1: numIter
for iIter = 1: numIter
    iIter
    % Set the substream index by the loop index. This ensures that
    % each interation uses its particular set of random numbers,
    % regardless of which workers run that iteration or what
    % sequence interations run in.
%         stream = RandStream.getGlobalSream();
%         stream.Substream = iPart+222;

    if iIter > numPart 
        if usePresetOrder
            randOrder = presetRandOrder(:,iIter);
        else
            randOrder = randperm(Ns);
        end
        labelIdx  = labels(randOrder);
        Y         = Yorig(randOrder,:);
    else
        Y         = Yorig;
        labelIdx  = labels;
    end
    cv_label = cvpartition(labelIdx,CVtype,numFolds);  %% partitions with equal numbersdrawn from each decile
    cv_label = repartition(cv_label);
   
    
    
    randOrderFW = [];
    if usePresetOrder
%         randIdx     = randperm(size(presetRandOrder,2),numPermFW);
        randOrderFW = presetRandOrder(:,1:numPermFW);
        
%         % Recreate cv_partition for spin rotations, where each rotation is
%         % split in 5 folds without shuffling of nodes. N.B. the shuffling 
%         % comes from the spin rotations
        trainSize = cv_label.TrainSize;
        testSize  = cv_label.TestSize;
        nfolds    = cv_label.NumTestSets;
        dummyvec  = zeros(Ns,1);
        for ifold = 1:nfolds
            istart = 1+sum(testSize(1:ifold-1));
            iend   = sum(testSize(1:ifold));
            dummyvec(istart:iend) = ifold;
        end
        dummies = logical(dummyvar(dummyvec));
        dummymatrix = zeros(size(dummies));
        dummyvec = zeros(Ns,1);
        randorder = presetRandOrder(:,end-iIter);
        cv_label = [];
        for ifold = 1:nfolds
            idxfold = randorder(dummies(:,ifold));
            dummyvec(idxfold) = ifold;
            dummymatrix(idxfold,ifold) = 1;
%             cv_label.test{ifold}  = dummymatrix(:,ifold) ;
%             cv_label.training{ifold} = ~dummymatrix(:,ifold);
        end
        cv_label = cvpartition(dummyvec,CVtype,numFolds);  %% partitions with equal numbersdrawn from each decile
   
%         cv_label.NumTestSets = nfolds;
%         cv_label.N     = Ns;
    else
        randOrderFW = palm_quickperms(Ns,[],numPermFW);
    end
    % ---------------------------------------------------------------------
    % Save Nulls for each partition in separate directories
    dirOut = [];
    if doSaveNullsFW
        dirOut = fullfile(CVA.dirOutNulls,sprintf('partition%04d',iIter));
        mkdir(dirOut);
    end
    
    % --------------------------------------------------------------------
    % Inner loop for Cross-Validations (with Fold-wise Permutations)
    cfg                 = [];
    cfg.randOrder       = randOrderFW;
    cfg.numFolds        = cv_label.NumTestSets;% numFolds;
    cfg.numComp         = numComp;
    cfg.lambdaX         = lambdaX;
    cfg.lambdaY         = lambdaY;
    cfg.Ns              = Ns;
    cfg.doCovs          = 0;
    cfg.doCovFW         = doCovFW;
    cfg.cv_label        = cv_label;
    cfg.X               = X;
    cfg.Y               = Y;
    cfg.permutePW       = permutePW;
    cfg.numPermFW       = numPermFW;
    cfg.doSplitHalfFW   = doSplitHalfFW;
    cfg.varNames        = varNames;
    cfg.doSaveNullsFW   = doSaveNullsFW;
    cfg.dirOut          = dirOut;
  [CCAtemp,resultsNull] = csa_stats_rCVA_CV_innerloop(cfg);

    
    % ---------------------------------------------------------------------
    % Split-half Stability Data-wise, similar to Kovachevic et al 2013
    % ---------------------------------------------------------------------
    if doSplitHalfDW
        
        % Estimate ucorr and vcorr based on half splits
        XwSSr = [];
        XwSSp = [];
        YwSSr = [];
        YwSSp = [];

        if iIter > numPart
            randOrder = randperm(Ns);
            labelIdx  = labels(randOrder);
            Y         = Yorig(randOrder,:);
            numBoot = 30;
        else
            Y         = Yorig;
            labelIdx  = labels;
            numBoot   = numBootDW;
        end

        split_label = cvpartition(labelIdx,'kFold',2);  %% partitions with equal numbersdrawn from each decile
%         numFolds = split_label.NumTestSets;
%         usePresetRandOrder
%     randIdx      = randperm(size(presetRandOrder,2),numIter);
    
%             randOrderOut = presetRandOrder(:,1:numIter);

        XW = cell(numBoot,1);
        YW = cell(numBoot,1);
        for iboot = 1:numBoot
            if usePresetOrder
                numRandSplit = floor(Ns/2) + round(mod(Ns,2)*rand); %  Choose random split for odd Ns 
                idx_1 = randOrderOut(1:numRandSplit,iboot);
                idx_2 = randOrderOut(numRandSplit+1:end,iboot);
            else
                split_label_in= repartition(split_label);
                idx_1   = split_label_in.training(1); % Index for subjects in the first half split
                idx_2   = split_label_in.test(1);     % Index for subjects in the second half split 
            end
            % Select subjects in each split
            X1  = X(idx_1,:);
            X2  = X(idx_2,:);
            Y1  = Y(idx_1,:);
            Y2  = Y(idx_2,:);
            XS1 = CCAtemp.valXs(idx_1,:);
            XS2 = CCAtemp.valXs(idx_2,:);
            YS1 = CCAtemp.valYs(idx_1,:);
            YS2 = CCAtemp.valYs(idx_2,:);

            % Compute Weights (for more reference see csa_stats_rCVA_permSplithalf_inner) 
            Xw1 = X1' * YS1;% Can replace Y1 * Yw with Ys for samples from this split
            Xw2 = X2' * YS2;
            Yw1 = Y1' * XS1;
            Yw2 = Y2' * XS2;

            
            % Correlations b/n projected left and right split-half patterns
            % These are taken as measres of the correspondence b/n X data
            % and Y weights on one hand, and Y data and X weights, on the
            % other hand
            [r,p] = corr(Xw1,Xw2);
            XwSSr(:,iboot) = diag(r); % Component by interation array
            XwSSp(:,iboot) = diag(p);
            [r,p] = corr(Yw1,Yw2);
            YwSSr(:,iboot) = diag(r); % Component by interation array
            YwSSp(:,iboot) = diag(p);
            
%             XW{iPermIn} = (zscore(Xw1)+zscore(Xw2))./2;
%             YW{iPermIn} = (zscore(Yw1)+zscore(Yw2))./2;
            % --------------------------------
            % Could consider Loadings
%             % Xs = X * Xw; XL = X' * Xs;
%             XL1 = X1' * X1 * Xw1;
%             XL2 = X2' * X2 * Xw2;
%             YL1 = Y1' * Y1 * Yw1;
%             YL2 = Y2' * Y2 * Yw2;

%             [r,p] = corr(XL1,XL2);
%             XLSSr(:,iPermIn) = diag(r); % Component by interation array
%             XLSSp(:,iPermIn) = diag(p);
%             [r,p] = corr(YL1,YL2);
%             YLSSr(:,iPermIn) = diag(r); % Component by interation array
%             YLSSp(:,iPermIn) = diag(p);
            
            
            
        end % iPermIn 

        CCAtemp.RbootXw = median(abs(XwSSr),2);
        CCAtemp.pPermRbootXw = median(XwSSp,2);
        CCAtemp.RbootYw = median(abs(YwSSr),2);
        CCAtemp.pPermRbootYw = median(YwSSp,2);

    end
    
    
%     if doSavePermParams
%         CCAtemp.permXL = reshape([resultsNull.valXL],[numVarX,numComp,numPermFW-1]);
%         CCAtemp.permXLfw = reshape([resultsNull.valXLfw],[numVarX,numComp,numPermFW-1]);
%     end
%     varNames = CCAtemp.varNames;
    if doSaveNullsFW
        CCAtemp.dirOutNulls = dirOut; 
    end
   
    if permuteFW
        resultsNull(1)= []; % Remove the first entry, as this is from the original order data
        
        for ivar = 1:numel(varNames)
            varName = varNames{ivar};
            varNamePerm = ['pPerm' varName];
            varNameNull = ['null' varName];
            
            dataNull = [resultsNull.(varNameNull)];
%             eval(sprintf('null%s(1,:) = [];',varNames{ivar}));% First value is with true labelling, so remove
            if strcmp(varNames{ivar}(1),'P') % For Pvalues search left side of the distribution
%                 eval(sprintf('CCApart(iIter).pPerm%s = (sum(abs(null%s) < abs(repmat(CCApart(iIter).%s'',numPermFW-1,1)))/numPermFW)'';',varNames{ivar},varNames{ivar},varNames{ivar}));
                
                 CCAtemp.(varNamePerm) = (sum(abs(dataNull) < abs(repmat(CCAtemp.(varName),1,numPermFW-1)),2)/numPermFW);
            else
%                 eval(sprintf('CCApart(iIter).pPerm%s = (sum(abs(null%s) > abs(repmat(CCApart(iIter).%s'',numPermFW-1,1)))/numPermFW)'';',varNames{ivar},varNames{ivar},varNames{ivar}));

                CCAtemp.(varNamePerm) = (sum(abs(dataNull) > abs(repmat(CCAtemp.(varName),1,numPermFW-1)),2)/numPermFW);
            end
        end
    end
    
    iCCA(iIter)      = CCAtemp;

    if doSaveNullsDW
        if iIter <= numPart % Use differnt names for partitions and permutations
            fout = fullfile(CVA.dirOutNulls,sprintf('results_part%06d.mat',iIter));
        else 
            fout = fullfile(CVA.dirOutNulls,sprintf('results_null%06d.mat',iIter));
        end
        kat_parfor_save(fout,{'CCAout'},CCAtemp);
    end
    
    CCAtemp = [];
    pNullDW = [];
    rNullDW = [];
    pNullFW = [];
    rNullFW = [];
    resultsNull=[];

end % End Iter loop
% -------------------------------------------------------------------------
CVA.varNames = varNames;
CCAperm = [];
CCApart = [];

% -------------------------------------------------------------------------
% Assemble results across permutations/pertitions, unless they are saved
% on request (e.g. for high number of perms/partitions CCAtemp could be 
% large in size causing memory issues).
if ~doSaveNullsDW

    % ------------------------------------------------------------
    % For partition-wise permutations save permutations to CCAperm
    try
        CCAperm = iCCA(numPart+1:end);
        CCApart = iCCA(1:numPart);
    end
    iCCA = [];
    % ---------------------------------------------------------------------------------------
    % There can be sign-flips across the Pertitionins/Iterations, correct for that (taking Iteration 1 as reference)
    % ---------------------------------------------------------------------------------------
    % warning('---- There can be sign-flips across the folds, correct for that (taking fold 1 as reference)\n');
    if numPart>1
        for nn=2:numPart
            XLsign = diag(corr(CCApart(1).valXL,CCApart(nn).valXL));
            YLsign = diag(corr(CCApart(1).valYL,CCApart(nn).valYL));
            XSsign = diag(corr(CCApart(1).valXs,CCApart(nn).valXs));
            YSsign = diag(corr(CCApart(1).valYs,CCApart(nn).valYs));
            signFlip = sign(mean([XLsign YLsign XSsign YSsign],2));

            CCApart(nn).valXL = CCApart(nn).valXL.*signFlip';
            CCApart(nn).valYL = CCApart(nn).valYL.*signFlip';
            CCApart(nn).valXs = CCApart(nn).valXs.*signFlip';
            CCApart(nn).valYs = CCApart(nn).valYs.*signFlip';

            if doSplitHalfFW
                CCApart(nn).valXLhsfw = CCApart(nn).valXLhsfw.*signFlip';
                CCApart(nn).valYLhsfw = CCApart(nn).valYLhsfw.*signFlip';
                CCApart(nn).valXWhsfw = CCApart(nn).valXWhsfw.*signFlip';
                CCApart(nn).valYWhsfw = CCApart(nn).valYWhsfw.*signFlip';
            end
        end
    end

    %% Assemble output CVA

    % CVA.opt    = CCApart(1).opt;
    % -------------------------------------------------------------------------
    % Variable names to average across partitions (Loadings, Scores, Weights),
    % not statistic measures 
    % -------------------------------------------------------------------------
    % paramNames = {'valXs','valXL','valXLp','valXLhsfw','valXWhsfw','valYs','valYL','valYLp','valYLhsfw','valYWhsfw'};
    % % if doSavePermParams
    % %     paramNames = [paramNames 'permXL' 'permXLfw'];
    % % end
    % for ivar = 1:numel(paramNames)
    %     parname = paramNames{ivar};
    %     dims = [size(CCApart(1).(parname)) numPart];
    %     dim2average = numel(dims)+1;
    %     eval(sprintf('CVA.%s = median(reshape([CCApart.%s],dims),dim2average);',...
    %                   parname,parname));
    % end

    % -------------------------------------------------------------------------
    % Averaging across partitions. 
    % use traditional naming convention for out-of-sample params, i.e. remove val* prefix 
    CVA.XL  = median(reshape([CCApart.valXL],[size(CCApart(1).valXL,1),size(CCApart(1).valXL,2),numPart]),3);
    CVA.XLp = median(reshape([CCApart.valXLp],[size(CCApart(1).valXLp,1),size(CCApart(1).valXLp,2),numPart]),3);
    CVA.YL  = median(reshape([CCApart.valYL],[size(CCApart(1).valYL,1),size(CCApart(1).valYL,2),numPart]),3);
    CVA.YLp = median(reshape([CCApart.valYLp],[size(CCApart(1).valYLp,1),size(CCApart(1).valYLp,2),numPart]),3);
    CVA.XS  = median(reshape([CCApart.valXs],[size(CCApart(1).valXs,1),size(CCApart(1).valXs,2),numPart]),3);
    CVA.YS  = median(reshape([CCApart.valYs],[size(CCApart(1).valYs,1),size(CCApart(1).valYs,2),numPart]),3);

    if doSplitHalfFW
        CVA.XLhsfw  = median(reshape([CCApart.valXLhsfw],[size(CCApart(1).valXLhsfw,1),size(CCApart(1).valXLhsfw,2),numPart]),3);
        CVA.YLhsfw  = median(reshape([CCApart.valYLhsfw],[size(CCApart(1).valYLhsfw,1),size(CCApart(1).valYLhsfw,2),numPart]),3);
        CVA.XWhsfw  = median(reshape([CCApart.valXWhsfw],[size(CCApart(1).valXWhsfw,1),size(CCApart(1).valXWhsfw,2),numPart]),3);
        CVA.YWhsfw  = median(reshape([CCApart.valYWhsfw],[size(CCApart(1).valYWhsfw,1),size(CCApart(1).valYWhsfw,2),numPart]),3);
    end

    % CVA.XL  = CVA.valXL;  % median(reshape([CCApart.valXL],[size(CCApart(1).valXL,1),size(CCApart(1).valXL,2),numPart]),3);
    % CVA.XLp = CVA.valXLp; %median(reshape([CCApart.valXLp],[size(CCApart(1).valXLp,1),size(CCApart(1).valXLp,2),numPart]),3);
    % CVA.YL  = CVA.valYL;  %median(reshape([CCApart.valYL],[size(CCApart(1).valYL,1),size(CCApart(1).valYL,2),numPart]),3);
    % CVA.YLp = CVA.valYLp; %median(reshape([CCApart.valYLp],[size(CCApart(1).valYLp,1),size(CCApart(1).valYLp,2),numPart]),3);
    % CVA.XS  = CVA.valXs;  %median(reshape([CCApart.valXs],[size(CCApart(1).valXs,1),size(CCApart(1).valXs,2),numPart]),3);
    % CVA.YS  = CVA.valYs;  %median(reshape([CCApart.valYs],[size(CCApart(1).valYs,1),size(CCApart(1).valYs,2),numPart]),3);

    % ---------------------------------------------------------------------------------
    % Average statistical measures (CVs Rval and pvals) across all paritions/iterations
    % ---------------------------------------------------------------------------------

    % Additional variables are defined here (not at the start), as they are
    % only needed for the loops below
    if doSplitHalfDW
        varNames = [varNames 'RbootXw' 'RbootYw'];% 'XwCorrR','XwCorrP','YwCorrR','YwCorrP' ];
    end
    for ivar = 1:numel(varNames)
        varname = varNames{ivar};
        CVA.(varname) = median([CCApart.(varname)],2);
        eval(sprintf('CVA.pPermFW%s = median([CCApart.pPerm%s],2)'';',varname,varname));

        if permutePW
            eval(sprintf('CVA.pPermPW%s = sum([CCAperm.%s]'' > abs(repmat(CVA.%s'',numPermPW,1)))/numPermPW;',varname,varname,varname));
        end
    end
end

% 
% if permutePW
%     
%     % Permutation-based null distribution across CVs and Folds
%     for ivar = 1:numel(varNames)
%         componentNull    = cell2mat(arrayfun(@(x) abs(x.(varNames{ivar})),CCAperm,'UniformOutput',0));
%         globalNull  = max(componentNull);
%         % Calculate the percentiles for each statistic based on permutation-based null dist
%         p_perm = zeros(numComp,1);
%         for ii = 1:numComp
%             if strcmp(varNames{ivar}(1),'P') % For Pvalues search left side of the distribution
%                 eval(sprintf('pPermGlobalDist%s(ii,1)    = (sum(globalNull <= abs(CVA.%s(ii)))/numPermPW)'';',varNames{ivar},varNames{ivar})); 
%                 eval(sprintf('pPermComponentDist%s(ii,1) = (sum(componentNull(ii,:) <= abs(CVA.%s(ii)))/numPermPW)'';',varNames{ivar},varNames{ivar})); 
%             else
%                 eval(sprintf('pPermGlobalDist%s(ii,1)    = (sum(globalNull >= abs(CVA.%s(ii)))/numPermPW)'';',varNames{ivar},varNames{ivar})); 
%                 eval(sprintf('pPermComponentDist%s(ii,1) = (sum(componentNull(ii,:) >= abs(CVA.%s(ii)))/numPermPW)'';',varNames{ivar},varNames{ivar})); 
%             end
% 
%         end
%         eval(sprintf('CVA.pPermGlobalDist%s = pPermGlobalDist%s'';',varNames{ivar},varNames{ivar}));
%         eval(sprintf('CVA.pPermComponentDist%s = pPermComponentDist%s'';',varNames{ivar},varNames{ivar}));
%     end
% end


%     if doSplitHalfDW
%         CCAtemp.pPermXwCorrR = sum(abs(nullXwCorrRmedian)' >= abs(repmat(XwCorrRmedian',numPermPW,1)))/numPermPW;
%         CCAtemp.pPermXwCorrP = sum(abs(nullXwCorrPmedian)' <= abs(repmat(XwCorrPmedian',numPermPW,1)))/numPermPW;
%         CCAtemp.pPermYwCorrR = sum(abs(nullYwCorrRmedian)' >= abs(repmat(YwCorrRmedian',numPermPW,1)))/numPermPW;
%         CCAtemp.pPermYwCorrP = sum(abs(nullYwCorrPmedian)' <= abs(repmat(YwCorrPmedian',numPermPW,1)))/numPermPW; 
%     end