function [CCAtemp,resultsNull] = csa_stats_rCVA_CV_innerloop(cfg)
% Permute subjects/samples within training sample (fold-wise)

% ----------------------------
% Unpack cfg structure
randOrder       = cfg.randOrder;
numPermFW       = cfg.numPermFW;
numFolds        = cfg.numFolds;
permutePW       = cfg.permutePW;
numComp         = cfg.numComp;
cv_label        = cfg.cv_label;
varNames        = cfg.varNames;
doSplitHalfFW   = cfg.doSplitHalfFW;
doSaveNullsFW   = cfg.doSaveNullsFW;

resultsNull(numPermFW).nullRfw = [];


if numPermFW == 1
    parforArg = Inf;
else 
    parforArg = 0;
end

for iPerm = 1:numPermFW
% parfor (iPerm = 1:numPermFW,parforArg)%
%     Y       = YorigPerm(randOrderFW(:,iPerm),:);   % for iPerm =1 original order, otherwise reodered randomly
    Y       = cfg.Y(randOrder(:,iPerm),:);   % for iPerm =1 original order, otherwise reodered randomly
    % Initate some variables
    Xs      = cell(numFolds,1);
    Ys      = cell(numFolds,1);
    Xw      = cell(numFolds,1);
    Yw      = cell(numFolds,1);
    XL      = cell(numFolds,1);
    YL      = cell(numFolds,1);
    XSte    = cell(numFolds,1);
    YSte    = cell(numFolds,1);
    XLte    = cell(numFolds,1);
    YLte    = cell(numFolds,1);
    XWte    = cell(numFolds,1);
    YWte    = cell(numFolds,1);
    r_fw    = nan(numComp,numFolds);
    p_fw    = nan(numComp,numFolds);
    c_fw    = nan(numComp,numFolds);
    RAtr    = nan(numComp,numFolds);
    NRAtr   = nan(numComp,numFolds);
    RAte    = nan(numComp,numFolds);
    NRAte   = nan(numComp,numFolds);
    RAraw   = nan(numComp,numFolds);
    NRAraw  = nan(numComp,numFolds);

%     opt          = struct();
%     opt.lambdaX  = lambdaX;
%     opt.lambdaY  = lambdaY;
%     opt.svd      = 0;
%     opt.numComp  = numComp;
%     opt.doCovs   = 0; % N.B. do not apply regression within csa_proc_CCA, as has been done above (Different for cross-validated approaches)
%     opt.verbose  = 0;
%     opt.Ns       = Ns;
%     opt.doCovFW  = doCovFW;
    
    opt = cfg;

    % Training and testing model with optimal lambda
%         warning('*** Training and testing model with optimal lambda\n\n','verbose');
    for nn = 1:numFolds
        % --------------------------------------------------------------------
        % ------------- XVALIDATION FOR OUT-OF-SAMPLE SCORES -----------------
        % --------------------------------------------------------------------
        try 
            train_idx = find(cv_label.training(nn)==1);
            test_idx  = find(cv_label.test(nn)==1);
        catch 
            train_idx = find(cv_label.training{nn}==1);
            test_idx  = find(cv_label.test{nn}==1);
        end
        Ytrain    = Y(train_idx,:);  
        Ytest     = Y(test_idx,:);    
        Xtrain    = cfg.X(train_idx,:);
        Xtest     = cfg.X(test_idx,:);
        opt.Xte   = Xtest;
        opt.Yte   = Ytest;
        opt.numComp = numComp;
        ccaFold = csa_proc_CCA(opt,Xtrain,Ytrain);

        Xw{nn} = ccaFold.XW;    Yw{nn} = ccaFold.YW;
        Xs{nn} = ccaFold.XS;    Ys{nn} = ccaFold.YS;
        XL{nn} = ccaFold.XL;    YL{nn} = ccaFold.YL;
        XSte{nn} = ccaFold.XSte;YSte{nn} = ccaFold.YSte;
        XLte{nn} = ccaFold.XLte;YLte{nn} = ccaFold.YLte;
        XWte{nn} = ccaFold.XWte;YWte{nn} = ccaFold.YWte;

        % Fold-wise out-of-sample correlation between test samples
        r_fw(:,nn)  = ccaFold.Rte;
        p_fw(:,nn)  = ccaFold.Pte;
        c_fw(:,nn)  = ccaFold.Cte;
    end

%         Rfw = mean(r_fw,2);
%         Pfw = mean(p_fw,2);
%         Cfw = mean(c_fw,2);


    % Create variable names and assign all outputs to 'results' structure
    results     = [];
    results.Rfw = mean(r_fw,2);
    results.Pfw = mean(p_fw,2);
    results.Cfw = mean(c_fw,2);

    % ---------------------------------------------------------------------------------------
    % There can be sign-flips across the folds, correct for that (taking fold 1 as reference)
    % ---------------------------------------------------------------------------------------
    % warning('---- There can be sign-flips across the folds, correct for that (taking fold 1 as reference)\n');
    [nx ny] = size(XL{1});
    datXL=reshape([XL{:}],nx,ny,numFolds);
    for comp=1:size(YL{1},2)
        [coeff, score, latent, tsquared, explained] = pca(squeeze(datXL(:,comp,:)));
        signflip = sign(coeff(:,1));
        
        for nn = 1:numFolds
            XSte{nn}(:,comp) = XSte{nn}(:,comp).*signflip(nn);
            YSte{nn}(:,comp) = YSte{nn}(:,comp).*signflip(nn);
            XLte{nn}(:,comp) = XLte{nn}(:,comp).*signflip(nn);
            YLte{nn}(:,comp) = YLte{nn}(:,comp).*signflip(nn);
            XWte{nn}(:,comp) = XWte{nn}(:,comp).*signflip(nn);
            YWte{nn}(:,comp) = YWte{nn}(:,comp).*signflip(nn);
            XL{nn}(:,comp)   = XL{nn}(:,comp).*signflip(nn);
            YL{nn}(:,comp)   = YL{nn}(:,comp).*signflip(nn);
        end
    end
%     for nn=2:numFolds
%         for comp=1:size(YL{1},2)
%             %%% check if scalar product of loadings is negative
%             if  (YL{1}(:,comp)' * YL{nn}(:,comp) < 0)
%                 XSte{nn}(:,comp)= -XSte{nn}(:,comp);
%                 YSte{nn}(:,comp)= -YSte{nn}(:,comp);
%                 YL{nn}(:,comp)  = -YL{nn}(:,comp);
%                 XL{nn}(:,comp)  = -XL{nn}(:,comp);
%                 YLte{nn}(:,comp)  = -YLte{nn}(:,comp);
%                 XLte{nn}(:,comp)  = -XLte{nn}(:,comp);
%                 YWte{nn}(:,comp)  = -YWte{nn}(:,comp);
%                 XWte{nn}(:,comp)  = -XWte{nn}(:,comp);
% 
%                 Output.VERBOSE('X+Y sign flip in fold=%d, comp=%d\n',nn,comp);
%             end
%         end
%     end

    % valXs and valYs how much each subject is loading on the pattern (linear
    % combination) of variables in each dataset, where datasets are optimised
    % to produce maximal correlation (seen in a correlation b/n valXs and valYs)
    % warning('---- Get valXS and valYS (out of sample) scores and collect data across folds\n');
    valXs= nan( cv_label.N, size(XSte{1},2));
    valYs= nan( cv_label.N, size(YSte{1},2));
    for nn=1:cv_label.NumTestSets
        try
            valXs(cv_label.test(nn)==1, :) = XSte{nn}-mean(XSte{nn}); % Demeaning values within each fold
            valYs(cv_label.test(nn)==1, :) = YSte{nn}-mean(YSte{nn});
        catch
             valXs(cv_label.test{nn}==1, :) = XSte{nn}-mean(XSte{nn}); % Demeaning values within each fold
            valYs(cv_label.test{nn}==1, :) = YSte{nn}-mean(YSte{nn});
        end
    end


    %%% N.B. small hack (TEMPORARY!!!), while taking real values in valXs
    %%% mt_proc_CCA.m
%         if any(isnan(valXs(:))) | any(isinf(valXs(:)))
%             continue
%         end

    %%% -----------------------------------------------------------------------
    %%% This is important: The xscores are designed to be orthogonal, but the OUT-OF-SAMPLE
    %%% xscores can be correlated: hence, we must orthogonalise each CCA component
    %%% w.r.t to the preceding components so that each component only explains a unique
    %%% variance !!!
    %%% Also orthogonalise w.r.t. TMOV!!
    %%% -----------------------------------------------------------------------
    % warning('---- Orthogonalising valXs and valYs\n');
    valXso= valXs;
    valYso= valYs;

    for nn=2:size(valXs,2)
        valXs(:,nn)= orthog( valXso(:,nn), [valXso(:,1:nn-1)]);
        valYs(:,nn)= orthog( valYso(:,nn), [valYso(:,1:nn-1)]);
    end
    valXs= valXs * diag(sign(diag(corr(valXs,valXso)) ));
    valYs= valYs * diag(sign(diag(corr(valYs,valYso)) ));

    % What is the (out-of-sample) correlation b/t datasets i.e. which
    % components/variates are significant
    [r2,p2] = corr(valXs, valYs);
%         Rdw  = diag(r2);
%         Pdw  = diag(p2);
%         Cdw  = diag(cov([valXs,valYs]),numComp);

    if permutePW
        results.Rdw  = diag(r2);
        results.Pdw  = diag(p2);
        results.Cdw  = diag(cov([valXs,valYs]),numComp);
    end

%         % -----------------------------------------------------------------
%         % REGULARISED ASSOCIATION on Out-of-sample scores
%         % -----------------------------------------------------------------
%         cfg     = [];
%         cfg.X   = X;
%         cfg.Y   = Y;
%         cfg.lambdaX  = val{1}.lambdaX; 
%         cfg.lambdaY  = val{1}.lambdaY;
% 
%         for ipc = 1:numComp
%             % Weights of the canonical variates for the Nth X component (u) 
%             % and the Nth Y component (v). Here we need to consider 
%             % average/SVD of weights across folds in the Nth compnents. The
%             % latter might be robust to sign-flips across folds
%             datXw = cell2mat(arrayfun(@(x) x{1}(:,ipc)',Xw,'UniformOutput',0))';
%             [coeff, scoreX, latent, tsquared, explained] = pca(datXw);
%             datYw = cell2mat(arrayfun(@(x) x{1}(:,ipc)',Yw,'UniformOutput',0))';
%             if size(datYw,1)>1 
%                 [coeff, scoreY, latent, tsquared, explained] = pca(datYw);
%             else % MLR mode
%                 scoreY = median(datYw);
%             end
% 
%             cfg.Xw  = scoreX(:,1);
%             cfg.Yw  = scoreY(:,1);
%             [RAoos(ipc,1), NRAoos(ipc,1)] = csa_stats_regularizedAssociation(cfg);
%         end
%         RAoos = abs(RAoos); % PCA on weights may flip signs, so take absolute for the moment. (ideally match with fold-wise RAraw)
%         NRAoos = abs(NRAoos);

    %%% --------------------------------------------
    %%% Out-of-sample loadings for visualization
    %%% --------------------------------------------
    % warning('---- Calculating out-of-sample loadings valXL and valYL in original space\n');
    [valXL,valXLp] = corr(cfg.X,valXs);
    [valYL,valYLp] = corr(Y,valYs);


    %%%% ----------------------------------------------------------------------
    % Fold-split Reproducibilty - Calculate average correlation of
    % Weights acoross folds, i.e. how reproducible are the weights across
    % datasets/fold-wise data. also referred to Half-split Reproducibility 
    % testing as in Churchill et al 2013
    %%% -----------------------------------------------------------------------
    if doSplitHalfFW % Not complete
        RhsfwXL  = [];
        RhsfwYL  = [];
        valXLfw  = [];
        valYLfw  = [];
        datXLte  = cell2mat(XLte');
        datYLte  = cell2mat(YLte');

        RhsfwXW  = [];
        RhsfwYW  = [];
        valXWfw  = [];
        valYWfw  = [];
        datXWte  = cell2mat(XWte');
        datYWte  = cell2mat(YWte');


%             datXLCorr = abs(corr(datXLte)); % All correlations
%             datYLCorr = abs(corr(datYLte)); % All correlations
        for ii = 1:numComp
            % ------------------
            % Loadings
            [~, Xscore,~,~,rvalX]    = pca(datXLte(:,ii:numComp:end),'NumComponents',1);
            [~, Yscore,~,~,rvalY]    = pca(datYLte(:,ii:numComp:end),'NumComponents',1);
            valXLfw(:,ii)   = Xscore;
            valYLfw(:,ii)   = Yscore;
            RhsfwXL(ii,1)   = rvalX(1)/100;%mean(abs(nonzeros(triu(corr(datXLte(:,ii:numComp:end)),1)))); 
            RhsfwYL(ii,1)   = rvalY(1)/100;%mean(abs(nonzeros(triu(corr(datYLte(:,ii:numComp:end)),1)))); 
           
            % ------------------
            % Weights
            [~, Xscore,~,~,rvalX]     = pca(datXWte(:,ii:numComp:end),'NumComponents',1);
            [~, Yscore,~,~,rvalY]     = pca(datYWte(:,ii:numComp:end),'NumComponents',1);
            valXWfw(:,ii)   = Xscore;
            valYWfw(:,ii)   = Yscore;
            RhsfwXW(ii,1)   = rvalX(1)/100;%mean(abs(nonzeros(triu(corr(datXWte(:,ii:numComp:end)),1)))); 
            RhsfwYW(ii,1)   = rvalY(1)/100;%mean(abs(nonzeros(triu(corr(datYWte(:,ii:numComp:end)),1)))); 
            
            dat = datXWte(:,ii:numComp:end);
            Ns = size(dat,1);
            R  = eye(Ns) - Xscore*pinv(Xscore);
            rW = R*dat;
            valXWfwStd(:,ii) = std(rW,0,2);
            
            dat = datYWte(:,ii:numComp:end);
            Ns = size(dat,1);
            R  = eye(Ns) - Yscore*pinv(Yscore);
            rW = R*dat;
            valYWfwStd(:,ii) = std(rW,0,2);        
        end
        % ------------------------------------------------
        % Compute scores with average Weights across folds
        valXSfw = cfg.X * valXWfw;
        valYSfw = cfg.Y * valYWfw;
        
        
        results.RhsfwXL = RhsfwXL;
        results.RhsfwYL = RhsfwYL;

        results.RhsfwXW = RhsfwXW;
        results.RhsfwYW = RhsfwYW;
        
        results.RhsfwS  = diag(corr(valXSfw,valYSfw));
       
    end


    %%% -----------------------------------
    %%% Assemble the data in CCA structure
    %%% -----------------------------------

    % Create variable names and assign all outputs to 'results' structure
%         results     = [];
%         results.Rfw = Rfw;
%         results.Pfw = Pfw;
%         results.Cfw = Cfw;

%         if permutePW
%             results.Rdw = Rdw;
%             results.Pdw = Pdw;
%             results.Cdw = Cdw;
%         end
%         if doSplitHalfFW
%             results.RsplitXL = RsplitXL;
%             results.RsplitYL = RsplitYL;
%         end
%         if doSplitHalfDW
%             varNames = [varNames 'XwCorrR','XwCorrP','YwCorrR','YwCorrP' ];
%             results.XwCorrR = XwCorrR;
%             results.XwCorrP = XwCorrP;
%             results.YwCorrR = YwCorrR;
%             results.YwCorrP = YwCorrP;
%         end


    if iPerm ==1
%             CCAtemp       = [];
%             CCAtemp(iPerm).CVA    = CVA;
%             CCAtemp(iPerm).opt    = opt;
%             CCAtemp(iPerm).val    = val;
        % Out-of-sample correlation b/t datasets i.e. which components/variates are significant

        % Out-of-sample loadings for visualization
        CCAtemp(iPerm).valXL  = valXL;
        CCAtemp(iPerm).valXLp = valXLp;
        CCAtemp(iPerm).valYL  = valYL;
        CCAtemp(iPerm).valYLp = valYLp;
        CCAtemp(iPerm).valXs  = valXs;
        CCAtemp(iPerm).valYs  = valYs;
%         CCAtemp(iPerm).varNames = varNames;

        try
            CCAtemp(iPerm).valXLhsfw = valXLfw;
            CCAtemp(iPerm).valYLhsfw = valYLfw;
            
            CCAtemp(iPerm).valXWhsfw = valXWfw;
            CCAtemp(iPerm).valYWhsfw = valYWfw;
            
            CCAtemp(iPerm).valXShsfw   = valXSfw;
            CCAtemp(iPerm).valYShsfw   = valYSfw;
            CCAtemp(iPerm).valXWhsfwStd = valXWfwStd;
            CCAtemp(iPerm).valYWhsfwStd = valYWfwStd;
        end


        for ivar = 1:numel(varNames)
            varName = varNames{ivar};

%                  try
%                     eval(sprintf('CCAtemp.%s = median(%s,2);',varNames{ivar},varNames{ivar}));
                CCAtemp(iPerm).(varName) = median(results.(varName),2); 
%                 end
%                     eval(sprintf('CCAtemp.null%s = [];',varNames{ivar}));
%                 eval(sprintf('CCAtemp.pPerm%s = [];',varNames{ivar}));
            varName = ['pPerm' varName];
            CCAtemp(iPerm).(varName) = [];
%                 
        end
    else
        for ivar = 1:numel(varNames)
            varName = varNames{ivar};
%                 eval(sprintf('null%s(iPerm,:) = mean(%s,2);',varNames{ivar},varNames{ivar}));
            varNull = ['null' varName];
            resultsNull(iPerm).(varNull) = median(results.(varName),2);

%                 if doSavePermParams
%                     resultsNull(iPerm).valXLfw = valXLfw;
%                     resultsNull(iPerm).valXL   = valXL;
%                 end
        end
    end

    % -----------------------------------------------------------------
    % Save Permutation-based parameters, e.g. Loadings, Weights etc.  
    if doSaveNullsFW
        fout = fullfile(cfg.dirOut,sprintf('nullResults_perm%06d.mat',iPerm));
        XLfw = valXLfw;
        XWfw = valXWfw;
        XL   = valXL;
        XS   = valXs;

        kat_parfor_save(fout,{'XL','XS','XLfw','XWfw','results'},XL,XS,XLfw,XWfw,results);
    end
%         Y = [];
%         opt = [];
%         RsplitYL = [];
%         RsplitXL = [];
end % fold-wise permutations