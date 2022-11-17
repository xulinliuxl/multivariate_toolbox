function [CCA,CCApart,CCAperm] = csa_stats_rCVA_permBootstr(CVA)
 % Half-split permutations (Strother et al 2002 NI and Kovacevic et al 2013)
   
% ----------
% Unpack CVA
% ----------

% varnames = fieldnames(CVA);
% for iVar = 1:numel(varnames)
%     varname = varnames{iVar};
%     eval(sprintf('%s = CVA.%s;',varname,varname));
% end

cfg = CVA.mode.permBootstr;

% -------------------------------------------------------------------------
% Unpack CVA
numPerm             = cfg.numPerm;
labels              = cfg.labels;
Xorig               = CVA.X;% Need this for the permutations
Yorig               = CVA.Y;
Ns                  = CVA.Ns;
numComp             = CVA.numComp;
dirOut              = CVA.dirOutNulls;
doSaveNulls         = CVA.doSaveNulls;
presetRandOrder     = CVA.randOrder;
usePresetRandOrder  = CVA.usePresetRandOrder;

% Labelling samples to ensure partitions with equal numbersdrawn from each group/decile
if isempty(labels)
    labels  = ones(Ns,1);
end




varnames = {'XwSSr','YwSSr','XwSSp','YwSSp','XwSSc','YwSSc',...
            'XlSSr','YlSSr','XlSSp','YlSSp','XlSSc','YlSSc',...,
            'XwSRr','YwSRr','XwSRp','YwSRp','XwSRc','YwSRc',...
            'XlSRr','YlSRr','XlSRp','YlSRp','XlSRc','YlSRc',...,
            'Rval','Cval',...
            };

% Predifine variabels for outer permutations
for ivar = 1:numel(varnames)%         
    eval(sprintf('resultsOut(numPerm).%sMedian = nan(numComp,1);',varnames{ivar}));
    eval(sprintf('resultsOut(numPerm).%s       = nan(numComp,numPerm);',varnames{ivar}));
    ...eval(sprintf('resultsOutput.%sMedian = nan(numComp,(numPerm));',varnames{ivar}));    
    eval(sprintf('resultsOutput.pPerm%s = nan(numComp,1);',varnames{ivar}));
    eval(sprintf('resultsOutput.pPareto%s = nan(numComp,1);',varnames{ivar}));
end       
varnames(end-1:end) = []; % Remove Rval and Cval 


% -------------------------------------------------------------------------
% Predefine order for outter and inner permutations
randOrderOut = [];
if usePresetRandOrder
%         randIdx     = randperm(size(presetRandOrder,2),numPermFW);
    randOrderOut = presetRandOrder(:,1:numPerm);
    numMid     = round(Ns/2);
    idxPermIn1 = randOrderOut(1:numMid,:);
    idxPermIn2 = randOrderOut(numMid+1:end,:);
else
    randOrderOut = palm_quickperms(Ns,[],numPerm);
    % Predifine order for inner permutations so that its kept the same for
    % all outer permutations
    for iP = 1:numPerm
        cv_label = cvpartition(labels,'kFold',2);  %% partitions with equal numbersdrawn from each decile
        idxPermIn1(:,iP) = cv_label.training(1); % Index for subjects in the first half split
        idxPermIn2(:,iP) = cv_label.test(1);
    end
end



if CVA.runInSerial
    parforArg = 0;
else
    parforArg = Inf;
end


opt             = CVA; 
opt.svd         = 0;
opt.verbose     = 0;
opt.numBoot     = numPerm;
opt.varnames    = varnames;
opt.idxPermIn1  = idxPermIn1;    
opt.idxPermIn2  = idxPermIn2;
opt.parforArg   = parforArg;
[CCA,resultsCCA] = csa_stats_rCVA_bootstr(opt,Xorig,Yorig);


% Set numPermIn to a lower number for inner loop, since parpool is on
% the outer loop
numPermIn = numPerm/2;
opt.numBoot = 100;%numPermIn;

XLabs   = abs(CCA.XL);
XLssAbs = abs(CCA.XLss);
XLsrAbs = abs(CCA.XLsr);
YLabs   = abs(CCA.YL);
YLssAbs = abs(CCA.YLss);
YLsrAbs = abs(CCA.YLsr);
XL      = zeros(CVA.numVarX,numComp);
XLss    = XL;
XLsr    = XL;
YL      = zeros(CVA.numVarY,numComp);
YLss    = YL;
YLsr    = YL;

% parfor (iPermOut = 1:numPerm, parforArg)
for iPerm = 1:numPerm % Permutation of subject labels. First iteration uses the correct labels
    
    if iPerm == 1    
        resultsIn = resultsCCA;
    else
        % Get ordering for X,Y and Covs
        tempOrderOut    = randOrderOut(:,iPerm);
        Y               = Yorig(tempOrderOut,:);
        [C,resultsIn]   = csa_stats_rCVA_bootstr(opt,Xorig,Y);
        
        % ----------------------------------------------------
        % Are permuted Loadings larger than observed loadings?
        % Used later to compute the significance of loadings based on
        % permutations.
        % ----------------------------------------------------
        xl = XLabs < abs(C.XL);
        XL = XL + xl;
        yl = YLabs < abs(C.YL);
        YL = YL + yl;

        xl   = XLssAbs < abs(C.XLss);
        XLss = XLss + xl;
        yl   = YLssAbs < abs(C.YLss);
        YLss = YLss + yl;

        xl   = XLsrAbs < abs(C.XLsr);
        XLsr = XLsr + xl;
        yl   = YLsrAbs < abs(C.YLsr);
        YLsr = YLsr + yl;
        
        % Write variables to resultsOut for later
        resultsOut(iPerm).RvalMedian = C.R;
        resultsOut(iPerm).CvalMedian = C.Cov;
        
         % Save null results
        if doSaveNulls
           fout = fullfile(dirOut,sprintf('nullResults_perm%06d.mat',iPerm));
           kat_parfor_save(fout,{'XLss','YLss','XLsr','YLsr',},C.XLss,C.YLss,C.XLsr,C.YLsr);
        end
        
    end
    
    for ivar = 1:numel(varnames)
        varname = varnames{ivar};
%             resultsOut(iPermOut).(varname) = resultsIn.(varname);%nanmedian(resultsIn.(varname),2);  % Use nanmedian as in some cases SVD results in very small (near zero) eigenvalues, which result in nans
        resultsOut(iPerm).([varname 'Median']) = nanmedian(resultsIn.(varname),2);%nanmedian(resultsIn.(varname),2);  % Use nanmedian as in some cases SVD results in very small (near zero) eigenvalues, which result in nans
    end
    
end % iPermOut
    

resultsOut(1).RvalMedian = CCA.R;
resultsOut(1).CvalMedian = CCA.Cov;

% Save Permutation-based parameters, e.g. Loadings, Weights etc.  


varnames = [varnames 'Rval' 'Cval']; % Add Rval and Cval to variable names
for ivar = 1:numel(varnames) % calculate stats based on null distribution

    varname = varnames{ivar};
    varnameMedian = [varname 'Median'];
    varnamePerm = ['pPerm' varname];
    varnamePareto = ['pPareto' varname];
    origdata    = resultsOut(1).(varnameMedian);
    nulldata    = [resultsOut(2:end).(varnameMedian)];

%         % significance based on identican splits across permutations
%         alldata = reshape([resultsOut.(varname)],numComp,numPermIn,numPerm);
%         alldata = permute(alldata, [3 1 2]);
%         alldata = alldata(:,:);
%         Xdat = repmat(squeeze(alldata(2,:,1))',1,numPerm-1);
%         Ydat = squeeze(alldata(2,:,2:end));
%         
%         [H,P,CI,STATS] = ttest(Xdat,Ydat);

    % N.B. Hack!!! Check for NaNs in data
    idxnan = isnan(nulldata);
    idxnanEnd1 = idxnan(end-1,:);
    nulldata(end-1,idxnanEnd1) = min(nulldata(:,idxnanEnd1));
    idxnanEnd = idxnan(end,:);
    nulldata(end,idxnanEnd) = min(nulldata(:,idxnanEnd));


    if strcmp(varnames{ivar}(end),'p') % If p-values, look for tests having lower value than null data
        resultsOutput.(varnamePerm) = sum(nulldata' < repmat(origdata',numPerm-1,1))/numPerm; 
        resultsOutput.(varnamePareto) = palm_pareto(origdata',nulldata',1,0,0);

%             a=reshape(sum(alldata(2:end,:) < repmat(alldata(1,:),numPerm-1,1))/numPerm,numComp,numPerm); 


%             eval(sprintf('pPerm%s = sum(null%sMedian'' < repmat(%sMedian'',numPerm-1,1))/numPerm;',varnames{ivar},varnames{ivar},varnames{ivar})); 
%             eval(sprintf('pPareto%s = palm_pareto(%sMedian'',null%sMedian'',1,0,0);',varnames{ivar},varnames{ivar},varnames{ivar})); 
    else % If r-values (or other non-pvalues), look for tests having higher value than null data
        resultsOutput.(varnamePerm) = sum(nulldata' > repmat(origdata',numPerm-1,1))/numPerm; 
        resultsOutput.(varnamePareto) = palm_pareto(origdata',nulldata',0,0,0);

%             a=reshape(sum(alldata(2:end,:) > repmat(alldata(1,:),numPerm-1,1))/numPerm,numComp,numPerm); 

%             eval(sprintf('pPerm%s = sum(null%sMedian'' > repmat(%sMedian'',numPerm-1,1))/numPerm;',varnames{ivar},varnames{ivar},varnames{ivar})); 
%             eval(sprintf('pPareto%s = palm_pareto(%sMedian'',null%sMedian'',0,0,0);',varnames{ivar},varnames{ivar},varnames{ivar})); 
    end
end

% try
for ivar = 1:numel(varnames)
   eval(sprintf('CCA.pPareto%s = resultsOutput.pPareto%s;',varnames{ivar},varnames{ivar})); 
   eval(sprintf('CCA.pPerm%s = resultsOutput.pPerm%s;',varnames{ivar},varnames{ivar})); 
%        eval(sprintf('CVA.null%sMedian = null%sMedian;',varnames{ivar},varnames{ivar})); 
%    eval(sprintf('CVA.%sMedian = %sMedian;',varnames{ivar},varnames{ivar})); 
end
% end

% % Estimate significance of Loadings based on permutations
CCA.XLp   = XL  ./(numPerm-1);
CCA.XLpSS = XLss./(numPerm-1);
CCA.XLpSR = XLsr./(numPerm-1);
CCA.YLp   = YL  ./(numPerm-1);
CCA.YLpSS = YLss./(numPerm-1);
CCA.YLpSR = YLsr./(numPerm-1);

doAdjustDOF = 0;
if doAdjustDOF % based on canoncor controlling for DFs
    R = CVA.r';
    [n,p1] = size(X);   
    [Q1,T11,perm1] = qr(X,0);
    rankX = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));

    p2 = size(Y,2);
    [Q2,T22,perm2] = qr(Y,0);
    rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));

    d = min(rankX,rankY);
    k = 0:(d-1);
    d1k = (rankX-k);
    d2k = (rankY-k);
    nondegen = find(R < 1);
    logLambda = -Inf( 1, d);
    logLambda(nondegen) = cumsum(log(1 - R(nondegen).^2), 'reverse');

    CCA.df1 = d1k .* d2k;

    % Lawley's modification to Bartlett's chi-squared statistic
    CCA.chisq = -(n - k - .5*(rankX+rankY+3) + cumsum([0 1./R(1:(d-1))].^2)) .* logLambda(1:d); % kat added '(1:d)' to take the first d R-values 
    a = d1k .* d2k;
    c = -(n - k - .5*(rankX+rankY+3) + cumsum([0 1./R(1:(d-1))].^2)) .* logLambda(1:d); % kat added '(1:d)' to take the first d R-values 
    CCA.pChisq = csa_stats_chi2pval(CCA.chisq, CCA.df1);
end

