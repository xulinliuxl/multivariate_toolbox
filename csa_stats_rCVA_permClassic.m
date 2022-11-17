function [CCA,CCApart,CCAperm] = csa_stats_rCVA_permClassic(CVA)
% ----------
% Unpack CVA
% ----------

% varnames = fieldnames(CVA);
% for iVar = 1:numel(varnames)
%     varname = varnames{iVar};
%     eval(sprintf('%s = CVA.%s;',varname,varname));
% end

cfg = CVA.mode.permClassic;

numPerm = cfg.numPerm;
Xorig   = CVA.X;% Need this for the permutations
Yorig   = CVA.Y;
Ns      = CVA.Ns;
numComp = CVA.numComp;
doSaveNulls = CVA.doSaveNulls;

if ~isempty(CVA.dirOut)
    dirOut = CVA.dirOut;
    
%     CVA.dirOut = dirOut;
end

varnames = {
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

% Predefine order for outter permutations
randOrderOut = palm_quickperms(Ns,[],numPerm);


if CVA.runInSerial
    parforArg = 0;
else
    parforArg = Inf;
end


opt             = CVA; 
opt.svd         = 0;
opt.verbose     = 0;
% opt.numPermIn   = numPerm;
opt.varnames    = varnames;
opt.parforArg   = parforArg;
opt.numComp     = numComp;
CCA             = csa_proc_CCA(opt,Xorig,Yorig);

XLabs   = abs(CCA.XL);
YLabs   = abs(CCA.YL);
XLsum   = zeros(CVA.numVarX,numComp);
YLsum   = zeros(CVA.numVarY,numComp);
% permXL  = zeros(CVA.numVarX,numComp,numPerm); 

% permYL
dirOutNulls = [];
if doSaveNulls
    dirOutNulls = CVA.dirOutNulls;
    fout = fullfile(dirOutNulls,sprintf('nullResults_perm_%06d.mat',1));
    XL = CCA.XL;
    XS = CCA.XS;
%     kat_save_parfor(fout,XLperm,'XL');
    kat_parfor_save(fout,{'XL','XS'},XL,XS);
end

% pctRunOnAll warning('off','all')
% ppm = ParforProgMon('Progress: ', numPerm, max((numPerm)/100,1));

parfor (iPermOut = 2:numPerm, parforArg)
%     for iPermOut = 1:numPerm % Permutation of subject labels. First iteration uses the correct labels
    
%     if iPermOut == 1    
%         resultsIn = resultsCCA;
%     else
        % Get ordering for X,Y and Covs
        tempOrderOut    = randOrderOut(:,iPermOut);
        Y               = Yorig(tempOrderOut,:);
        [CCAiperm] = csa_proc_CCA(opt,Xorig,Y);
        
       

        
        % Write variables to resultsOut for later
        resultsOut(iPermOut).RvalMedian = abs(CCAiperm.R);
        resultsOut(iPermOut).CvalMedian = abs(CCAiperm.Cov);
        
        if doSaveNulls
            fout = fullfile(dirOutNulls,sprintf('nullResults_perm_%06d.mat',iPermOut));
            XL = CCAiperm.XL;
            XS = CCAiperm.XS;
%             kat_save_parfor(fout,XLperm,'XLperm');
            kat_parfor_save(fout,{'XL','XS'},XL,XS);
            
             % ----------------------------------------------------
            % Are permuted Loadings larger than observed loadings?
            % Used later to compute the significance of loadings based on
            % permutations.
            % Here using abs statistic, since the sign of the loadings is
            % artibitrary. 
            % ----------------------------------------------------
            xl = XLabs < abs(CCAiperm.XL);
            XLsum = XLsum + xl;
            yl = YLabs < abs(CCAiperm.YL);
            YLsum = YLsum + yl;
            
        end
%     end
    
%     for ivar = 1:numel(varnames)
%         varname = varnames{ivar};
% %             resultsOut(iPermOut).(varname) = resultsIn.(varname);%nanmedian(resultsIn.(varname),2);  % Use nanmedian as in some cases SVD results in very small (near zero) eigenvalues, which result in nans
%         resultsOut(iPermOut).([varname 'Median']) = nanmedian(resultsIn.(varname),2);%nanmedian(resultsIn.(varname),2);  % Use nanmedian as in some cases SVD results in very small (near zero) eigenvalues, which result in nans
%     end



%     if (mod(iPermOut, 10) == 0)
%         ppm.increment();
%     end
end % iPermOut

% pctRunOnAll warning('on','all'); 
% clear ppm; 

resultsOut(1).RvalMedian = abs(CCA.R);
resultsOut(1).CvalMedian = abs(CCA.Cov);

% varnames = [varnames 'Rval' 'Cval']; % Add Rval and Cval to variable names
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

    % N.B. Hack!!! Check for NaNs in the data
%     idxnan = isnan(nulldata);
%     idxnanEnd1 = idxnan(end-1,:);
%     nulldata(end-1,idxnanEnd1) = min(nulldata(:,idxnanEnd1));
%     idxnanEnd = idxnan(end,:);
%     nulldata(end,idxnanEnd) = min(nulldata(:,idxnanEnd));


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
if doSaveNulls
    CCA.XLp   = XLsum  ./(numPerm-1);
    CCA.YLp   = YLsum  ./(numPerm-1);
end

if CVA.doAdjustDOF % based on canoncor controlling for DFs
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

