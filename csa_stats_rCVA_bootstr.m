function [CCA,Results] = csa_stats_rCVA_bootstr(opt,X,Y)
% Fonction for Inner Bootstrapping of Split-half permuations

% Where X and Y has been provided, do not consider opt.X and opt.Y
if ~exist('X','var')
    X = opt.X;
    Y = opt.Y;
end


CCA     = csa_proc_CCA(opt,X,Y);
CCA.Cov = diag(cov([CCA.XS CCA.YS]),opt.numComp);%diag(CCA.XS' * CCA.YS); % Calculate covariance between subject scores

numComp    = opt.numComp;
numBoot    = opt.numBoot;
idxPermIn1 = opt.idxPermIn1;
idxPermIn2 = opt.idxPermIn2;
varnames   = opt.varnames;

% Define variables in inner permutation
for ivar = 1:numel(varnames)
    varname = varnames{ivar};
    resultsIn.(varname) = nan(numComp,numBoot); 
end  

XlSS = zeros(size(X,2),numComp);
YlSS = zeros(size(Y,2),numComp);
XlSR = zeros(size(X,2),numComp);
YlSR = zeros(size(Y,2),numComp);

% ccaboot = cell(numBoot,1);

%         cv_label = cvpartition(labelIdx,'kFold',2);  %% partitions with equal numbersdrawn from each decile

% -----------------------------------------------------------------
% Inner Permutations
% -----------------------------------------------------------------
% for iboot = 1:numBoot
parfor (iboot = 1:numBoot, opt.parforArg) 
    
    iCCA = CCA;
%             cv_label= repartition(cv_label); % Could replace with palm perms, but here we ensure equal representation of each group (labelIdx), e.g. age, or gender
    idx_1   = idxPermIn1(:,iboot);%cv_label.training(1); % Index for subjects in the first half split
    idx_2   = idxPermIn2(:,iboot);%cv_label.test(1);     % Index for subjects in the second half split 

    % Select subjects in each split
    iCCA.X1  = iCCA.X(idx_1,:);
    iCCA.X2  = iCCA.X(idx_2,:);
    iCCA.Y1  = iCCA.Y(idx_1,:);
    iCCA.Y2  = iCCA.Y(idx_2,:);

    iCCA.Xs1 = iCCA.XS(idx_1,:);
    iCCA.Xs2 = iCCA.XS(idx_2,:);
    iCCA.Ys1 = iCCA.YS(idx_1,:);
    iCCA.Ys2 = iCCA.YS(idx_2,:);

%                 CCA.Xw1 = CCA.XW(idx_1,:);
%                 CCA.Xw2 = CCA.XW(idx_2,:);
%                 CCA.Yw1 = CCA.YW(idx_1,:);
%                 CCA.Yw2 = CCA.YW(idx_2,:);

    % ---------------------------------------------------------
    % Compute Split-half Stability (SS) Kovacevic et al 2013
    % ---------------------------------------------------------
    % Correlations b/n projected left and right split-half patterns
    % These are taken as measres of the correspondence b/n X data
    % and Y weights on one hand, and Y data and X weights, on the
    % other hand. 
    % N.B. Unlike Split-half Reliability (below), in split-half stability
    % the association is not based on the ABSOLUTE correlation, since W
    % and L is expected to be the same across splits.

%                 CCA.Xw1 = corr(CCA.X1,CCA.Ys1);% Can replace Y1 * Yw with Ys for samples from this split
%                 CCA.Xw2 = corr(CCA.X2,CCA.Ys2);
%                 CCA.Yw1 = corr(CCA.Y1,CCA.Xs1);
%                 CCA.Yw2 = corr(CCA.Y2,CCA.Xs2);
%                                 
%                 CCA.Xw1 = CCA.X1' * CCA.Ys1;% Can replace Y1 * Yw with Ys for samples from this split
%                 CCA.Xw2 = CCA.X2' * CCA.Ys2;
%                 CCA.Yw1 = CCA.Y1' * CCA.Xs1;
%                 CCA.Yw2 = CCA.Y2' * CCA.Xs2;
%                
    % Use Weights, which might be more sensitive for low
    % correlation and high covariance
    iCCA.Xw1 = iCCA.X1' * iCCA.Y1 * iCCA.YW;% Can replace Y1 * Yw with Ys for samples from this split
    iCCA.Xw2 = iCCA.X2' * iCCA.Y2 * iCCA.YW;
    iCCA.Yw1 = iCCA.Y1' * iCCA.X1 * iCCA.XW;
    iCCA.Yw2 = iCCA.Y2' * iCCA.X2 * iCCA.XW;

%             CCA.yw1 = CCA.Y1' * CCA.Y1 * CCA.YW;% Can replace Y1 * Yw with Ys for samples from this split
%             CCA.yw2 = CCA.Y2' * CCA.Y2 * CCA.YW;
%             CCA.xw1 = CCA.X1' * CCA.X1 * CCA.XW;
%             CCA.xw2 = CCA.X2' * CCA.X2 * CCA.XW;
% 
%             
%             yw(:,iboot) = diag(corr( CCA.yw1, CCA.yw2));
%             xw(:,iboot) = diag(corr( CCA.xw1, CCA.xw2));
%         end
    [r,p] = corr(iCCA.Xw1,iCCA.Xw2);
    resultsIn(iboot).XwSSr = diag(r); % Component by interation array
    resultsIn(iboot).XwSSp = diag(p);
    resultsIn(iboot).XwSSc = diag(cov([iCCA.Xw1 iCCA.Xw2]),numComp);

    [r,p] = corr(iCCA.Yw1,iCCA.Yw2);
    resultsIn(iboot).YwSSr = diag(r); % Component by iteration array
    resultsIn(iboot).YwSSp = diag(p);  
    resultsIn(iboot).YwSSc = diag(cov([iCCA.Yw1 iCCA.Yw2]),numComp);

    % Weights as in Kovacevic paper, which are the same as
    % above
%                 cm1=CCA.X1' * CCA.Y1; % Covariance in split 1
%                 cm2=CCA.X2' * CCA.Y2; % Covariance in split 2
%                 CCA.Xw1k = cm1 * CCA.YW;
%                 CCA.Xw2k = cm2 * CCA.YW;
%                 CCA.Yw1k = cm1' * CCA.XW;
%                 CCA.Yw2k = cm2' * CCA.XW;
%                 [r1,p1] = corr(CCA.Xw1k,CCA.Xw2k);
%                 [r1,p1] = corr(CCA.Xw1k,CCA.Xw2k);


    % Loadings
    iCCA.Xl1 = iCCA.X1' * iCCA.X1 * iCCA.Xw1; % Xs = X * Xw; XL = X' * Xs;
    iCCA.Xl2 = iCCA.X2' * iCCA.X2 * iCCA.Xw2;
    iCCA.Yl1 = iCCA.Y1' * iCCA.Y1 * iCCA.Yw1;
    iCCA.Yl2 = iCCA.Y2' * iCCA.Y2 * iCCA.Yw2;

    [r,p] = corr(iCCA.Xl1,iCCA.Xl2);
    resultsIn(iboot).XlSSr = diag(r); % Component by iteration array
    resultsIn(iboot).XlSSp = diag(p);
    resultsIn(iboot).XlSSc = diag(cov([iCCA.Xl1 iCCA.Xl2]),numComp);


    [r,p] = corr(iCCA.Yl1,iCCA.Yl2);
    resultsIn(iboot).YlSSr = diag(r); % Component by iteration array
    resultsIn(iboot).YlSSp = diag(p);  
    resultsIn(iboot).YlSSc = diag(cov([iCCA.Yl1 iCCA.Yl2]),numComp);


    % Sum Loadings across permutations to estimate average Loadings
    % afterwards
    tempdat = (iCCA.Xl1+iCCA.Xl2)./2;
    signdat = sign(diag(corr(CCA.XL,tempdat)));
    tempdat = tempdat .* signdat';
    XlSS    = XlSS + tempdat;

    tempdat = (iCCA.Yl1+iCCA.Yl2)./2;
    signdat = sign(diag(corr(CCA.YL,tempdat)));
    tempdat = tempdat .* signdat';
    YlSS    = YlSS + tempdat;


    % ---------------------------------------------------------
    % Compute Split-half Reproducibility (SR) Churchill et al
    % 2013
    % Compares Ws (or Ls) estimated on independent data split-halves
    % ---------------------------------------------------------
    [CCA1]= csa_proc_CCA(opt,iCCA.X1,iCCA.Y1);
    [CCA2]= csa_proc_CCA(opt,iCCA.X2,iCCA.Y2);

    [r,p] = corr(CCA1.XW,CCA2.XW);
    resultsIn(iboot).XwSRr = abs(diag(r));
    resultsIn(iboot).XwSRp = diag(p);
    resultsIn(iboot).XwSRc = abs(diag(cov([CCA1.XW CCA2.XW]),numComp));%                 XwSRc(:,iboot) = diag(CCA1.XW' * CCA2.XW);

    [r,p] = corr(CCA1.YW,CCA2.YW);
    resultsIn(iboot).YwSRr = abs(diag(r));
    resultsIn(iboot).YwSRp = diag(p);
    resultsIn(iboot).YwSRc = abs(diag(cov([CCA1.YW CCA2.YW]),numComp));%                 YwSRc(:,iboot) = diag(CCA1.YW' * CCA2.YW);

    [r,p] = corr(CCA1.XL,CCA2.XL);
    resultsIn(iboot).XlSRr = abs(diag(r));
    resultsIn(iboot).XlSRp = diag(p);
    resultsIn(iboot).XlSRc = abs(diag(cov([CCA1.XL CCA2.XL]),numComp));%                 XlSRc(:,iboot) = diag(CCA1.XL' * CCA2.XL);

    [r,p] = corr(CCA1.YL,CCA2.YL);
    resultsIn(iboot).YlSRr = abs(diag(r));
    resultsIn(iboot).YlSRp = diag(p);   
    resultsIn(iboot).YlSRc = abs(diag(cov([CCA1.YL CCA2.YL]),numComp));%                 YlSRc(:,iboot) = diag(CCA1.YL' * CCA2.YL);


    % Sum Loadings across permutations to estimate average Loadings
    % afterwards
    signdat1 = sign(diag(corr(CCA.XL,CCA1.XL)));
    signdat2 = sign(diag(corr(CCA.XL,CCA2.XL)));
    tempdat = ((CCA1.XL .* signdat1') + (CCA2.XL .* signdat2'))./2;
    XlSR    = XlSR + tempdat;

    signdat1 = sign(diag(corr(CCA.YL,CCA1.YL)));
    signdat2 = sign(diag(corr(CCA.YL,CCA2.YL)));
    tempdat = ((CCA1.YL .* signdat1') + (CCA2.YL .* signdat2'))./2;
    YlSR    = YlSR + tempdat;
    
%     
%     ccaboot{iboot}.XLsr1 = CCA1.XL;
%     ccaboot{iboot}.YLsr1 = CCA1.YL;
%     ccaboot{iboot}.XWsr1 = CCA1.XW;
%     ccaboot{iboot}.YWsr1 = CCA1.YW;
%     
%     ccaboot{iboot}.XLsr2 = CCA2.XL;
%     ccaboot{iboot}.YLsr2 = CCA2.YL;
%     ccaboot{iboot}.XWsr2 = CCA2.XW;
%     ccaboot{iboot}.YWsr2 = CCA2.YW;
%     
%     ccaboot{iboot}.XLss1 = iCCA.Xl1;
%     ccaboot{iboot}.YLss1 = iCCA.Yl1;
%     ccaboot{iboot}.XWss1 = iCCA.Xw1;
%     ccaboot{iboot}.YWss1 = iCCA.Yw1;
%     
%     ccaboot{iboot}.XLss2 = iCCA.Xl2;
%     ccaboot{iboot}.YLss2 = iCCA.Yl2;
%     ccaboot{iboot}.XWss2 = iCCA.Xw2;
%     ccaboot{iboot}.YWss2 = iCCA.Yw2;
    
    iCCA = [];
end % iboot 

for ivar = 1:numel(varnames)
    varname = string(varnames(ivar));
    Results.(varname) = [resultsIn.(varname)];    
end

% Get mean SS and SR Loadings across all splits
CCA.XLss = XlSS./numBoot;
CCA.XLsr = XlSR./numBoot;

CCA.YLss = YlSS./numBoot;
CCA.YLsr = YlSR./numBoot;





