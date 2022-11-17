function [CCA]=  csa_proc_CCA(opt,X,Y)
% This is based on mt_proc_CCA
% Performs a regularized canonical correlation analysis (rCCA).
%Regularization is implemented as a convex combination the empirical
%covariance matrix and an identity matrix:
%C_reg = (1-lambda) * C_orig + lambda * trace(C_orig) * I
%
%Regularization is implemented separately for
%
%Usage:
% mt_proc_CCA(opt,X,Y)
%
%Input
% opt - configuration struct (see below)
% X   - [N x P] matrix with N observations of P variables
% Y   - [N x Q] matrix with N observations of Q variables
%
% More options:
% .lambdaX   - regularization paramter in [0,1] for X variable (default 0
%              = no regularization). Set to 'auto' to use automatic
%              Ledoit-Wolf regularisation.
% .lambdaY   - regularization paramter in [0,1] for Y variable (default 0
%              = no regularization). Set to 'auto' to use automatic
%              Ledoit-Wolf regularisation.
%              We write the regularized covariance matrix as a convex combination of
% the empirical covariance matrix Cxx and the identity matrix I. To bring
% the identity matrix on equal footing with Cxx, we multiply I by
% trace(Cxx).
% The regularization parameters also be vectors. In this case the return
% arguments are cell arrays of size [numel(lambdaX), numel(lambdaY)].
%
% .nComp      - number of CCA components, default min(P,Q)
% .verbose    - if 1 prints some extra output (default 0)
%
% .valX , .valY - if unseen validation [=test] data is given [M x P] resp. [M x Q] then the
%                 canonical variates for the validation data and the canonical
%                 correlation are calculated. Note that this only works
%                 when lambdaX and lambdaY are scalars.
% .covariates   - [N x C] matrix of C covariates. The data is
%                 orthogonalised wrt the covariates after centereing and
%                 before running the CCA [ *** TODO ***]
%Output:
% [Xw, Yw, R, Xs, Ys, more] with
%
% Xw,Yw - canonical weights specifying linear combinations to obtain the
%         latent canonical variates
% R     - vector with canonical correlations
% Xs,Ys - canonical variate scores, or representation of the latent
%         variables in CCA space
%
% more  - struct with extra fields
% XL,YL - canonical loadings of the latent variable on the original
%         variables
% .valXs, valYs - canonical variates for the provided validation data
% .valR   - out-of-sample canonical correlations for the provided validation data
% .valP   - corresponding p-values 
% .lambdaX, lambdaY
%
%Reference:
% HD Vinod. "Canonical ridge and econometrics of joint production".
% Journal of Econometrics, Volume 4, Issue 2, May 1976, Pages 147-166
%
% ------------------------------
% Edit Log
% -----------------------
% 2018-12-21 - Ensure the covariance matrix after chol is symmetric and
%              positive definite using nearestSPD
% 2021-01-21 - Variation of csa_proc_CCA_fast where outputs are assembled
%              at the bottom of the function


% CCA = opt;

N = size(X,1);
P = size(X,2);
Q = size(Y,2);
lambdaX = opt.lambdaX;
lambdaY = opt.lambdaY;

% Do we need mt_set
% csa_setDefault(opt,'nComp', min( [N,P,Q] ));
% csa_setDefault(opt,'lambdaX', 0 );
% csa_setDefault(opt,'lambdaY', 0 );
% csa_setDefault(opt,'Xte', [] );
% csa_setDefault(opt,'Yte', []);
% csa_setDefault(opt,'verbose', 0);
% csa_setDefault(opt,'svd', 0);
% csa_setDefault(opt,'nSvd', max([N,P,Q]));

% CCA.nComp = opt.numComp;

%% Regress covariates of no interest (if requested)
% --------------------------------
% Regress Covariates out from data
% --------------------------------
if opt.doCovs 
    
     if isfield(opt,'doCovFW')
            % covary_out_in_each_fold
        COV_X = CVA.Covs;     % Covariates
        COV_Y = COV_X; 

        COV_Xtrain = COV_X(train_idx,:);
        COV_Xtest  = COV_X(test_idx,:);
        COV_Ytrain = COV_Y(train_idx,:);
        COV_Ytest  = COV_Y(test_idx,:);

        beta_x = nan( size(COV_X,2), size(Xtest,2));
        beta_y = nan( size(COV_Y,2), size(Ytest,2));
        for xx=1:size(Xtest,2)
            beta_x(:,xx)= regress( Xtrain(:,xx), COV_Xtrain);
        end
        for yy=1:size(Ytest,2)
            beta_y(:,yy)= regress( Ytrain(:,yy), COV_Ytrain);
        end

        % Regress out from training data
        Xtrain= Xtrain - COV_Xtrain * beta_x;
        Ytrain= Ytrain - COV_Ytrain * beta_y;

        % Regress out from test data
        Xtest= Xtest - COV_Xtest * beta_x;
        Ytest= Ytest - COV_Ytest * beta_y;
            
    else
        X0 = opt.Covs;
        R  = eye(opt.Ns) - X0*pinv(X0);
        X  = R*X;
    end

%             % In Permutated variantions Y has different ordering and so
%             % Covs needs reordering too.  
%             Y0 = CVA.Covs(tempOrderOut,:);
%             R  = eye(Ns) - Y0*pinv(Y0);
%             Y  = R*Y;
%             

%                 beta_x = nan( size(CVA.Covs,2), size(X,2));
%                 beta_y = nan( size(CVA.Covs,2), size(Y,2));
%                 for xx=1:size(X,2)
%                     beta_x(:,xx) = regress(X(:,xx),CVA.Covs);
%                 end
%                 Xr = X - (CVA.Covs * beta_x);%                 
%                 for yy=1:size(Y,2)
%                     beta_y(:,yy) = regress(Y(:,yy),CVA.Covs);
%                 end
%                 Yr = Y - (CVA.Covs * beta_y);%                 

%                 beta_y = mvregress(CVA.Covs, Y);
%                 
end


%%

% if size(X,1) ~= size(Y,1)
%     error('X and Y must have the same number of observations')
% end
% if CCA.nComp > min([N,P,Q])
%     error('nComp cannot be larger than the dimensions of X or Y')
% end

if ischar(lambdaX) && strcmp(lambdaX,'auto')
    [~, lambdaX] = cov1para(X);
%     [~, opt.lambdaX] = stats_shrinkage_cov(X,'rblw');
end

if ischar(lambdaY) && strcmp(lambdaY,'auto')
    [~, lambdaY] = cov1para(Y);
%     [~, opt.lambdaY] = stats_shrinkage_cov(Y,'rblw');
end


%% Preprocessing
% Center the variables
Xorig= X;
Yorig= Y;

meanX= mean(X,1);
meanY= mean(Y,1);

X = X - repmat(meanX, N, 1);
Y = Y - repmat(meanY, N, 1);


%% Get (cross-)covariance matrices of X and Y
Cxx= (X'*X)/N; %cov(X);%
Cyy= (Y'*Y)/N; %cov(Y);

% Regularization matrices (scaled identity matrices)
Ip= trace(Cxx)/P * eye(P);
Iq= trace(Cyy)/Q * eye(Q);

% Cross-covariance matrices
Cxy= (X' * Y)/N;
Cyx= Cxy';

%% Run CCA 
% warning('Start looping over lambdas\n');
% warning('Iteration x=%d, y=%d [dim(X)=%d, dim(Y)=%d, lambdax=%0.6f, lambday=%0.6f]\n',ix,iy,P,Q,lambdaX,lambdaY);
Cxx_reg= (1-lambdaX) * Cxx + lambdaX * Ip;
Cyy_reg= (1-lambdaY) * Cyy + lambdaY * Iq;



%% Cholesky decomposition of Cyy to ensure symmetric EV problem
%% (see Hardoon et al p. 2644)
[Ryy,p]= chol(Cyy_reg);

%         iRyy= inv(Ryy);

% We perform a change in coordinates Uy = Ryy * Wy and thus obtain
% a symmetric eigenvalue problem  A*Uy = lambda * Uy
%         A= iRyy' * Cyx * (Cxx_reg \ Cxy) * iRyy;
A= (Ryy' \ Cyx) * (Cxx_reg \ Cxy) / Ryy;


% symmetric and positive definite  (can be non-symmetric due to numerical problems)
A= 0.5*(A + A');

% Ensure symmetric and positive matrix with eig
% https://uk.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
%         A = nearestSPD(A); It may enter in endless loop 
% Using https://github.com/higham/modified-cholesky, which seems to
% work similarly to nearestSPD, but does not enter in endless loop
[L, D, P, D0] = modchol_ldlt(A); 
A = P'*L*D*L'*P;      % Modified matrix: symmetric pos def.        


% try chol(A);
% %     disp('Matrix is symmetric positive definite.')
% catch
%     A = csa_util_nearestSPD(A);
% end


% 
% opts     = [];
% opts.tol = 1e-300; % 'Tolerance'
% opts.cholB = true; %'IsCholesky'
% opts.spdB  = 1;
% opts.issym = 1;
% opts.disp  = 1;
% if isempty(CCA.numComp) 
%     [Uy, D]= eig(double(A));  
% else
    [Uy, D]= eigs(double(A), opt.numComp,'largestabs');%,opts);%  eig returns ALL EVs and it returns them in no particular order
% end

% if any(any(D<0))
%     disp('N.B.!!! Negative value in matrix D of eigenvalues!!!!  Flipping the sign!!!' );
%     D(D<0)=D(D<0)*-1;
% end

% % Sort D and Uy in descending order, Reordering accorind to Eigenvalues, needed for eig
% [D,idx] = sort(diag(D),'descend');
% D = diag(D);
% Uy = Uy(:, idx); 

% Change coordinates back to Wy
YW = real(Ryy \ Uy);

%% Get CCA coefficients (following formula from the Vinod paper)
%%% Get Yw's as eigenvectors
%         MAT=  (Cyy_reg \ Cyx) * (Cxx_reg \ Cxy);
%         MAT= 0.5 * ( MAT + MAT');  % ensure that it's symmetric (can be slightly asymemtric due to neumerical issues I think)
%         [Uy,D]= eigs(MAT,opt.nComp);
%         
%         % Sort D and Uy in descending order
%         [D,idx] = sort(diag(D),'descend');
%         D = diag(D);
%         Uy = Uy(:, idx); 
%         Yw_eigs{ix,iy} =Uy;



% Add small values along the diagonal zeros
% diagZero = diag(diag(D==0));
% diagNonZero = diag(diag(D~=0));
% D(diagZero) = min([1e-30; D(diagNonZero)]);

%%% Xw's can be obtained from Yws as (Vinod et al paper)
XW =  real((Cxx_reg \ Cxy) * YW / sqrtm(D));

% -----------------------------
% Calculate canonical Variates
% -----------------------------
XS = X*XW;       % Canonical Variate X
YS = Y*YW;       % Canonical Variate Y

% -----------------------------
% Canonical Loadings
% -----------------------------
XL = corr(X,XS); % Canonical Loadings X | CCA.X' * CCA.XS; %Xorig' * Xs{ix,iy};
YL = corr(Y,YS); % Canonical Loadings Y | CCA.Y' * CCA.YS; %Yorig' * Ys{ix,iy};

% -------------------------------------
% Canonical Correlation and Covariance
% -------------------------------------
[rval,pval] = corr(XS,YS);
Cov = diag(XS' * YS)/N; %diag(cov([CCA.XS CCA.YS]),opt.numComp);% Calculate covariance between subject scores


%% Assemble output
CCA   = opt;
CCA.X  = X; % Mean centered X data
CCA.Y  = Y; % Mean centered Y data
CCA.XW = XW;
CCA.YW = YW;
CCA.XS = XS; % Canonical Variate X
CCA.YS = YS; % Canonical Variate Y
CCA.XL = XL; % Loadings X
CCA.YL = YL; % Loadings Y

% ----------------------
% Canonical correlations
% ----------------------
CCA.R  = diag(rval);
CCA.P  = diag(pval);

% ----------------------
% Canonical Covariance
% ----------------------
CCA.Cov = Cov;


% % If correlation are negative, flip the respective X variables to make them
% % positive
% for ii=1:numel(CCA.R)
%     if CCA.R(ii)<0
%         %         warning('Canonical pair #%d yields a negative correlation, flipping the correlation and the respective X components',ii)
%         CCA.R(ii)    = -CCA.R(ii);
%         CCA.XS(:,ii) = -CCA.XS(:,ii);
%         CCA.XW(:,ii) = -CCA.XW(:,ii);
%     end
% end




%% Apply to new validation data?

if isfield(opt,'Xte')
    
%     CCA.XSte= cell(nLambdaX,nLambdaY);
%     CCA.YSte= cell(nLambdaX,nLambdaY);
%     CCA.XLte= cell(nLambdaX,nLambdaY);
%     CCA.YLte= cell(nLambdaX,nLambdaY);
%     CCA.Rte= cell(nLambdaX,nLambdaY);
%     CCA.Pte= cell(nLambdaX,nLambdaY);
    
    Xte= opt.Xte;
    Yte= opt.Yte;
    
    % Center the validation data (using means on the training data)
    Nval= size(Xte,1);
    Xte = Xte - repmat(meanX, Nval, 1);
    Yte = Yte - repmat(meanY, Nval, 1);
    
    CCA.XSte   = Xte * CCA.XW;
    CCA.YSte   = Yte * CCA.YW;
    [Rte,Pte]  = corr(CCA.XSte,CCA.YSte); % The traanspose perhpas?
    CCA.Rte    = diag(Rte);
    CCA.Pte    = diag(Pte);
    CCA.Cte    = diag(cov([CCA.XSte,CCA.YSte]),CCA.numComp);
    
    CCA.XWte = Xte' * Yte * CCA.YW;% As in Kovacevic et al, See csa_stats_rCVA_bootstr for more information
    CCA.YWte = Yte' * Xte * CCA.XW;% 
    
    % Canonical loadings (non-normalised)
%     CCA.XLte   = Xte' * CCA.XSte;
%     CCA.YLte   = Yte' * CCA.YSte;
    
    CCA.XLte   = corr(Xte,CCA.XSte);
    CCA.YLte   = corr(Yte,CCA.YSte);
    
    
end
