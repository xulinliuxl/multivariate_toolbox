function CCA = csa_stats_rCVA_params()
% Default parameters for rCVA

CCA.X               = [];
CCA.Y               = [];
CCA.cfg             = [];
CCA.doCovs          = 0;
CCA.Covs            = [];
CCA.lambdaX         = 0;% 'auto'
CCA.lambdaY         = 0;% 'auto'
CCA.runInSerial     = 0;
CCA.numComp         = [];
CCA.doAdjustDOF     = 0;
CCA.nameAnalysis    = '';
CCA.dirOut          = '';
CCA.dirOutNulls     = '';
CCA.doSaveNulls     = 0;
CCA.usePresetRandOrder = 0;
CCA.randOrder       = [];



%% Cross-validation settings
CCA.mode.cv.do        = 0;
CCA.mode.cv.name      = 'Cross-Validation';
% CCA.mode.cv.nameShort = 'cv';
CCA.mode.cv.CVtype    = 'kFold';% 'splitHalf' | 'kFold' | 'LeaveOut'
CCA.mode.cv.numFolds  = 5;
CCA.mode.cv.labels    = [];     % Labelling samples to ensure partitions with equal numbersdrawn from each group/decile     
CCA.mode.cv.numPart   = 1;      % Number of repatitionings
CCA.mode.cv.partSeed  = 1;
CCA.mode.cv.doCovFW   = 1;      % Regress out covariates in each fold

CCA.mode.cv.permutePW = 0;      % Permute partition-wise - it permutes the entire data
                                % before cross-validations, i.e. crossvalidations on permuted data
CCA.mode.cv.numPermPW = 5;

CCA.mode.cv.permuteFW = 0;      % Permute fold-wise, i.e. each fold is permuted (takes time to run!)
CCA.mode.cv.numPermFW = 1;

% CCA.mode.cv.permuteDW = 0;      % Permute data-wise, i.e. perumatation on the entire sampe
% CCA.mode.cv.numPermDW = 10;

CCA.mode.cv.doSplitHalfFW = 0;  % CCA on entire sample then half-split; 
                        % Split-half stability (SS) approach, which measures how consistently PLS 
                        % reconstructs U on data split=halves, using V estimated on the full 
                        % dataset (and vice-versa), see Kovacevic et al 2013
CCA.mode.cv.doSplitHalfDW = 0;  % Half-split sample, then CCA within each split
                        % Split-half reproducibility (SR) approach, which compares U and V
                        % estimated on independent data split-halves, see Churchill et al 2013 
CCA.mode.cv.numBootDW = 1000; 
CCA.mode.cv.doSaveNullsFW = 0;
CCA.mode.cv.doSaveNullsDW = 0;

%% Classical permuations
CCA.mode.permClassic.do      = 0;
CCA.mode.permClassic.name    = 'Classic permutations';
CCA.mode.permClassic.numPerm = 500;

%% Half-Split Permutations
CCA.mode.permBootstr.do       = 0;
CCA.mode.permBootstr.name     = 'Bootstrap Permutations'; 
CCA.mode.permBootstr.numPerm  = 500;
CCA.mode.permBootstr.labels   = [];

%% Standard CCA (No CV or permutaions)
CCA.mode.standard.do        = 0;
CCA.mode.standard.name      = 'Standard CCA';