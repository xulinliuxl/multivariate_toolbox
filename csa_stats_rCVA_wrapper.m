function [CVA,CCApart,CCAperm] = csa_stat_rCVA_wrapper(CCA)
%% ---------------------------------------------------------------------
%  Reguralized CCA for full- and partial-FC predicting performance using
%  DW - data-wise permutations
%  FW - fold-wise permutations
%  PW - partition-wise permutations
%  SH - half-split reliability (Kovacevic et al 2013)
% ----------------------------------------------------------------------

CVA        = csa_stats_rCVA_params();

% ---------------------------------------------------
% Update default settings with user specific settings
% ---------------------------------------------------
fieldsL1User    = fieldnames(CCA);
fieldsL1Default = fieldnames(CVA);
 
%% Go through all fields on first level
for iL1 = 1:numel(fieldsL1User)
    fnL1U = fieldsL1User{iL1};

    % Check whether fieldname matches default fieldnames
    if sum(ismember(fieldsL1Default,fnL1U)) == 0
        fprintf('Field name: %s \n',fnL1U);
        error('Non-existent field has been specified. Please check params function for detials. \n');
    end

    if isstruct(CCA.(fnL1U)) % Second-level structure
        fieldsL2User = fieldnames(CCA.(fnL1U));
        fieldsL2Default = fieldnames(CVA.(fnL1U));
        
        for iL2 = 1:numel(fieldsL2User)
            fnL2U = fieldsL2User{iL2};
        
            % Check whether fieldname matches default fieldnames
            if sum(ismember(fieldsL2Default,fnL2U)) == 0
                fprintf('Field name: %s > %s \n',fnL1U,fnL2U);
                error('Non-existent field has been specified. Please check params function for detials. \n');
            end
            
            if isstruct(CCA.(fnL1U).(fnL2U)) % Third-level structure
                fieldsL3User = fieldnames(CCA.(fnL1U).(fnL2U));
                fieldsL3Default = fieldnames(CVA.(fnL1U).(fnL2U));
                
                for iL3 = 1:numel(fieldsL3User)
                    fnL3U = fieldsL3User{iL3};
                    
                    % Check whether fieldname matches default fieldnames
                    if sum(ismember(fieldsL3Default,fnL3U)) == 0
                        fprintf('Field name: %s > %s > %s \n',fnL1U,fnL2U,fnL3U);
                        error('Non-existent field has been specified. Please check params function for detials. \n');
                    end
                    CVA.(fnL1U).(fnL2U).(fnL3U) = CCA.(fnL1U).(fnL2U).(fnL3U);
                end
                
            else
                CVA.(fnL1U).(fnL2U) = CCA.(fnL1U).(fnL2U);
            end
        end
    else
        CVA.(fnL1U) = CCA.(fnL1U);
    end
end
    
    %%

[CVA.Ns,CVA.Ng] = size(CVA.Y);
CVA.numVarX  = size(CVA.X,2);
CVA.numVarY  = size(CVA.Y,2);

if isempty(CVA.Covs)
    CVA.Covs = ones(CVA.Ns,1);
end

% Set default number of components if it not prespecified 
if isempty(CVA.numComp)
    CVA.numComp  = min([CVA.Ns,CVA.numVarX,CVA.numVarY]);
end


% Predefine dummy output structures, which are not present for some CCA
% modes
CCApart = [];
CCAperm = [];

% rng(CVA.partSeed,'combRecursive');

%% Determine in what mode to run CCA
% cross-validated CCA
% Standard permuations
% Split-half permuations
% Classical CCA (neither cross-validations, nor permuations)

namemode = string(fieldnames(CVA.mode));
for m = 1:numel(namemode)
    doMode(m) = CVA.mode.(namemode{m}).do;
end
doMode = logical(doMode);
% Check that only one mode has been requested
if sum(doMode)>1
    error('More than one CCA modes has been requested at the same time!. Select only one mode');
elseif sum(doMode) == 0
    doMode(4) = 1;
    CVA.mode.(namemode(doMode)).do = 1;
end
% doMode  = find(doMode);
ccaMode = CVA.mode.(namemode(doMode)).name;

fprintf('\nRunning CCA in mode: %s. \n \n',ccaMode);
switch ccaMode
    case  'Cross-Validation'
        [CVA,CCApart,CCAperm] = csa_stats_rCVA_CV(CVA);
        
    case  'Bootstrap Permutations'
        CVA = csa_stats_rCVA_permBootstr(CVA);

    case 'Classic permutations'
        CVA = csa_stats_rCVA_permClassic(CVA);

    case 'Standard CCA'
        CVA = csa_proc_CCA(CVA,CVA.X,CVA.Y);

end

CVA.ccaMode = ccaMode;
