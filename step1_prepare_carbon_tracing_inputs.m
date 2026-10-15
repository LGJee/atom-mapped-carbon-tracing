%% Prepare carbon-tracing inputs with six explicit COBRA model blocks
% Edit paths/analysis settings here and the COBRA calls in each named model block.
% Each block uses modelWT/modelIUPD as in COBRA tutorials and saves a named result.
% Checks, medium filtering and export details are in the helper section below.
clear; clc;
projectFolder = fileparts(mfilename('fullpath'));

%% Common paths and numerical settings
cfg.resultsDirectory = 'results';
cfg.cobraFolder = ''; % Optional; normally initialize COBRA before running.
cfg.solver = 'glpk';
cfg.saccharomycesModelFile = 'data/yeast9/yeast-GEM.mat';
cfg.escherichiaModelFile = 'data/ecoli/iML1515.mat';
cfg.escherichiaAnnotationFile = 'data/ecoli/iML1515.json'; % Optional identifiers.
cfg.yarrowiaModelFile = 'data/yarrowia/iYali5.mat';
cfg.yarrowiaGlucoseUptake = 1; % iYali5 ships with glucose closed; glucose-medium limit.
cfg.nonzeroTolerance = 1e-8;

%% Ordered carbon sources for the six model blocks
cfg.carbonSources = struct();
cfg.carbonSources.Saccharomyces_cerevisiae_WT = {'D-glucose[e]'};
cfg.carbonSources.Saccharomyces_cerevisiae_IUPD = {'D-glucose[e]','isoprenol[e]'};
cfg.carbonSources.Escherichia_coli_WT = {'D-glucose[e]'};
cfg.carbonSources.Escherichia_coli_IUPD = {'D-glucose[e]','isoprenol[e]'};
cfg.carbonSources.Yarrowia_lipolytica_WT = {'D-glucose[e]'};
cfg.carbonSources.Yarrowia_lipolytica_IUPD = {'D-glucose[e]','isoprenol[e]'};

%% Initialize COBRA and prepare output bookkeeping
if ~isempty(cfg.cobraFolder)
    addpath(cfg.cobraFolder);
    initCobraToolbox(false);
end
% Select the LP solver.
changeCobraSolver(cfg.solver,'LP');
[cfg,plan,planFile,names]=initializeAnalysis(cfg,projectFolder);

%% 1. Saccharomyces_cerevisiae_WT
% Load the original model independently; no already-constrained WT is reused.
data = load(fullfile(projectFolder,cfg.saccharomycesModelFile));
modelWT = data.model;
checkModelInput(modelWT,'Saccharomyces_cerevisiae_WT');
% Read the native glucose bound; keep both bounds unchanged for growth optimization.
glucoseLowerBoundWT = modelWT.lb(strcmp(modelWT.rxns,'r_1714')); % D-glucose exchange.
assert(isscalar(glucoseLowerBoundWT) && isfinite(glucoseLowerBoundWT) && glucoseLowerBoundWT<0, ...
    'The source model must allow glucose uptake.');

% Culture conditions. Negative exchange flux denotes uptake.
% Inherit native oxygen and remaining Yeast9 medium bounds.

% Maximize growth, then fix 50% of this strain's own maximum growth rate.
% The loaded model already uses the checked biomass objective.
growthWT = optimizeCbModel(modelWT,'max');
checkSolution(modelWT,growthWT,'Saccharomyces_cerevisiae_WT');
maxGrowthSaccharomycesWT = growthWT.f; % Growth reference shared by WT and IUPD.
modelWT = changeRxnBounds(modelWT,'r_2111',0.5*growthWT.f,'b'); % Growth (biomass objective).
% Double the native glucose uptake limit only after fixing growth.
modelWT = changeRxnBounds(modelWT,'r_1714',2*glucoseLowerBoundWT,'l'); % D-glucose exchange.

% Maximize squalene, then minimize total absolute flux ('one' = L1 pFBA).
modelWT = changeObjective(modelWT,'r_1012'); % Squalene synthase.
solutionWT = optimizeCbModel(modelWT,'max','one');
checkSolution(modelWT,solutionWT,'Saccharomyces_cerevisiae_WT');

% Keep this model and solution before the next block reuses the working variables.
modelSaccharomycesWT = modelWT;
solutionSaccharomycesWT = solutionWT;

%% 2. Saccharomyces_cerevisiae_IUPD
% Load the original model independently; no already-constrained WT is reused.
data = load(fullfile(projectFolder,cfg.saccharomycesModelFile));
modelIUPD = data.model;
checkModelInput(modelIUPD,'Saccharomyces_cerevisiae_IUPD');
% Keep the source model's native glucose bounds throughout the IUPD analysis.
glucoseLowerBoundIUPD = modelIUPD.lb(strcmp(modelIUPD.rxns,'r_1714')); % D-glucose exchange.
assert(isscalar(glucoseLowerBoundIUPD) && isfinite(glucoseLowerBoundIUPD) && glucoseLowerBoundIUPD<0, ...
    'The source model must allow glucose uptake.');

% Culture conditions. Negative exchange flux denotes uptake.
% Inherit native oxygen and remaining Yeast9 medium bounds.
modelIUPD = generateRules(modelIUPD);

% Knock out ERG13 and introduce the two-step isoprenol utilization pathway.
modelBeforeKnockout = modelIUPD;
modelIUPD = deleteModelGenes(modelIUPD,{'YML126C'});
checkKnockout(modelBeforeKnockout,modelIUPD);
modelIUPD = addMetabolite(modelIUPD,'isoprenol[e]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isoprenol[c]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isopentenyl_phosphate[c]', ...
    'metName','isopentenyl phosphate','metFormula','C5H11O4P');
% Ex_isopentenol_e: isoprenol uptake limit = 2.09 times the native glucose limit.
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e','reactionFormula', ...
    'isoprenol[e] -> ','lowerBound',2.09*glucoseLowerBoundIUPD,'upperBound',1000);
% Ex_isopentenol_e_t: isoprenol transport from extracellular space to cytosol.
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e_t','reactionFormula', ...
    'isoprenol[e] -> isoprenol[c]','lowerBound',0,'upperBound',1000);
% isopentenolK: phosphorylate isoprenol to isopentenyl phosphate (ATP -> ADP).
% s_0434: ATP; s_0394: ADP; s_0943: isopentenyl diphosphate (IPP), all cytosolic.
modelIUPD = addReaction(modelIUPD,'isopentenolK','reactionFormula', ...
    'isoprenol[c] + s_0434 -> isopentenyl_phosphate[c] + s_0394', ...
    'geneRule','YLR133W','lowerBound',0,'upperBound',1000);
% PisopentenolK: phosphorylate isopentenyl phosphate to IPP (ATP -> ADP).
modelIUPD = addReaction(modelIUPD,'PisopentenolK','reactionFormula', ...
    'isopentenyl_phosphate[c] + s_0434 -> s_0943 + s_0394', ...
    'geneRule','AtIPK','lowerBound',0,'upperBound',1000);

% Fix the same growth rate as WT: 50% of the matching WT maximum.
assert(exist('maxGrowthSaccharomycesWT','var')==1,'Run the matching WT model block first.');
modelIUPD = changeRxnBounds(modelIUPD,'r_2111',0.5*maxGrowthSaccharomycesWT,'b'); % Growth (biomass objective).

% Maximize squalene, then minimize total absolute flux ('one' = L1 pFBA).
modelIUPD = changeObjective(modelIUPD,'r_1012'); % Squalene synthase.
solutionIUPD = optimizeCbModel(modelIUPD,'max','one');
checkSolution(modelIUPD,solutionIUPD,'Saccharomyces_cerevisiae_IUPD');

% Keep this model and solution before the next block reuses the working variables.
modelSaccharomycesIUPD = modelIUPD;
solutionSaccharomycesIUPD = solutionIUPD;

%% 3. Escherichia_coli_WT
% Load the original model independently; no already-constrained WT is reused.
data = load(fullfile(projectFolder,cfg.escherichiaModelFile));
modelWT = data.iML1515;
checkModelInput(modelWT,'Escherichia_coli_WT');
% Read the native glucose bound; keep both bounds unchanged for growth optimization.
glucoseLowerBoundWT = modelWT.lb(strcmp(modelWT.rxns,'EX_glc__D_e')); % D-glucose exchange.
assert(isscalar(glucoseLowerBoundWT) && isfinite(glucoseLowerBoundWT) && glucoseLowerBoundWT<0, ...
    'The source model must allow glucose uptake.');

% Both E. coli comparators need SQS and an output for accumulated squalene.
% RHEA:32295; this is a functional addition, not a specified enzyme construct.
modelWT = addMetabolite(modelWT,'squalene[c]', ...
    'metName','squalene','metFormula','C30H50','Charge',0);
% SQS: introduced squalene synthase; converts two FPP molecules to squalene.
modelWT = addReaction(modelWT,'SQS', ...
    'reactionName','Heterologous squalene synthase','reactionFormula', ...
    '2 frdp_c + nadph_c + h_c -> squalene[c] + 2 ppi_c + nadp_c', ...
    'lowerBound',0,'upperBound',1000);
% DM_squalene: squalene demand reaction (output sink).
modelWT = addReaction(modelWT,'DM_squalene', ...
    'reactionName','Squalene output','reactionFormula','squalene[c] -> ', ...
    'lowerBound',0,'upperBound',1000);

% Culture conditions. Negative exchange flux denotes uptake.
% Close other carbon uptake; inherit native oxygen, inorganic medium and maintenance.
modelWT = closeOtherCarbonUptake(modelWT,'EX_glc__D_e'); % D-glucose exchange.

% Maximize growth, then fix 50% of this strain's own maximum growth rate.
% The loaded model already uses the checked biomass objective.
growthWT = optimizeCbModel(modelWT,'max');
checkSolution(modelWT,growthWT,'Escherichia_coli_WT');
maxGrowthEscherichiaWT = growthWT.f; % Growth reference shared by WT and IUPD.
modelWT = changeRxnBounds(modelWT,'BIOMASS_Ec_iML1515_core_75p37M',0.5*growthWT.f,'b'); % Core biomass objective (iML1515).
% Double the native glucose uptake limit only after fixing growth.
modelWT = changeRxnBounds(modelWT,'EX_glc__D_e',2*glucoseLowerBoundWT,'l'); % D-glucose exchange.

% Maximize squalene, then minimize total absolute flux ('one' = L1 pFBA).
modelWT = changeObjective(modelWT,'SQS'); % Introduced squalene synthase.
solutionWT = optimizeCbModel(modelWT,'max','one');
checkSolution(modelWT,solutionWT,'Escherichia_coli_WT');

% Keep this model and solution before the next block reuses the working variables.
modelEscherichiaWT = modelWT;
solutionEscherichiaWT = solutionWT;

%% 4. Escherichia_coli_IUPD
% Load the original model independently; no already-constrained WT is reused.
data = load(fullfile(projectFolder,cfg.escherichiaModelFile));
modelIUPD = data.iML1515;
checkModelInput(modelIUPD,'Escherichia_coli_IUPD');
% Keep the source model's native glucose bounds throughout the IUPD analysis.
glucoseLowerBoundIUPD = modelIUPD.lb(strcmp(modelIUPD.rxns,'EX_glc__D_e')); % D-glucose exchange.
assert(isscalar(glucoseLowerBoundIUPD) && isfinite(glucoseLowerBoundIUPD) && glucoseLowerBoundIUPD<0, ...
    'The source model must allow glucose uptake.');

% Both E. coli comparators need SQS and an output for accumulated squalene.
% RHEA:32295; this is a functional addition, not a specified enzyme construct.
modelIUPD = addMetabolite(modelIUPD,'squalene[c]', ...
    'metName','squalene','metFormula','C30H50','Charge',0);
% SQS: introduced squalene synthase; converts two FPP molecules to squalene.
modelIUPD = addReaction(modelIUPD,'SQS', ...
    'reactionName','Heterologous squalene synthase','reactionFormula', ...
    '2 frdp_c + nadph_c + h_c -> squalene[c] + 2 ppi_c + nadp_c', ...
    'lowerBound',0,'upperBound',1000);
% DM_squalene: squalene demand reaction (output sink).
modelIUPD = addReaction(modelIUPD,'DM_squalene', ...
    'reactionName','Squalene output','reactionFormula','squalene[c] -> ', ...
    'lowerBound',0,'upperBound',1000);

% Culture conditions. Negative exchange flux denotes uptake.
% Close other carbon uptake; inherit native oxygen, inorganic medium and maintenance.
modelIUPD = closeOtherCarbonUptake(modelIUPD,'EX_glc__D_e'); % D-glucose exchange.
modelIUPD = generateRules(modelIUPD);

% Knock out ispG and introduce the two-step isoprenol utilization pathway.
modelBeforeKnockout = modelIUPD;
modelIUPD = deleteModelGenes(modelIUPD,{'b2515'});
checkKnockout(modelBeforeKnockout,modelIUPD);
modelIUPD = addMetabolite(modelIUPD,'isoprenol[e]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isoprenol[c]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isopentenyl_phosphate[c]', ...
    'metName','isopentenyl phosphate','metFormula','C5H11O4P');
% Ex_isopentenol_e: isoprenol uptake limit = 2.09 times the native glucose limit.
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e','reactionFormula', ...
    'isoprenol[e] -> ','lowerBound',2.09*glucoseLowerBoundIUPD,'upperBound',1000);
% Ex_isopentenol_e_t: isoprenol transport from extracellular space to cytosol.
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e_t','reactionFormula', ...
    'isoprenol[e] -> isoprenol[c]','lowerBound',0,'upperBound',1000);
% isopentenolK: phosphorylate isoprenol to isopentenyl phosphate (ATP -> ADP).
% atp_c: ATP; adp_c: ADP; ipdp_c: isopentenyl diphosphate (IPP), all cytosolic.
modelIUPD = addReaction(modelIUPD,'isopentenolK','reactionFormula', ...
    'isoprenol[c] + atp_c -> isopentenyl_phosphate[c] + adp_c', ...
    'geneRule','YLR133W','lowerBound',0,'upperBound',1000);
% PisopentenolK: phosphorylate isopentenyl phosphate to IPP (ATP -> ADP).
modelIUPD = addReaction(modelIUPD,'PisopentenolK','reactionFormula', ...
    'isopentenyl_phosphate[c] + atp_c -> ipdp_c + adp_c', ...
    'geneRule','AtIPK','lowerBound',0,'upperBound',1000);

% Fix the same growth rate as WT: 50% of the matching WT maximum.
assert(exist('maxGrowthEscherichiaWT','var')==1,'Run the matching WT model block first.');
modelIUPD = changeRxnBounds(modelIUPD,'BIOMASS_Ec_iML1515_core_75p37M',0.5*maxGrowthEscherichiaWT,'b'); % Core biomass objective (iML1515).

% Maximize squalene, then minimize total absolute flux ('one' = L1 pFBA).
modelIUPD = changeObjective(modelIUPD,'SQS'); % Introduced squalene synthase.
solutionIUPD = optimizeCbModel(modelIUPD,'max','one');
checkSolution(modelIUPD,solutionIUPD,'Escherichia_coli_IUPD');

% Keep this model and solution before the next block reuses the working variables.
modelEscherichiaIUPD = modelIUPD;
solutionEscherichiaIUPD = solutionIUPD;

%% 5. Yarrowia_lipolytica_WT
% Load the original model independently; no already-constrained WT is reused.
data = load(fullfile(projectFolder,cfg.yarrowiaModelFile));
modelWT = data.model;
checkModelInput(modelWT,'Yarrowia_lipolytica_WT');
% The published iYali5 medium uses glycerol; explicitly enable glucose instead.
assert(isscalar(cfg.yarrowiaGlucoseUptake) && isfinite(cfg.yarrowiaGlucoseUptake) && ...
    cfg.yarrowiaGlucoseUptake>0,'Set a positive Yarrowia glucose uptake limit.');
modelWT = changeRxnBounds(modelWT,'1714',-cfg.yarrowiaGlucoseUptake,'l'); % D-glucose exchange.
modelWT = changeRxnBounds(modelWT,'1714',1000,'u'); % Permit secretion; do not force uptake.
% Read the glucose-medium bound before growth optimization.
glucoseLowerBoundWT = modelWT.lb(strcmp(modelWT.rxns,'1714')); % D-glucose exchange.
assert(isscalar(glucoseLowerBoundWT) && isfinite(glucoseLowerBoundWT) && glucoseLowerBoundWT<0, ...
    'The source model must allow glucose uptake.');

% Culture conditions. Negative exchange flux denotes uptake.
% Close other carbon uptake; inherit native oxygen, inorganic medium and maintenance.
modelWT = closeOtherCarbonUptake(modelWT,'1714'); % D-glucose exchange.

% Maximize growth, then fix 50% of this strain's own maximum growth rate.
% The loaded model already uses the checked biomass objective.
growthWT = optimizeCbModel(modelWT,'max');
checkSolution(modelWT,growthWT,'Yarrowia_lipolytica_WT');
maxGrowthYarrowiaWT = growthWT.f; % Growth reference shared by WT and IUPD.
modelWT = changeRxnBounds(modelWT,'xBIOMASS',0.5*growthWT.f,'b'); % Biomass production objective.
% Double the native glucose uptake limit only after fixing growth.
modelWT = changeRxnBounds(modelWT,'1714',2*glucoseLowerBoundWT,'l'); % D-glucose exchange.

% Maximize squalene, then minimize total absolute flux ('one' = L1 pFBA).
modelWT = changeObjective(modelWT,'1012'); % Squalene synthase.
solutionWT = optimizeCbModel(modelWT,'max','one');
checkSolution(modelWT,solutionWT,'Yarrowia_lipolytica_WT');

% Keep this model and solution before the next block reuses the working variables.
modelYarrowiaWT = modelWT;
solutionYarrowiaWT = solutionWT;

%% 6. Yarrowia_lipolytica_IUPD
% Load the original model independently; no already-constrained WT is reused.
data = load(fullfile(projectFolder,cfg.yarrowiaModelFile));
modelIUPD = data.model;
checkModelInput(modelIUPD,'Yarrowia_lipolytica_IUPD');
% The published iYali5 medium uses glycerol; explicitly enable glucose instead.
assert(isscalar(cfg.yarrowiaGlucoseUptake) && isfinite(cfg.yarrowiaGlucoseUptake) && ...
    cfg.yarrowiaGlucoseUptake>0,'Set a positive Yarrowia glucose uptake limit.');
modelIUPD = changeRxnBounds(modelIUPD,'1714',-cfg.yarrowiaGlucoseUptake,'l'); % D-glucose exchange.
modelIUPD = changeRxnBounds(modelIUPD,'1714',1000,'u'); % Permit secretion; do not force uptake.
% Keep the configured glucose-medium bounds throughout the IUPD analysis.
glucoseLowerBoundIUPD = modelIUPD.lb(strcmp(modelIUPD.rxns,'1714')); % D-glucose exchange.
assert(isscalar(glucoseLowerBoundIUPD) && isfinite(glucoseLowerBoundIUPD) && glucoseLowerBoundIUPD<0, ...
    'The source model must allow glucose uptake.');

% Culture conditions. Negative exchange flux denotes uptake.
% Close other carbon uptake; inherit native oxygen, inorganic medium and maintenance.
modelIUPD = closeOtherCarbonUptake(modelIUPD,'1714'); % D-glucose exchange.
modelIUPD = generateRules(modelIUPD);

% Knock out ERG13 and introduce the two-step isoprenol utilization pathway.
modelBeforeKnockout = modelIUPD;
modelIUPD = deleteModelGenes(modelIUPD,{'YALI0F30481g'});
checkKnockout(modelBeforeKnockout,modelIUPD);
modelIUPD = addMetabolite(modelIUPD,'isoprenol[e]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isoprenol[c]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isopentenyl_phosphate[c]', ...
    'metName','isopentenyl phosphate','metFormula','C5H11O4P');
% Ex_isopentenol_e: isoprenol uptake limit = 2.09 times the native glucose limit.
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e','reactionFormula', ...
    'isoprenol[e] -> ','lowerBound',2.09*glucoseLowerBoundIUPD,'upperBound',1000);
% Ex_isopentenol_e_t: isoprenol transport from extracellular space to cytosol.
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e_t','reactionFormula', ...
    'isoprenol[e] -> isoprenol[c]','lowerBound',0,'upperBound',1000);
% isopentenolK: phosphorylate isoprenol to isopentenyl phosphate (ATP -> ADP).
% m141: ATP; m143: ADP; m491: isopentenyl diphosphate (IPP), all cytosolic.
modelIUPD = addReaction(modelIUPD,'isopentenolK','reactionFormula', ...
    'isoprenol[c] + m141 -> isopentenyl_phosphate[c] + m143', ...
    'geneRule','YLR133W','lowerBound',0,'upperBound',1000);
% PisopentenolK: phosphorylate isopentenyl phosphate to IPP (ATP -> ADP).
modelIUPD = addReaction(modelIUPD,'PisopentenolK','reactionFormula', ...
    'isopentenyl_phosphate[c] + m141 -> m491 + m143', ...
    'geneRule','AtIPK','lowerBound',0,'upperBound',1000);

% Fix the same growth rate as WT: 50% of the matching WT maximum.
assert(exist('maxGrowthYarrowiaWT','var')==1,'Run the matching WT model block first.');
modelIUPD = changeRxnBounds(modelIUPD,'xBIOMASS',0.5*maxGrowthYarrowiaWT,'b'); % Biomass production objective.

% Maximize squalene, then minimize total absolute flux ('one' = L1 pFBA).
modelIUPD = changeObjective(modelIUPD,'1012'); % Squalene synthase.
solutionIUPD = optimizeCbModel(modelIUPD,'max','one');
checkSolution(modelIUPD,solutionIUPD,'Yarrowia_lipolytica_IUPD');

% Keep this model and solution before the next block reuses the working variables.
modelYarrowiaIUPD = modelIUPD;
solutionYarrowiaIUPD = solutionIUPD;

%% Export the six inputs and generate the Python analysis plan
plan.models.Saccharomyces_cerevisiae_WT = exportModel( ...
    modelSaccharomycesWT,solutionSaccharomycesWT,'Saccharomyces_cerevisiae_WT', ...
    cfg.carbonSources.Saccharomyces_cerevisiae_WT,'',projectFolder,cfg);

plan.models.Saccharomyces_cerevisiae_IUPD = exportModel( ...
    modelSaccharomycesIUPD,solutionSaccharomycesIUPD,'Saccharomyces_cerevisiae_IUPD', ...
    cfg.carbonSources.Saccharomyces_cerevisiae_IUPD,'',projectFolder,cfg);

plan.models.Escherichia_coli_WT = exportModel( ...
    modelEscherichiaWT,solutionEscherichiaWT,'Escherichia_coli_WT', ...
    cfg.carbonSources.Escherichia_coli_WT,cfg.escherichiaAnnotationFile,projectFolder,cfg);

plan.models.Escherichia_coli_IUPD = exportModel( ...
    modelEscherichiaIUPD,solutionEscherichiaIUPD,'Escherichia_coli_IUPD', ...
    cfg.carbonSources.Escherichia_coli_IUPD,cfg.escherichiaAnnotationFile,projectFolder,cfg);

plan.models.Yarrowia_lipolytica_WT = exportModel( ...
    modelYarrowiaWT,solutionYarrowiaWT,'Yarrowia_lipolytica_WT', ...
    cfg.carbonSources.Yarrowia_lipolytica_WT,'',projectFolder,cfg);

plan.models.Yarrowia_lipolytica_IUPD = exportModel( ...
    modelYarrowiaIUPD,solutionYarrowiaIUPD,'Yarrowia_lipolytica_IUPD', ...
    cfg.carbonSources.Yarrowia_lipolytica_IUPD,'',projectFolder,cfg);

% Every carbon-source entry must have a matching model export.
assert(isequal(sort(fieldnames(plan.models)),sort(names(:))), ...
    'Each carbon-source entry needs a matching model block and export.');

%% Save portable paths only after every input has been exported successfully
temporaryFile = [planFile '.tmp'];
fid=fopen(temporaryFile,'w','n','UTF-8');
assert(fid~=-1,'Cannot write the analysis plan.');
try
    fprintf(fid,'%s\n',jsonencode(plan));
    fclose(fid);
catch exception
    fclose(fid); rethrow(exception);
end
movefile(temporaryFile,planFile,'f');
fprintf('Analysis plan saved: %s\n',planFile);
fprintf('In Linux/WSL, run:\npython3 step2_atom_mapping.py\npython3 step3_carbon_flow.py\n');

%% Shared Excel export and provenance helpers
function entry=exportModel(model,solution,name,sources,annotationFile,projectFolder,cfg)
    relativeFile=strrep(fullfile(cfg.resultsDirectory,[name '_1_pFBA_nonzero_reactions.xlsx']),'\','/');
    filename=fullfile(projectFolder,relativeFile);
    tolerance=cfg.nonzeroTolerance;
    rows=find(abs(solution.x)>tolerance);
    participating=find(any(model.S(:,rows)~=0,2));
    labels=metLabels(model);
    named=model; named.metNames=cellstr(labels);
    equation=string(printRxnFormula(named,'rxnAbbrList',model.rxns(rows),'printFlag',false,'metNameFlag',true));
    equation=replace(equation," -> "," => ");
    reactions=table(string(model.rxns(rows)),string(model.rxnNames(rows)),solution.x(rows),equation, ...
        'VariableNames',{'ReactionID','ReactionName','Flux','Equation'});
    metabolites=table(string(model.mets(participating)),labels(participating), ...
        'VariableNames',{'MetaboliteID','MetaboliteName'});
    metabolites.MolecularFormula=string(model.metFormulas(participating));
    if isfield(model,'metCharge'); metabolites.Charge=model.metCharge(participating);
    elseif isfield(model,'metCharges'); metabolites.Charge=model.metCharges(participating); end
    if ~isempty(annotationFile)
        annotations=jsondecode(fileread(fullfile(projectFolder,annotationFile)));
        keys=strings(numel(model.mets),1);
        for entry=annotations.metabolites(:)'
            position=find(strcmp(model.mets,entry.id));
            if ~isempty(position) && isfield(entry.annotation,'inchi_key')
                key=string(entry.annotation.inchi_key);
                if numel(key)==1; keys(position)=key; end
            end
        end
        metabolites.InChIKey=keys(participating);
    elseif isfield(model,'inchis')
        metabolites.InChI=string(model.inchis(participating));
    end
    % Write the whole workbook before replacing the previous input.
    temporaryFile=[filename(1:end-5) '.tmp.xlsx'];
    if isfile(temporaryFile); delete(temporaryFile); end
    writetable(reactions,temporaryFile,'Sheet','Nonzero');
    writetable(metabolites,temporaryFile,'Sheet','Metabolites');
    movefile(temporaryFile,filename,'f');
    fprintf('%s: objective %.9g; %d reactions, %d metabolites -> %s\n', ...
        name,solution.f,height(reactions),height(metabolites),filename);
    entry=struct('input_file',relativeFile,'carbon_sources',{sources},'input_sha256',fileSHA256(filename));
end

function labels=metLabels(model)
    labels=strings(numel(model.mets),1);
    for k=1:numel(model.mets)
        token=regexp(model.mets{k},'\[([^\]]+)\]$','tokens','once');
        if ~isempty(token); compartment=string(token{1});
        elseif isfield(model,'metComps') && isfield(model,'comps') && model.metComps(k)>=1 && model.metComps(k)<=numel(model.comps)
            compartment=string(model.comps{model.metComps(k)});
        else
            token=regexp(model.mets{k},'_([^_]+)$','tokens','once');
            assert(~isempty(token),'No identifiable compartment for %s.',model.mets{k});
            compartment=string(token{1});
        end
        name=string(model.metNames{k});
        if name=="D-Glucose"; name="D-glucose"; end
        labels(k)=name+"["+compartment+"]";
    end
end

function value=fileSHA256(filename)
    fid=fopen(filename,'rb');
    assert(fid~=-1,'Cannot read input for provenance: %s.',filename);
    cleanup=onCleanup(@()fclose(fid));
    bytes=fread(fid,Inf,'*uint8');
    digest=java.security.MessageDigest.getInstance('SHA-256');
    digest.update(typecast(bytes,'int8'));
    value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end

function [cfg,plan,planFile,names]=initializeAnalysis(cfg,projectFolder)
    % Normalize portable paths once, including Windows-style nested output folders.
    cfg.resultsDirectory=strrep(cfg.resultsDirectory,'\','/');
    assert(isempty(regexp(cfg.resultsDirectory,'^[A-Za-z]:|^/','once')), ...
        'Use a project-relative results directory, not a Windows absolute path.');
    outputFolder = fullfile(projectFolder,cfg.resultsDirectory);
    if ~isfolder(outputFolder); mkdir(outputFolder); end
    planFile = fullfile(projectFolder,'data','analysis_plan.json');
    if ~isfolder(fileparts(planFile)); mkdir(fileparts(planFile)); end
    names = fieldnames(cfg.carbonSources);
    names = names(:)';
    assert(~isempty(names) && numel(unique(names))==numel(names),'Analysis names must be nonempty and unique.');
    assert(all(cellfun(@isvarname,names)),'Use distinct MATLAB identifier names for analyses.');
    % An incomplete export must not leave the previous run's plan usable.
    if isfile(planFile); delete(planFile); end
    plan.version = 3;
    plan.models = struct();
    plan.output_dir = cfg.resultsDirectory;
end

function checkModelInput(model,name)
    % Reference-model compatibility checks; no model construction is performed.
    model=generateRules(model);
    switch name
        case 'Saccharomyces_cerevisiae_WT'
            assert(nnz(model.c)==1 && model.c(strcmp(model.rxns,'r_2111'))==1, ...
                'The source model must use its reference biomass objective.');
            assert(strcmp(string(model.version),'9.1.1'), ...
                'Use the pinned Yeast9 9.1.1 reference model.');
        case 'Saccharomyces_cerevisiae_IUPD'
            assert(strcmp(string(model.version),'9.1.1'), ...
                'Use the pinned Yeast9 9.1.1 reference model.');
            assert(ismember('YML126C',model.genes),'The knockout gene is absent from the source model.');
            assert(all(ismember({'s_0434','s_0394','s_0943'},model.mets)), ...
                'The source model is missing cytosolic ATP, ADP or IPP.');
        case 'Escherichia_coli_WT'
            assert(nnz(model.c)==1 && model.c(strcmp(model.rxns,'BIOMASS_Ec_iML1515_core_75p37M'))==1, ...
                'The source model must use its reference biomass objective.');
            assert(~ismember('SQS',model.rxns),'The source model already contains SQS.');
            assert(all(ismember({'frdp_c','nadph_c','h_c','ppi_c','nadp_c'},model.mets)), ...
                'The source model is missing the configured SQS substrates.');
        case 'Escherichia_coli_IUPD'
            assert(~ismember('SQS',model.rxns),'The source model already contains SQS.');
            assert(all(ismember({'frdp_c','nadph_c','h_c','ppi_c','nadp_c'},model.mets)), ...
                'The source model is missing the configured SQS substrates.');
            assert(ismember('b2515',model.genes),'The knockout gene is absent from the source model.');
            assert(all(ismember({'atp_c','adp_c','ipdp_c'},model.mets)), ...
                'The source model is missing cytosolic ATP, ADP or IPP.');
        case 'Yarrowia_lipolytica_WT'
            assert(nnz(model.c)==1 && model.c(strcmp(model.rxns,'xBIOMASS'))==1, ...
                'The source model must use its reference biomass objective.');
        case 'Yarrowia_lipolytica_IUPD'
            assert(ismember('YALI0F30481g',model.genes),'The knockout gene is absent from the source model.');
            assert(all(ismember({'m141','m143','m491'},model.mets)), ...
                'The source model is missing cytosolic ATP, ADP or IPP.');
    end
end

function checkKnockout(original,deleted)
    assert(any(original.lb~=deleted.lb | original.ub~=deleted.ub), ...
        'The knockout did not constrain any reaction.');
end

function checkSolution(model,solution,name)
    assert(solution.stat==1 && isfinite(solution.f) && solution.f>0 && all(isfinite(solution.x)), ...
        'Optimization failed for %s.',name);
    assert(norm(model.S*solution.x-model.b,inf)<1e-6, ...
        'Steady-state residual exceeds tolerance for %s.',name);
end

function model=closeOtherCarbonUptake(model,glucoseExchange)
    % Retain the same inorganic medium/maintenance bounds and secretion capacity.
    boundaryReactions=find(sum(model.S~=0,1)==1);
    for reactionIndex=boundaryReactions
        if strcmp(model.rxns{reactionIndex},glucoseExchange); continue; end
        metaboliteIndex=find(model.S(:,reactionIndex)~=0);
        if isempty(regexp(model.metFormulas{metaboliteIndex},'C(?![a-z])','once')); continue; end
        if model.S(metaboliteIndex,reactionIndex)<0 && model.lb(reactionIndex)<0
            model=changeRxnBounds(model,model.rxns{reactionIndex},0,'l');
        elseif model.S(metaboliteIndex,reactionIndex)>0 && model.ub(reactionIndex)>0
            model=changeRxnBounds(model,model.rxns{reactionIndex},0,'u');
        end
    end
end
