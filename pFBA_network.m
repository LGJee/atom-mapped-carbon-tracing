%% Export WT and IUPD pFBA-active reactions from Yeast9 9.1.1
% This is a study-specific reference model, not a universal culture setup.
clear; clc;

%% User settings
scriptFolder = fileparts(mfilename('fullpath'));
cfg.modelFolder = fullfile(scriptFolder,'model');
cfg.cobraFolder = ''; % Optional COBRA Toolbox folder; otherwise use the MATLAB path.
cfg.resultsFolder = fullfile(scriptFolder,'results');
cfg.solver = 'glpk';
cfg.expectedModelVersion = '9.1.1';
cfg.glucoseForGrowth = -1;
cfg.oxygenLowerBound = -1000;
cfg.isoprenolLowerBound = -2.09;
cfg.wtGlucoseForProduction = -2;
cfg.growthFraction = 0.5;
cfg.nonzeroTolerance = 1e-8;

%% 1. Load the model and set the reference culture conditions
if ~isempty(cfg.cobraFolder)
    addpath(cfg.cobraFolder);
    initCobraToolbox(false);
end
assert(exist('optimizeCbModel','file')==2, ...
    'Install and initialize COBRA Toolbox following its official instructions.');
assert(changeCobraSolver(cfg.solver,'LP'), 'The configured LP solver is unavailable.');
if ~isfolder(cfg.resultsFolder); mkdir(cfg.resultsFolder); end
data = load(fullfile(cfg.modelFolder,'yeast-GEM.mat'));

modelWT = generateRules(data.model);
assert(strcmp(string(modelWT.version),cfg.expectedModelVersion),'Use the pinned Yeast9 9.1.1 reference model.');
modelWT = changeRxnBounds(modelWT,'r_1714',cfg.glucoseForGrowth,'l');
modelWT = changeRxnBounds(modelWT,'r_1992',cfg.oxygenLowerBound,'l');
modelWT = changeObjective(modelWT,'r_2111');

%% 2. Build IUPD: delete ERG13 and add the isoprenol utilization pathway
modelIUPD = deleteModelGenes(modelWT,{'YML126C'});
modelIUPD = addMetabolite(modelIUPD,'isoprenol[e]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isoprenol[c]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addMetabolite(modelIUPD,'isopentenyl_phosphate[c]', ...
    'metName','isopentenyl phosphate','metFormula','C5H11O4P');
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e','reactionFormula', ...
    'isoprenol[e] -> ','lowerBound',cfg.isoprenolLowerBound,'upperBound',1000);
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e_t','reactionFormula', ...
    'isoprenol[e] -> isoprenol[c]','lowerBound',0,'upperBound',1000);
modelIUPD = addReaction(modelIUPD,'isopentenolK','reactionFormula', ...
    'isoprenol[c] + s_0434 -> isopentenyl_phosphate[c] + s_0394', ...
    'geneRule','YLR133W','lowerBound',0,'upperBound',1000);
modelIUPD = addReaction(modelIUPD,'PisopentenolK','reactionFormula', ...
    'isopentenyl_phosphate[c] + s_0434 -> s_0943 + s_0394', ...
    'geneRule','AtIPK','lowerBound',0,'upperBound',1000);

%% 3. Fix a fraction of maximum growth, maximize squalene, then minimize total flux
modelWT = changeObjective(modelWT,'r_2111');
modelIUPD = changeObjective(modelIUPD,'r_2111');

growthWT = optimizeCbModel(modelWT,'max');
growthIUPD = optimizeCbModel(modelIUPD,'max');
assert(growthWT.stat==1 && isfinite(growthWT.f), 'WT growth optimization failed.');
assert(growthIUPD.stat==1 && isfinite(growthIUPD.f), 'IUPD growth optimization failed.');
maxGrowthWT = growthWT.f;
maxGrowthIUPD = growthIUPD.f;

modelWT = changeRxnBounds(modelWT,'r_2111',cfg.growthFraction*maxGrowthWT,'b');
modelIUPD = changeRxnBounds(modelIUPD,'r_2111',cfg.growthFraction*maxGrowthIUPD,'b');

% Reference study: growth is optimized at -1, then WT glucose availability increases to -2.
% Negative exchange flux denotes uptake; changing the lower bound does not force uptake.
modelWT = changeRxnBounds(modelWT,'r_1714',cfg.wtGlucoseForProduction,'l');

modelWT = changeObjective(modelWT,'r_1012');
modelIUPD = changeObjective(modelIUPD,'r_1012');

% 'one' preserves the optimal objective and minimizes the sum of absolute fluxes.
solutionWT = optimizeCbModel(modelWT,'max','one');
solutionIUPD = optimizeCbModel(modelIUPD,'max','one');
assert(solutionWT.stat==1 && all(isfinite(solutionWT.x)), 'WT pFBA failed.');
assert(solutionIUPD.stat==1 && all(isfinite(solutionIUPD.x)), 'IUPD pFBA failed.');

%% 4. Collect nonzero reactions and their participating metabolites
nonzeroWT = makeNonzeroTable(modelWT,solutionWT.x,cfg.nonzeroTolerance);
nonzeroIUPD = makeNonzeroTable(modelIUPD,solutionIUPD.x,cfg.nonzeroTolerance);
[yamlIDs,yamlSmiles] = readYAMLSmiles(fullfile(cfg.modelFolder,'yeast-GEM.yml'));
metabolitesWT = makeMetaboliteTable(modelWT,solutionWT.x,yamlIDs,yamlSmiles,cfg.nonzeroTolerance);
metabolitesIUPD = makeMetaboliteTable(modelIUPD,solutionIUPD.x,yamlIDs,yamlSmiles,cfg.nonzeroTolerance);

%% 5. Write a separate workbook for each strain
wtFile = fullfile(cfg.resultsFolder,'WT_pFBA_nonzero_reactions.xlsx');
iupdFile = fullfile(cfg.resultsFolder,'IUPD_pFBA_nonzero_reactions.xlsx');
if isfile(wtFile); delete(wtFile); end
if isfile(iupdFile); delete(iupdFile); end
writetable(nonzeroWT,wtFile,'Sheet','WT_nonzero');
writetable(metabolitesWT,wtFile,'Sheet','Metabolites');
writetable(nonzeroIUPD,iupdFile,'Sheet','IUPD_nonzero');
writetable(metabolitesIUPD,iupdFile,'Sheet','Metabolites');


fprintf('Results saved in:\n%s\n',cfg.resultsFolder);

%% Nonzero reactions; equations are generated with COBRA printRxnFormula
function T = makeNonzeroTable(model,v,tolerance)
    rows = find(abs(v)>tolerance);
    namedModel = model;
    namedModel.metNames = cellstr(metLabels(model,(1:numel(model.mets))'));
    formula = string(printRxnFormula(namedModel,'rxnAbbrList',model.rxns(rows), ...
        'printFlag',false,'metNameFlag',true));
    formula = replace(formula," -> "," => ");
    T = table(string(model.rxns(rows)),string(model.rxnNames(rows)),v(rows),formula, ...
        'VariableNames',{'ReactionID','ReactionName','Flux','Equation'});
end


%% Metabolites participating in the exported nonzero reactions
function T = makeMetaboliteTable(model,v,yamlIDs,yamlSmiles,tolerance)
    activeRxns = find(abs(v)>tolerance);
    rows = find(any(model.S(:,activeRxns)~=0,2));
    metaboliteID = string(model.mets(rows));
    metaboliteName = metLabels(model,rows);

    smilesFromYeast9 = strings(numel(rows),1);
    [found,location] = ismember(metaboliteID,yamlIDs);
    smilesFromYeast9(found) = yamlSmiles(location(found));

    manualSmiles = strings(numel(rows),1);
    manualSmiles(ismember(metaboliteID,["isoprenol[e]","isoprenol[c]"])) = "CC(=C)CCO";
    manualSmiles(metaboliteID=="isopentenyl_phosphate[c]") = "CC(=C)CCOP(=O)(O)O";
    manualSmiles(metaboliteID=="s_0555") = ...
        "C([C@@H]1[C@H]([C@@H]([C@](O1)(COP(=O)([O-])[O-])O)O)O)OP(=O)([O-])[O-]";

    T = table(metaboliteID,metaboliteName,smilesFromYeast9,manualSmiles, ...
        'VariableNames',{'MetaboliteID','MetaboliteName','SMILES_Yeast9','SMILES_Manual'});
end

%% Read metabolite IDs and the first SMILES annotation from the pinned Yeast9 YAML
function [ids,smiles] = readYAMLSmiles(yamlFile)
    lines = readlines(yamlFile);
    firstMet = find(lines=="- metabolites:",1);
    firstRxn = find(lines=="- reactions:",1);
    lines = lines(firstMet+1:firstRxn-1);

    ids = strings(0,1); smiles = strings(0,1); current = 0;
    for k = 1:numel(lines)
        idToken = regexp(lines(k),'^    - id: (.+)$','tokens','once');
        if ~isempty(idToken)
            ids(end+1,1) = string(idToken{1}); %#ok<AGROW>
            smiles(end+1,1) = ""; %#ok<AGROW>
            current = numel(ids);
        elseif current>0 && strcmp(lines(k),'      - smiles:') && k<numel(lines)
            smileToken = regexp(lines(k+1),'^        - (.+)$','tokens','once');
            if ~isempty(smileToken)
                smiles(current) = strrep(string(smileToken{1}),'"','');
            end
        end
    end
end
%% Compartment-qualified metabolite names, for example UDP-D-glucose[c]
function label = metLabels(model,rows)
    rows = rows(:); compartment = strings(numel(rows),1);
    valid = model.metComps(rows)>=1 & model.metComps(rows)<=numel(model.comps);
    compartment(valid) = string(model.comps(model.metComps(rows(valid))));
    ids = string(model.mets(rows));
    for k = find(~valid)'
        token = regexp(char(ids(k)),'\[([^\]]+)\]$','tokens','once');
        if ~isempty(token); compartment(k)=string(token{1}); else; compartment(k)="?"; end
    end
    label = string(model.metNames(rows))+"["+compartment+"]";
end
