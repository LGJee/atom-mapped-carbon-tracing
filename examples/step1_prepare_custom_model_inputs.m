%% Example: aerobic glucose growth pFBA with no input SMILES
% Initialize COBRA Toolbox before running this script. Other models need their
% own exchange IDs, medium and objective; the defaults below are for iML1515.
clear; clc;
projectFolder = fileparts(fileparts(mfilename('fullpath')));

%% User settings
cfg.modelFile = fullfile(projectFolder,'data','ecoli','iML1515.mat');
cfg.annotationFile = fullfile(projectFolder,'data','ecoli','iML1515.json'); % Optional BiGG JSON.
cfg.modelName = 'Escherichia_coli';
cfg.outputFolder = fullfile(projectFolder,'results');
cfg.solver = 'glpk';
cfg.glucoseExchange = 'EX_glc__D_e';
cfg.oxygenExchange = 'EX_o2_e';
cfg.growthReaction = 'BIOMASS_Ec_iML1515_core_75p37M';
cfg.glucoseLowerBound = -10;
cfg.oxygenLowerBound = -1000;
cfg.nonzeroTolerance = 1e-8;

%% 1. Load the model, preserving its inorganic medium and maintenance settings
assert(exist('optimizeCbModel','file')==2,'Install and initialize COBRA Toolbox first.');
assert(changeCobraSolver(cfg.solver,'LP'),'The configured LP solver is unavailable.');
data = load(cfg.modelFile);
fields = fieldnames(data);
models = fields(cellfun(@(name)isstruct(data.(name)) && isfield(data.(name),'S') ...
    && isfield(data.(name),'rxns') && isfield(data.(name),'mets'),fields));
assert(numel(models)==1,'The MAT file must contain exactly one COBRA model.');
model = data.(models{1});
assert(all(ismember({cfg.glucoseExchange,cfg.oxygenExchange,cfg.growthReaction},model.rxns)), ...
    'Configure exchange and objective IDs for the selected model.');

% Close other carbon-containing boundary uptake; secretion remains possible.
boundary = find(sum(model.S~=0,1)==1);
for reaction = boundary
    metabolite = find(model.S(:,reaction)~=0);
    if strcmp(model.rxns{reaction},cfg.glucoseExchange); continue; end
    if ~isempty(regexp(model.metFormulas{metabolite},'C(?![a-z])','once'))
        if model.S(metabolite,reaction)<0
            model.lb(reaction) = max(0,model.lb(reaction));
        else
            model.ub(reaction) = min(0,model.ub(reaction));
        end
    end
end
model = changeRxnBounds(model,cfg.glucoseExchange,cfg.glucoseLowerBound,'l');
model = changeRxnBounds(model,cfg.oxygenExchange,cfg.oxygenLowerBound,'l');
model = changeObjective(model,cfg.growthReaction);

%% 2. Maximize growth, then minimize total absolute flux
solution = optimizeCbModel(model,'max','one');
assert(solution.stat==1 && all(isfinite(solution.x)),'Growth pFBA failed.');
assert(norm(model.S*solution.x-model.b,inf)<1e-6,'Steady-state residual exceeds tolerance.');
active = find(abs(solution.x)>cfg.nonzeroTolerance);
participating = find(any(model.S(:,active)~=0,2));

%% 3. Export reactions and metabolite identities, without retrieving structures
names = qualifiedNames(model);
namedModel = model;
namedModel.metNames = cellstr(names);
equations = string(printRxnFormula(namedModel,'rxnAbbrList',model.rxns(active), ...
    'printFlag',false,'metNameFlag',true));
reactions = table(string(model.rxns(active)),string(model.rxnNames(active)), ...
    solution.x(active),equations,'VariableNames',{'ReactionID','ReactionName','Flux','Equation'});
inchiKeys = strings(numel(model.mets),1);
if ~isempty(cfg.annotationFile)
    annotations = jsondecode(fileread(cfg.annotationFile));
    for entry = annotations.metabolites(:)'
        position = find(strcmp(model.mets,entry.id));
        if ~isempty(position) && isfield(entry.annotation,'inchi_key')
            keys = string(entry.annotation.inchi_key);
            if numel(keys)==1; inchiKeys(position)=keys; end
        end
    end
end
metabolites = table(string(model.mets(participating)),names(participating), ...
    string(model.metFormulas(participating)),model.metCharge(participating),inchiKeys(participating), ...
    'VariableNames',{'MetaboliteID','MetaboliteName','MolecularFormula','Charge','InChIKey'});
if ~isfolder(cfg.outputFolder); mkdir(cfg.outputFolder); end
output = fullfile(cfg.outputFolder,[cfg.modelName '_1_growth_pfba.xlsx']);
if isfile(output); delete(output); end
writetable(reactions,output,'Sheet','Nonzero');
writetable(metabolites,output,'Sheet','Metabolites');
fprintf('%s: growth %.9g; %d nonzero reactions, %d metabolites; no SMILES supplied.\n', ...
    cfg.modelName,solution.f,height(reactions),height(metabolites));
fprintf('Input saved: %s\n',output);

function names = qualifiedNames(model)
    names = strings(numel(model.mets),1);
    for k=1:numel(model.mets)
        if isfield(model,'metComps') && isfield(model,'comps')
            compartment = string(model.comps{model.metComps(k)});
        else
            token = regexp(model.mets{k},'\[([^\]]+)\]$','tokens','once');
            if isempty(token); token=regexp(model.mets{k},'_([^_]+)$','tokens','once'); end
            assert(~isempty(token),'A metabolite has no identifiable compartment: %s',model.mets{k});
            compartment = string(token{find(~cellfun(@isempty,token),1)});
        end
        name = string(model.metNames{k});
        if strcmp(name,'D-Glucose'); name="D-glucose"; end
        names(k) = name+"["+compartment+"]";
    end
end
