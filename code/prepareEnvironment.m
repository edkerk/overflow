% This function prepares the environment to run the analysis, by cloning
% the correct repositories and setting appropriate parameters

code = pwd();

%% Prepare software and model
% Clone GECKO and checkout version 2.0.0. GECKO will not be tracked in the
% overflow repository, and this should therefore be run when cloning the
% overflow repository.
if ~exist('GECKO', 'dir')
    git('clone https://github.com/SysBioChalmers/GECKO.git')
    cd GECKO
    git('checkout tags/v2.0.0 -b overflow')
    cd(code)
end

% Replace custom GECKO scripts
fileNames = struct2cell(dir('customGECKO'));
fileNames = fileNames(1,:);
fileNames(startsWith(fileNames,'.')) = [];
for i = 1:length(fileNames)
    GECKO_path = dir(['GECKO/**/' fileNames{i}]);
    copyfile(['customGECKO' filesep fileNames{i}],GECKO_path.folder)
    disp(['Replaced ' fileNames{i} ' at ' GECKO_path.folder '\'])
end

% % Get yeast-GEM 8.1.3 if ecYeast-GEM is first to be reconstructed.
% if ~exist('yeast-GEM', 'dir')
%     git('clone https://github.com/SysBioChalmers/yeast-GEM.git')
%     cd yeast-GEM
%     git('checkout tags/v8.1.3 -b overflow')
%     cd(code)
%     copyfile yeast-GEM/ModelFiles/xml/yeastGEM.xml ../models/yeastGEM.xml
% end

% % Reconstruct ecYeast-GEM 8.1.3
% model = importModel('../models/yeastGEM.xml');
% model.mets = regexprep(model.mets,'\[(.*)\]$','');
% cd GECKO/geckomat
% [ecModel,ecModel_batch] = enhanceGEM(model,'RAVEN');
% cd(code)
% save('../models/ecYeastGEM.mat',{'ecModel','ecModel_batch'});
load('../models/ecYeastGEM.mat');

% % % Instead of ecYeast-GEM 8.1.3, use 8.3.5 with the following code.
% % % Clone ecYeastGEM and checkout ecYeastGEM 8.3.5. The model will be kept in
% % % the overflow repository, so this code only needs to be run (and modified)
% % % when a new version of ecYeastGEM will be used.
% % git('clone https://github.com/SysBioChalmers/ecModels.git')
% % cd ecModels
% % git('checkout 574fa93d39c78f6206cee55770ad99d5c5828d7b -b overflow')
% % copyfile ecYeastGEM/model/ecYeastGEM_batch.xml ../../models/
% % copyfile ecYeastGEM/model/ecYeastGEM.mat ../../models/
% % copyfile ecYeastGEM/model/ecYeastGEM_batch.mat ../../models/
% % mkdir ../../models/prot_constrained/
% % cd(code)
% % load ('../models/ecYeastGEM.mat')
% % load ('../models/ecYeastGEM_batch.mat')
% % ecYeastGEM_batch.mets = regexprep(ecYeastGEM_batch.mets,'\[(.*)\]$','');
% % ecYeastGEM.mets = regexprep(ecYeastGEM.mets,'\[(.*)\]$','');
cd(code)
%% Prepare data
%Load proteomics data
fID       = fopen('../data/abs_proteomics.txt');
prot.cond = textscan(fID,['%s' repmat(' %s',1,17)],1);
prot.data = textscan(fID,['%s %s' repmat(' %f',1,16)],'TreatAsEmpty',{'NA','na','NaN'});
prot.cond = [prot.cond{3:end}];
prot.IDs  = prot.data{1};
prot.data = prot.data(3:end);
fclose(fID);

%Load total protein content and fermentation data
fID       = fopen('../data/fermentationData.txt');
byProds   = textscan(fID,['%s' repmat(' %s',1,9)],1,'Delimiter','\t');
data      = textscan(fID,['%s' repmat(' %f',1,9)],'TreatAsEmpty',{'NA','na','NaN'});
fclose(fID);

flux.conds     = data{1};
flux.Ptot      = data{2};
flux.Drate     = data{3};
flux.GUR       = data{4};
flux.CO2prod   = data{5};
flux.OxyUptake = data{6};
flux.byP_flux  = [data{7:end}];
flux.byP_flux(isnan(flux.byP_flux))=0;
flux.byProds   = [byProds{7:end}];

%% Load model and additional parameters
cd GECKO/geckomat
params = getModelParameters;
cd(code)

%Set some additional parameters
oxPhos = ecModel.rxns(startsWith(ecModel.rxns,params.oxPhos));
grouping=[3 3 3 3 4];
%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = find(strcmpi(ecModel.rxnNames,params.c_source));
positionsEC(2) = find(strcmpi(ecModel.rxns,params.bioRxn));
%Remove prot_abundance.txt and relative_proteomics.txt files
%(for f factor calculation)
cd GECKO/geckomat/limit_proteins/
f       = measureAbundance(ecModel.enzymes); %Protein mass in model/Total theoretical proteome
cd(code)
clear ans fID data byProds fileNames GECKO_path i fID