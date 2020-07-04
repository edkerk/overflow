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
% save('../models/ecYeastGEM.mat','ecModel','ecModel_batch');
load('../models/ecYeastGEM.mat');
cd(code)
%% Prepare data
%Load proteomics data
fID       = fopen('../data/abs_proteomics.txt');
prot.cond = textscan(fID,['%s' repmat(' %s',1,17)],1);
prot.data = textscan(fID,['%s %s' repmat(' %f',1,16)],'TreatAsEmpty',{'NA','na','NaN'});
prot.cond = [prot.cond{3:end}];
prot.IDs  = prot.data{1};
prot.data = cell2mat(prot.data(3:end));
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
clear repl
for i=1:length(grouping)
    try
        repl.first(i)=repl.last(end)+1;
        repl.last(i)=repl.last(end)+grouping(i);
    catch
        repl.first=1;
        repl.last=grouping(1);
    end
end
%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = find(strcmpi(ecModel.rxnNames,params.c_source));
positionsEC(2) = find(strcmpi(ecModel.rxns,params.bioRxn));
clear ans fID data byProds fileNames GECKO_path i fID