% Analyze enzymeUsage per subsystem
% Assumed is that models are already constructed by generateProtModels
% and ribosome.m is run to add ribosomal subunits.

%% Get enzymes usages to each reaction
for i=1:5
    disp(['Now testing: ' flux.conds{i}])
    sol{i} = solveLP(ecModels{i});
    [absUsage{i}, capUsage{i}, UB{i}, protName{i}]=enzymeUsage(ecModels{i},sol{i}.x,true);
    printFluxes(ecModels{i},sol{i}.x,false,0,fullfile('..','results','modelSimulation',['allFluxes_' flux.conds{i},'.txt']),'%rxnID\t%rxnName\t%eqn\t%flux\n');
end

%% Map enzymes to reactions and their subsystems
clear out usages
[Lia, Locb]=ismember(ecModels{1}.enzGenes,ecModels{1}.genes);
usages.enzymes = char.empty;
usages.subSys = char.empty;
usages.absUse = double.empty;
usages.capUse = double.empty;

for j=1:length(Locb);
    rxnIdx = find(ecModels{1}.rxnGeneMat(:,Locb(j)));
    subSys = ecModels{1}.subSystems(rxnIdx);
    subSys = flattenCell(subSys);
    subSys = reshape(subSys,numel(subSys),1);
    subSys(cellfun(@isempty,subSys)) = [];
    subSys = unique(subSys);
    if isempty(subSys); subSys = ' '; end;
    nSub = numel(subSys);
    enzymes = repmat(ecModels{1}.enzymes(j),nSub,1);
    genes = repmat(ecModels{1}.enzGenes(j),nSub,1);
    usages.enzymes = [usages.enzymes; enzymes];
    usages.subSys = [usages.subSys; subSys];
end
usages.subSys=regexprep(usages.subSys,'^sce\d{5}  ','');
if exist('ribo','var')
    usages.subSys(end+1:end+numel(ribo))={'Ribosome (main subunits)'};
    usages.enzymes(end+1:end+numel(ribo))=ribo;
end

%% Populate with usage information
for i=1:numel(ecModels)
    enz = protName{i};
    [Lia, Locb]=ismember(usages.enzymes,enz);
    usages.absUse(:,i) = absUsage{i}(Locb);
    usages.capUse(:,i) = capUsage{i}(Locb);
    usages.UB(:,i)     = UB{i}(Locb);
end

%% All usage per subSystem
head={'protID','subSys','CN4','CN22','CN38','CN78','hGR'};
all=cell2table([usages.enzymes,usages.subSys,num2cell(usages.capUse)],'VariableNames',head);
writetable(all,fullfile('..','results','enzymeUsage','enzymeCapUsages.txt'),'Delimiter','\t')