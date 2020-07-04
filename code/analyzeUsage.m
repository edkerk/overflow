% Analyze enzymeUsage per subsystem
% Assumed is that models are already constructed by generateProtModels
% and ribosome.m is run to add ribosomal subunits.

%% Get enzymes usages to each reaction
for i=1:5
    disp(['Now testing: ' flux.conds{i}])
    sol{i} = solveLP(ecModels{i});
    [capUsage{i}, absUsage{i}, protName{i}]=enzymeUsage(ecModels{i},sol{i}.x,true);
    printFluxes(ecModels{i},sol{i}.x,false,0,fullfile('..','results','modelSimulation',['allFluxes_' flux.conds{i},'.txt']),'%rxnID\t%rxnName\t%eqn\t%flux\n');
end

%% Map enzymes to reactions and their subsystems
clear out
[Lia, Locb]=ismember(ecModels{1}.enzGenes,ecModels{1}.genes);
out.enzymes = char.empty;
out.subSys = char.empty;
out.absUse = double.empty;
out.capUse = double.empty;

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
    out.enzymes = [out.enzymes; enzymes];
    out.subSys = [out.subSys; subSys];
end
out.subSys=regexprep(out.subSys,'^sce\d{5}  ','');
if exist('ribo','var')
    out.subSys(end+1:end+numel(ribo))={'Ribosome (main subunits)'};
    out.enzymes(end+1:end+numel(ribo))=ribo;
end

%% Populate with usage information
for i=1:numel(ecModels)
    enz = protName{i};
    [Lia, Locb]=ismember(out.enzymes,enz);
    out.absUse(:,i) = absUsage{i}(Locb);
    out.capUse(:,i) = capUsage{i}(Locb);
end

%% All usage per subSystem
head={'protID','subSys','CN04','CN22','CN38','CN78','hGR'};
all=cell2table([out.enzymes,out.subSys,num2cell(out.capUse)],'VariableNames',head);
writetable(all,fullfile('..','results','enzymeUsage','enzymeCapUsages.txt'),'Delimiter','\t')
all=cell2table([out.enzymes,out.subSys,num2cell(out.absUse)],'VariableNames',head);
writetable(all,fullfile('..','results','enzymeUsage','enzymeAbsUsages.txt'),'Delimiter','\t')

%% Summarize per subSystem
subSysUse.subSys = unique(out.subSys);
[Lia, Locb]=ismember(out.subSys,subSysUse.subSys);
for i=1:numel(ecModels)
    for j=1:numel(subSysUse.subSys)
        subSysUse.capUse(j,i) = mean(out.capUse(Locb==j,i));
        subSysUse.absUse(j,i) = mean(out.absUse(Locb==j,i));
        subSysUse.over50use(j,i) = numel(find(out.capUse(Locb==j,i)>0.50));
    end
end
subSysUse.subSys=regexprep(subSysUse.subSys,'^sce\d{5}  ','');

out=[subSysUse.subSys, num2cell(subSysUse.capUse), num2cell(subSysUse.absUse)];
head = {'subSys', 'capUse_CN04', 'capUse_CN22', 'capUse_CN38', 'capUse_CN78', ...
    'capUse_hGR', 'absUse_CN04', 'absUse_CN22', 'absUse_CN38', 'absUse_CN78', ...
    'absUse_hGR'};
out=cell2table(out,'VariableNames',head);
writetable(out,fullfile('..','results','enzymeUsage','subSysUsage.txt'),'Delimiter','\t')