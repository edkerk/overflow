% Analyze enzymeUsage per subsystem
% Assumed is that models are already constructed by generateProtModels
% and ribosome.m is run to add ribosomal subunits.

load('../models/ecModel_P_CN4.mat')
load('../models/ecModel_P_CN22.mat')
load('../models/ecModel_P_CN38.mat')
load('../models/ecModel_P_CN78.mat')
load('../models/ecModel_P_hGR.mat')

ecModels{1}=ecModelP_CN4;
ecModels{2}=ecModelP_CN22;
ecModels{3}=ecModelP_CN38;
ecModels{4}=ecModelP_CN78;
ecModels{5}=ecModelP_hGR;

%% Get enzymes usages to each reaction
for i=1:5
    disp(['Now testing: ' flux.conds{i}])
    sol{i} = solveLP(ecModels{i});
    [absUsage{i}, capUsage{i}, UB{i}, protName{i}]=enzymeUsage(ecModels{i},sol{i}.x,true);
    printFluxes(ecModels{i},sol{i}.x,false,0,fullfile('..','results','modelSimulation',['allFluxes_' flux.conds{i},'.txt']),'%rxnID\t%rxnName\t%eqn\t%flux\n');
end

%% Prepare output
clear out
out(:,1)=ecModels{1}.enzymes;
out(:,2)=ecModels{1}.enzGenes;
out(:,3)=ecModels{1}.enzNames;
for i=1:5
    out(:,3+i)=strtrim(cellstr(num2str(capUsage{i},3)));
end
for i=1:5
    out(:,8+i)=strtrim(cellstr(num2str(absUsage{i},3)));
end
for i=1:5
    out(:,13+i)=strtrim(cellstr(num2str(UB{i},3)));
end

%% All usage per subSystem
head={'protID','geneID','protName','capUse_CN4','capUse_CN22','capUse_CN38','capUse_CN78',...
    'capUse_hGR','absUse_CN4','absUse_CN22','absUse_CN38','absUse_CN78','absUse_hGR',...
    'UB_CN4','UB_CN22','UB_CN38','UB_CN78','UB_hGR'};
out=cell2table(out,'VariableNames',head);
writetable(out,fullfile('..','results','enzymeUsage','enzymeUsages.txt'),'Delimiter','\t')