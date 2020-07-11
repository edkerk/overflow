% This script modifies protein content, scales remaining biomass,
% recalculates GAM, and subsequentially runs random sampling of these
% condition specific samples. Selected fluxes are written, together with a
% table of all random sampling results.

prepareEnvironment
clear code ecModel* grouping oxPhos positionsEC repl params prot
model = importModel('../models/yeastGEM_v8.3.4.xml');
model.mets = regexprep(model.mets,'\[(.*)\]$','');

% Set some parameters
GAM = 34;
exchRxns = {'r_1714' 'r_1992' 'r_1672' 'r_1808' 'r_1634' 'r_1761' 'r_1793' 'r_2111'};
model = setParam(model,'eq','r_0445',0); % Block formate -> CO2
model = setParam(model,'ub','r_4046',Inf);
goodRxns = [];

% Modify models and perform random sampling
for i=1:length(flux.conds)
%     cd GECKO/geckomat/limit_proteins
%     disp(['Perform random sampling for condition: ' flux.conds{i}])
%     % Scale biomass to updated protein content
%     models{i}=scaleBioMass(model,flux.Ptot(i),GAM,true);
%     % Set exchange reactions to 95-105% of measured values
%     meas = [flux.GUR(i) flux.OxyUptake(i) flux.CO2prod(i) flux.byP_flux(i,:) flux.Drate(i)];
%     LBs = [-1.05 -1.05 0.95 0.95 0.95 0.95 0.95 0.95] .* meas;
%     UBs(UBs == 0) = 0.01;
%     UBs = [-0.95 -0.95 1.05 1.05 1.05 1.05 1.05 1.05] .* meas;
%     
%     models{i}=setParam(models{i},'lb',exchRxns,LBs);
%     models{i}=setParam(models{i},'ub',exchRxns,UBs);
%     % Constrain model to 95% of max NGAM
%     models{i}=setParam(models{i},'obj','r_4046',1);
%     sol=solveLP(models{i});
%     models{i}=setParam(models{i},'lb','r_4046',0.95*sol.x(getIndexes(models{i},'r_4046','rxns')));
%     % Random sampling
%     cd ../../../
%     [sols{i},goodRxns]=rs(models{i},5000,true,true,true,goodRxns,true);
%     
    out.mean(:,i)=full(mean(sols{i},2));
    out.std(:,i)=full(std(sols{i},0,2));
end

%% Write selected reaction fluxes, yields and turnover
fluxRxns={'r_0226','r_1022','r_0892','r_0962','r_4235','r_0886','r_0534',...
    'r_4046','r_4041','r_1166','r_0961','r_0658','r_0714','r_0713','r_0770'};
idxs=getIndexes(model,fluxRxns,'rxns');

% Numbers to be used as input for S4B
rGlu=out.mean(idxs(10),:);
fluxes(1,:)=rGlu;
fluxes(2,:)=sum(out.mean(idxs(1:2),:)); % ETC rATP
fluxes(3,:)=fluxes(2,:)./rGlu; % ETC YATP
fluxes(4,:)=sum(out.mean(idxs(3:4),:))-sum(out.mean(idxs(5:7),:)); % Glycolysis rATP
fluxes(5,:)=fluxes(4,:)./rGlu; % Glycolysis YATP
for i=1:numel(models); fluxes(6,i)=models{i}.S(getIndexes(models{i},'ATP[c]','metcomps'),idxs(9)); end;
fluxes(6,:)=-fluxes(6,:).*out.mean(idxs(9),:); % GAEC rATP
fluxes(7,:)=fluxes(6,:)./rGlu; % GAEC YATP
fluxes(8,:)=out.mean(idxs(8),:); % NGAM rATP
fluxes(9,:)=fluxes(8,:)./rGlu; % NGAM YATP
fluxes(10,:)=fluxes(6,:)+fluxes(8,:); % NGAM+GAEC rATP
fluxes(11,:)=fluxes(7,:)+fluxes(9,:); % NGAM+GAEC YATP
fluxes(12:16,:)=out.mean(idxs(11:15),:); % rPDH, rIDH, rMDHc, rMDHm, rNDE

% Rates related to NAD turnover
NADmets=find(strcmp(model.metNames,'NAD'));
[~,NADrxns]=find(model.S(NADmets,:));%>0);
for i=1:numel(models)
    NADtmp=full(models{i}.S(NADmets,NADrxns).*out.mean(NADrxns,i)');
    NADtmp(NADtmp<0)=0;
    NAD(:,i)=sum(NADtmp,2);
end
zeroFlux=sum(NAD,2)==0;
NAD=[NAD(~zeroFlux,:);sum(NAD,1)];
NADmets=[strcat('rNAD[',model.comps(model.metComps(NADmets(~zeroFlux))),']'); ...
    'rNAD[tot]']';

fluxes=[fluxes;NAD];

rowNames=[{'rGlu','ETC_rATP','ETC_YATP','glycolysis_rATP','glycolysis_YATP',...
    'GAEC_rATP','GAEC_YATP','NGAM_rATP','NGAM_YATP','GAEC+NGAM_rATP',...
    'GAEC+NGAM_YATP','rPDH','rIDH','rMDHc','rMDHm','rNDE'}, NADmets];
varNames={'CN4','CN22','CN38','CN78','hGR'};

fluxes=cell2table(num2cell(fluxes),'RowNames',rowNames,'VariableNames',varNames);   
writetable(fluxes,'../results/randomSampling/selectedFluxes.txt','WriteRowNames',true,'Delimiter','\t')

%% Write all random sampling results (means and standard deviations)
outTable(:,1)=model.rxns;
outTable(:,2)=model.rxnNames;
outTable(:,3)=constructEquations(model);
outTable(:,4)=model.eccodes;
outTable(:,5)=model.grRules;
outTable(:,[6,9,12,15,18])=num2cell(out.mean);
outTable(:,[7,10,13,16,19])=num2cell(out.std);
outTable(:,[8,11,14,17,20])=num2cell((out.std./out.mean)*100);

varNames={'ID','NAME','EQUATION','EC-NUMBER','GENE ASSOCIATION',...
    'CN4_AVERAGE','CN4_STDEV','CN4_STDEV, %','CN22_AVERAGE','CN22_STDEV',...
    'CN22_STDEV, %','CN38_AVERAGE','CN38_STDEV','CN38_STDEV, %',...
    'CN78_AVERAGE','CN78_STDEV','CN78_STDEV, %','hGR_AVERAGE',...
    'hGR_STDEV','hGR_STDEV, %'};

outTable=cell2table(outTable,'VariableNames',varNames);
writetable(outTable,'../results/randomSampling/allFluxes.txt','Delimiter','\t')
