% This script modifies protein content, scales remaining biomass,
% recalculates GAM, and subsequentially runs random sampling of these
% condition specific samples. Selected fluxes are written, together with a
% table of all random sampling results.

prepareEnvironment
clear ecModel* grouping oxPhos positionsEC repl prot
model = importModel('../models/yeastGEM_v8.3.4.xml');
model.mets = regexprep(model.mets,'\[(.*)\]$','');

% Set some parameters
GAM = 34;
exchRxns = {'r_1714' 'r_1992' 'r_1672' 'r_1808' 'r_1634' 'r_1761' 'r_1793' 'r_2111'};
model = setParam(model,'eq','r_0445',0); % Block formate -> CO2
model = setParam(model,'ub','r_4046',Inf);
goodRxns = [];
goodRxnsNoForm = [];

% Modify models and perform random sampling
for i=1:length(flux.conds)
    cd GECKO/geckomat/limit_proteins
    disp(['Perform random sampling for condition: ' flux.conds{i}])
    % Scale biomass to updated protein content
    models{i}=scaleBioMass(model,flux.Ptot(i),GAM,true);
    [~,Ptot,Ctot,Rtot,Dtot] = sumBioMass(models{i});
    GAECpol(1:4,i) = [Ptot Ctot Rtot Dtot] .* params.pol_cost;
    GAECpol(5,i) = sum(GAECpol(1:4,i));
    GAECpol(6:7,i) = [34; GAECpol(5,i) + 34];    
    % Set exchange reactions to 95-105% of measured values
    meas = [flux.GUR(i) flux.OxyUptake(i) flux.CO2prod(i) flux.byP_flux(i,:) flux.Drate(i)];
    LBs = [-1.05 -1.05 0.95 0.95 0.95 0.95 0.95 0.95] .* meas;
    UBs = [-0.95 -0.95 1.05 1.05 1.05 1.05 1.05 1.05] .* meas;
    UBs(UBs == 0) = 0.01; % Allow small amount of non detected overflow metabolite
    
    % First run random sampling without formate exchange
    models{i}=setParam(models{i},'lb',exchRxns([1:6,8]),LBs([1:6,8]));
    models{i}=setParam(models{i},'ub',exchRxns([1:6,8]),UBs([1:6,8]));
    % Constrain model to 95% of max NGAM
    models{i}=setParam(models{i},'obj','r_4046',1);
    sol=solveLP(models{i});
    models{i}=setParam(models{i},'lb','r_4046',0.95*sol.x(getIndexes(models{i},'r_4046','rxns')));
    % Random sampling. minFlux = false, as internal loops are fine.
    cd(code)
    [solsnoFrom{i},goodRxnsNoForm]=rs(models{i},5000,true,true,true,goodRxnsNoForm,false);
    
    outNoForm(:,i)=full(mean(solsnoFrom{i},2));
    
    % Rerun random sampling with measured formate
    models{i}=setParam(models{i},'lb',exchRxns,LBs);
    models{i}=setParam(models{i},'ub',exchRxns,UBs);
    % Constrain model to 95% of max NGAM
    models{i}=setParam(models{i},'lb','r_4046',0);
    models{i}=setParam(models{i},'obj','r_4046',1);
    sol=solveLP(models{i});
    models{i}=setParam(models{i},'lb','r_4046',0.95*sol.x(getIndexes(models{i},'r_4046','rxns')));
    % Random sampling
    [sols{i},goodRxns]=rs(models{i},5000,true,true,true,goodRxns,true);
    
    out.mean(:,i)=full(mean(sols{i},2));
    out.std(:,i)=full(std(sols{i},0,2));
end

%% Predicted overflow metabolites from random sampling (result: formate)
[~,idx]=getExchangeRxns(models{1});
exch=find(sum(outNoForm(idx,:)>0.001,2));
idx=idx(exch);

varNames={'ID','Reaction name','CN4','CN22','CN38','CN78','hGR'};
exch=cell2table([models{1}.rxns(idx), models{1}.rxnNames(idx), ...
    num2cell(outNoForm(idx,:))], 'VariableNames', varNames);
writetable(exch,'../results/randomSampling/altExchangeFlux.txt','Delimiter','\t')

%% Write selected reaction fluxes, yields and turnover
fluxRxns={'r_0226','r_1022','r_0892','r_0962','r_4235','r_0886','r_0534',...
    'r_4046','r_4041','r_1166','r_0961','r_0658','r_0714','r_0713','r_0770'};
idxs=getIndexes(model,fluxRxns,'rxns');

% Numbers to be used as input for S4B
% ATP production rates and yields
rGlu=out.mean(idxs(10),:);
fluxes(1,:)=rGlu;
fluxes(2,:)=sum(out.mean(idxs(1:2),:)); % ETC rATP
fluxes(3,:)=fluxes(2,:)./rGlu; % ETC YATP_glu
fluxes(4,:)=fluxes(2,:)./flux.Drate'; % ETC YATP_mu
fluxes(5,:)=sum(out.mean(idxs(3:4),:))-sum(out.mean(idxs(5:7),:)); % Glycolysis rATP
fluxes(6,:)=fluxes(5,:)./rGlu; % Glycolysis YATP_glu
fluxes(7,:)=fluxes(5,:)./flux.Drate'; % Glycolysis YATP_mu
for i=1:numel(models); fluxes(8,i)=models{i}.S(getIndexes(models{i},'ATP[c]','metcomps'),idxs(9)); end;
fluxes(8,:)=-fluxes(8,:).*out.mean(idxs(9),:); % GAEC rATP
fluxes(9,:)=fluxes(8,:)./rGlu; % GAEC YATP_glu
fluxes(10,:)=fluxes(8,:)./flux.Drate'; % GAEC YATP_mu
fluxes(11,:)=out.mean(idxs(8),:); % NGAM rATP
fluxes(12,:)=fluxes(11,:)./rGlu; % NGAM YATP_glu
fluxes(13,:)=fluxes(11,:)./flux.Drate'; % NGAM YATP_mu
fluxes(14,:)=fluxes(2,:)+fluxes(5,:)-fluxes(8,:)-fluxes(11,:); % Metabolism rATP
fluxes(15,:)=fluxes(3,:)+fluxes(6,:)-fluxes(9,:)-fluxes(12,:); % Metabolism YATP_glu
fluxes(16,:)=fluxes(4,:)+fluxes(7,:)-fluxes(10,:)-fluxes(13,:); % Metabolism YATP_mu
fluxes(17,:)=fluxes(8,:)+fluxes(11,:)+fluxes(14,:); % NGAM+GAEC+Metabolism rATP
fluxes(18,:)=fluxes(9,:)+fluxes(12,:)+fluxes(15,:); % NGAM+GAEC+Metabolism YATP_glu
fluxes(19,:)=fluxes(10,:)+fluxes(13,:)+fluxes(16,:); % NGAM+GAEC+Metabolism YATP_mu

% Primary NAD turnover reactions
fluxes(20:24,:)=out.mean(idxs(11:15),:); % rPDH, rIDH, rMDHc, rMDHm, rNDE

% Flux rates related to NAD turnover, per compartment
NADmets=find(strcmp(model.metNames,'NAD'));
[~,NADrxns]=find(model.S(NADmets,:)); clear NAD
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

% NADPH
NADPmets=find(strcmp(model.metNames,'NADP(+)'));
[~,NADPrxns]=find(model.S(NADPmets,:)); clear NADP
for i=1:numel(models)
    NADPtmp=full(models{i}.S(NADPmets,NADPrxns).*out.mean(NADPrxns,i)');
    NADPtmp(NADPtmp<0)=0;
    NADP(:,i)=sum(NADPtmp,2);
end
zeroFlux=sum(NADP,2)==0;
NADP=[NADP(~zeroFlux,:);sum(NADP,1)];
NADPmets=[strcat('rNADP[',model.comps(model.metComps(NADPmets(~zeroFlux))),']'); ...
    'rNADP[tot]']';
fluxes=[fluxes;NADP];

% Polymerization cost of biomass
fluxes=[fluxes;GAECpol];

% Define row and variable names
rowNames=[{'rGlu','ETC_rATP','ETC_YATP_glu','ETC_YATP_mu','glycolysis_rATP',...
    'glycolysis_YATP_glu','glycolysis_YATP_mu','GAEC_rATP','GAEC_YATP_glu',...
    'GAEC_YATP_mu','NGAM_rATP','NGAM_YATP_glu','NGAM_YATP_mu','Metabolism_rATP',...
    'Metabolism_YATP_glu','Metabolism_YATP_mu','GAEC+NGAM+Metabolism_rATP',...
    'GAEC+NGAM+Metabolism_YATP_glu','GAEC+NGAM+Metabolism_YATP_mu','rPDH',...
    'rIDH','rMDHc','rMDHm','rNDE'}, NADmets, NADPmets, {'GAECpol_Prot','GAECpol_Carb',...
    'GAECpol_RNA','GAECpol_DNA','GAECpol_total','GAEC_noPol','GAEC_total'}];
varNames={'CN4','CN22','CN38','CN78','hGR'};

for i=1:5
    fluxes2(:,i)=strtrim(cellstr(num2str(fluxes(:,i),4)));
end
fluxes=cell2table(fluxes2,'RowNames',rowNames,'VariableNames',varNames);   
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
