% This script takes the ecmodels reconstructed by generateProtModels
% and adds ribosome subunits to a translation pseudoreaction.

prepareEnvironment
load('../models/ecModel_P_CN4.mat')
load('../models/ecModel_P_CN22.mat')
load('../models/ecModel_P_CN38.mat')
load('../models/ecModel_P_CN78.mat')
load('../models/ecModel_P_hGR.mat')

clear oxPhos params positionsEC grouping ecModel ecModel_batch

ecModels{1}=ecModelP_CN4;
ecModels{2}=ecModelP_CN22;
ecModels{3}=ecModelP_CN38;
ecModels{4}=ecModelP_CN78;
ecModels{5}=ecModelP_hGR;

%% Load ribosome information
fid=fopen('../data/ribosome.txt');
data=textscan(fid,'%q %q %q %q %q','HeaderLines',1,'Delimiter','\t');
fclose(fid);
ribo=data{1};

%Keep only the most abundant subunite (>1e5, measured in all samples
[~,repRibo.dataIdx]=ismember(ribo,prot.IDs);
repRibo.protein=ribo(repRibo.dataIdx>0);
repRibo.dataIdx(repRibo.dataIdx==0)=[];
repRibo.avgLevel=mean(prot.data(repRibo.dataIdx,:),2);
%Density plot to see bi-modal distribution of subunit abundances
[f,xi]=ksdensity(log10(repRibo.avgLevel),'Bandwidth',0.1);
plot(xi,f);
xlabel('Subunit abundance (log10(mmol/gDCW))');
ylabel('Density');
title('Distribution of average ribosomal subunit abundances');
saveas(gca,fullfile('..','results','modelGeneration','average_riboSubunit_abundance.pdf'));
%Only include ribosomal subunite with abundance over 1e-5 mmol/gDCW 
rmRibo=repRibo.avgLevel<1e-5 | isnan(repRibo.avgLevel);
ribo=repRibo.protein(~rmRibo);

%Prepare subunit information to be added to each model
[~,idx]=ismember(ribo,data{1});
enzNames=data{2}(idx);
enzGenes=data{3}(idx);
MWs=str2double(data{4}(idx))/1000;
sequences=data{5}(idx);
pathways={'sce03010  Ribosome'};

clear repRibo rmRibo idx data fid ans
%% Modify reactions
adjusted=cell.empty();
for j=1:numel(ecModels);
    disp(['Add ribosome subunits to condition: ' flux.conds{j}])
    model=ecModels{j};
    %Add enzymes and genes
    model.enzymes=[model.enzymes;ribo];
    model.enzNames=[model.enzNames;enzNames];
    model.enzGenes=[model.enzGenes;enzGenes];
    genesToAdd.genes=enzGenes;
    genesToAdd.geneShortNames=enzNames;
    model=addGenesRaven(model,genesToAdd);
    model.MWs=[model.MWs;MWs];
    model.sequences=[model.sequences;sequences];
    model.pathways=[model.pathways;repmat(pathways,numel(ribo),1)];

    %Include new amino acid pseudometabolite
    metsToAdd.mets={'aminoAcids'};
    metsToAdd.metNames={'amino acids for protein'};
    metsToAdd.compartments={'c'};
    model=addMets(model,metsToAdd);
    
    %Modify protein pseudoreaction to produce amino acid pseudometabolite
    protMetIdx=getIndexes(model,'protein','metnames');
    aaMetIdx=getIndexes(model,'aminoAcids','mets');
    protRxnIdx=getIndexes(model,'r_4047','rxns');
    model.S([protMetIdx,aaMetIdx],protRxnIdx)=[0,1];
    
    %Include ribosome subunits as pseudometabolites
    riboToAdd.mets=strcat('prot_',sort(ribo));
    riboToAdd.metNames=riboToAdd.mets;
    riboToAdd.compartments='c';
    model=addMets(model,riboToAdd);
    
    %Determine stoichiometric coefficient of ribosomes
    mmolAA=full(model.S(:,protRxnIdx));
    mmolAA=-sum(mmolAA(mmolAA<0)); %mmol amino acids in biomass
    riboKcat=10.5*3600; %10.5 aa/s -> p hour
    riboKcat=mmolAA/riboKcat; %compensate for the amount of amino acids elongated
    % Include new reaction representing ribosomes (=translation)
    rxnsToAdd.rxns={'translation'};
    rxnsToAdd.mets=[model.mets(protMetIdx),model.mets(aaMetIdx),riboToAdd.mets'];
    rxnsToAdd.stoichCoeffs=[1,-1,repmat(-riboKcat,1,numel(riboToAdd.mets))];
    rxnsToAdd.subSystem={'sce03010  Ribosome'};
    rxnsToAdd.grRules={strjoin(enzGenes,' and ')};
    model=addRxns(model,rxnsToAdd);
    
    %Add exchange reactions or usage reactions
    riboExchId=numel(model.rxns)+1;
    model=addExchangeRxns(model,'in',riboToAdd.mets);
    riboExchId(2)=numel(model.rxns);
    %Add UB for enzyme exchange reactions based on measurements.
    %If UB is too low, then adjust to value predicted by model.
    cd GECKO/geckomat/utilities/integrate_proteomics
    abundances   = prot.data(:,repl.first(j):repl.last(j));
    [pIDs, filtAbundances] = filter_ProtData(prot.IDs,abundances,1.96,true);
    cd ../../../..
    sol=solveLP(model);
    for i=riboExchId(1):riboExchId(2)
        protId=regexprep(model.rxnNames{i},'prot_(......).*','$1');
        k=find(strcmp(protId,pIDs));
        if ~isempty(k)
            model.rxns{i}=['prot_' protId '_exchange'];
            model.rxnNames{i}=['prot_' protId '_exchange'];
            model.concs(strcmp(protId,model.enzymes))=filtAbundances(k);
            if sol.x(i)<filtAbundances(k)
                model.ub(i)=filtAbundances(k);
            else
                model.ub(i)=sol.x(i)*1.01;
                adjusted{end+1,1}=protId;
                adjusted{end,2}=filtAbundances(k);
                adjusted{end,3}=sol.x(i)*1.01;
                adjusted{end,4}=flux.conds{j};
                fprintf('%s abundance adjusted. Measured: %e / Adjusted: %e\n',adjusted{end,1:3})
            end
        else
            model.ub(i)=Inf;
            model.rxns{i}=['draw_prot_' protId];
            model.rxnNames{i}=['draw_prot_' protId];
        end
    end
    ecModels{j}=model;
end
adjusted=adjusted';
fid=fopen(fullfile('..','results','modelGeneration','modifiedRibosomeSubunits.txt'),'w');
fprintf(fid,'protein_IDs previous_values modified_values condition\n');
fprintf(fid,'%s %f %f %s\n',adjusted{:});
fclose(fid);

clear cond abundances aaMetIdx filtAbundances i j genesToAdd metsToAdd mmolAA
clear protId prot pIDs pathways MWs model enzGenes enzNames enzAdjust protMetIdx
clear protRxnIdx repl sequences riboExchId riboKcat riboToAdd rxnsToAdd
clear sequences sol fid adjusted