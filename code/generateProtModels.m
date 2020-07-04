% This script generates condition-specific proteome-integrated ecModels.
% Prior to running this script, run prepareEnvironment.
prepareEnvironment
for i=1:length(flux.conds)
    disp(['Construct ecModel for condition: ' flux.conds{i}])
    
    %% Preparing proteomics data
    %Extract data for the i-th condition
    abundances   = prot.data(:,repl.first(i):repl.last(i));
       
    %Calculate sample specific f-factor, before filtering data. While there
    %might be individual proteins with too much variability, this should
    %not affect f calculation too much, meanwhile ensuring higher coverage.
    cd GECKO/geckomat/limit_proteins
    f = measureAbundance(ecModel.enzymes,prot.IDs,mean(abundances,2,'omitnan'));
    sumP = sum(mean(abundances,2,'omitnan'),'omitnan');
    
    %Filter proteomics data, to only keep high quality measurements
    cd ../utilities/integrate_proteomics
    [pIDs, filtAbundances] = filter_ProtData(prot.IDs,abundances,1.96,true);
    sumP = sum(mean(filtAbundances,2,'omitnan'),'omitnan');

    disp(['Filtered out ' num2str(round((1-(numel(pIDs)/numel(prot.IDs)))*100,1)) '% of protein measurements due to low quality.'])
    cd ..
    
    %correct oxPhos complexes abundances
    for j=1:length(oxPhos)
       [filtAbundances,pIDs] = fixComplex(oxPhos{j},ecModel,filtAbundances,pIDs); end
    
    %% Scale biomass to new protein content and recalculate GAM
    %Set minimal medium
    cd ../kcat_sensitivity_analysis
    ecModelP  = changeMedia_batch(ecModel,params.c_source);
    tempModel = changeMedia_batch(ecModel_batch,params.c_source);
    cd ../limit_proteins
    %Use ecModel_batch to rescale the biomass and estimate the GAM
    [~,GAM] = scaleBioMass(tempModel,flux.Ptot(i),[],true);
    %Rescale the ecModel that will have proteomics integrated and reuse the
    %GAM calculated above
    ecModelP = scaleBioMass(ecModelP,flux.Ptot(i),GAM,true);
    disp(' ')
    
    %% Simulate condition in protein pool model to determine minimum enzyme levels
    %Block production of non-observed metabolites before data incorporation
    %and flexibilization
    expData  = [flux.GUR(i),flux.CO2prod(i),flux.OxyUptake(i)];
    flexGUR  = 1.05*flux.GUR(i); % Allow 5% deviation from measured glucose uptake
    cd ../../..
    ecModelP = DataConstrains(ecModelP,flux.byProds,flux.byP_flux(i,:),1.1);
    
    %Get a temporary model structure with the same constraints to be used
    %for minimal enzyme requirements analysis. For all measured enzymes
    %(present in the dataset) a minimal usage is obtained from a FBA
    %simulation with the ecModel_batch, then
    tempModel       = DataConstrains(tempModel,flux.byProds,flux.byP_flux(i,:),1.1);
    tempModel       = setParam(tempModel,'ub',positionsEC(1),flexGUR);
    [matchedEnz,iA] = intersect(pIDs,tempModel.enzymes);
    enzModel        = setParam(tempModel,'lb',positionsEC(2),flux.Drate(i));
    for j=1:length(matchedEnz)
        rxnIndex  = find(contains(tempModel.rxnNames,matchedEnz{j}));
        tempModel = setParam(enzModel,'obj',rxnIndex,-1);
        tempSol   = solveLP(tempModel);
        %Compare enzyme minimum usage with abundance value
        if (tempSol.x(rxnIndex)-filtAbundances(iA(j)))>0
            %Flexibilize limiting values
            disp(['Limiting abundance found for: ' matchedEnz{j} '/Previous value: ' num2str(filtAbundances(iA(j))) ' /New value: ' num2str(tempSol.x(rxnIndex))])
            filtAbundances(iA(j)) = 1.01*tempSol.x(rxnIndex);
        end
        enzIndex = find(contains(tempModel.enzymes,matchedEnz{j}));
    end
    
    %% Populate model with proteome data
    %Get model with proteomics
    cd GECKO/geckomat/limit_proteins
    disp(['Incorporation of proteomics constraints for ' flux.conds{i} ' condition'])
    sumPfilt = sum(filtAbundances);
    Ptot=flux.Ptot(i)*(sumPfilt/sumP);
    [ecModelP,usagesT,modificationsT,~,coverage] = constrainEnzymes(ecModelP,f,GAM,Ptot,pIDs,filtAbundances,flux.Drate(i),flexGUR);
    %% Fit NGAM
    cd ../utilities
    %NGAM interval for fitting
    ecModelP = setChemostatConstraints(ecModelP,positionsEC,flux.Drate(i),true,0.01);
    cd integrate_proteomics
    ecModelP = fitNGAM(ecModelP,'r_4046',expData,[0 5],true);
    cd ..
    
    %% Set all constraints to match experiment
    ecModelP = setChemostatConstraints(ecModelP,positionsEC,flux.Drate(i),true,0.01,flux.GUR(i));
    cd ../../../
    solution = solveLP(ecModelP,1);
    if ~isempty(solution.f)
        disp(['Succesful simulation of condition: ' flux.conds{i}])
    end
    eval(['ecModelP_' flux.conds{i} ' = ecModelP;']);
    save(fullfile('..','models',['ecModel_P_' flux.conds{i}]), ['ecModelP_' flux.conds{i}]);
    writetable(modificationsT,fullfile('..','results','modelGeneration',['modifiedEnzymes_' flux.conds{i} '.txt']),'Delimiter','\t')
end
