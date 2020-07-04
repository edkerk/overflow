function model = DataConstrains(model,compounds,bounds,flexBounds)
if ~isempty(compounds)
    disp('Constraining byproducts exchange fluxes with fermentation data')
    for i=1:length(compounds)
        %Get exchange rxn index
        if ~strcmpi(compounds{i},'oxygen')
            rxnName = [compounds{i} ' exchange'];
        else
            rxnName = [compounds{i} ' exchange (reversible)'];
        end
        BPindex = find(strcmpi(model.rxnNames,rxnName));
        if ~isempty(BPindex)
            disp([compounds{i} ' exchange has been constrained to: ' num2str(bounds(i)) ' [mmol/gDw h]'])
            %Allow some flexibility
            model = setParam(model,'ub',BPindex,flexBounds(1)*bounds(i));
            if numel(flexBounds)>1
                model = setParam(model,'lb',BPindex,flexBounds(2)*bounds(i));
            end
        else
            disp(['No exchange rxn for ' compounds{i} ' was found in ecModel'])
        end
    end
end
disp(' ')
end