function [capUsage,absUsage,protName] = enzymeUsage(ecModel,fluxes,zero)
% enzymeUsage
%
%   Calculate enzyme usage
%
%   Input:
%   ecModel         enzyme-constrained model
%   fluxes          vector of fluxes, for instance sol.x from sol=solveLP
%   zero            logical whether also zero enzyme usages should be
%                   included (opt, default true)
%
%   Output:
%   capUsage        ratio of enzyme usage per available enzyme, i.e.
%                   capacity usage (flux divided by ub)
%   absUsage        absolute enzyme usage (flux)
%
% Note: the order of the enzymes in capUsage and absUsage is identical to
% the ecModel.enzGenes and ecModel.enzymes that are used as inputs.
%
% Usage: [capUsage,absUsage] = enzymeUsage(ecModel,fluxes,zero)
%
% Eduard Kerkhoven  Last edited: 2019-01-25

if nargin<3
    zero=true;
end

if isfield(ecModel,'enzymes')
    protIdx     = find(contains(ecModel.rxnNames,ecModel.enzymes));
    matchProt   = regexprep(ecModel.rxnNames(protIdx),'(draw_)?prot_','');
    matchProt   = regexprep(matchProt,'_exchange.*','');
    [~,b]       = ismember(ecModel.enzymes, matchProt);

    absUsage    = fluxes(protIdx);
    capUsage    = absUsage./ecModel.ub(protIdx);

    absUsage    = absUsage(b);
    capUsage    = capUsage(b);
    protName    = ecModel.enzymes;
else
    expression  = '^(draw_prot_.*)|(prot_.*_exchange$)';
    protRxn     = regexp(ecModel.rxnNames,expression);
    protRxn     = find(~cellfun(@isempty,protRxn));
    expression  = '(^draw_prot_)|(^prot_)|(_exchange$)';
    protName    = regexprep(ecModel.rxnNames(protRxn),expression,'');
    
    absUsage    = fluxes(protRxn);
    capUsage    = absUsage./ecModel.ub(protRxn);
end
    
    
if ~zero
    nonzero     = absUsage>0;
    absUsage    = absUsage(nonzero);
    capUsage    = capUsage(nonzero);
    protName    = protName(nonzero);
end
end
