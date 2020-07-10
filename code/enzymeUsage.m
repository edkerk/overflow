function [absUsage,capUsage,UB,protId] = enzymeUsage(ecModel,fluxes,zero)
% enzymeUsage
%
%   Gives enzyme usages based on a provided flux distribution. It can give:
%   1)  absolute usage: the specific enzyme usage in mmol/gDCW/h, which can
%       be given for enzymes both with- and without abundance information;
%   2)  capacity usage: the ratio of available enzyme that is used, calcuted
%       by (absUsage/UB) (note that capacity usage is 0 if an enzyme
%       abundance was not constrained in the model);
%   3)  UB: the upper bound of each enzyme exchange reaction;
%   4)  protId: the protein identifiers for each enzyme (if the model has an
%       enzymes field than this order is used, otherwise it is given
%       alphabetically.
%
%   Input:
%   ecModel         enzyme-constrained model
%   fluxes          vector of fluxes, for instance sol.x
%   zero            logical whether also zero absolute enzyme usages should be
%                   included (opt, default true)
%
%   Output:
%   capUsage        vector of enzyme capacity usages
%   absUsage        vector of absolute enzyme usages
%   UB              vector of enzyme exchange reaction upper bounds
%   protId          string array of protein IDs matching the other output
%
% Usage: [absUsage,capUsage,UB,protId] = enzymeUsage(ecModel,fluxes,zero)

if nargin<3
    zero=true;
end

if isfield(ecModel,'enzymes')
    protIdx     = find(contains(ecModel.rxnNames,ecModel.enzymes));
    matchProt   = regexprep(ecModel.rxnNames(protIdx),'(draw_)?prot_','');
    matchProt   = regexprep(matchProt,'_exchange.*','');
    [~,b]       = ismember(ecModel.enzymes, matchProt);

    absUsage    = fluxes(protIdx);
    UB          = ecModel.ub(protIdx);
    capUsage    = absUsage./UB;
   
    absUsage    = absUsage(b);
    capUsage    = capUsage(b);
    UB          = UB(b);
    protId      = ecModel.enzymes;
else
    expression  = '^(draw_prot_.*)|(prot_.*_exchange$)';
    protRxn     = regexp(ecModel.rxnNames,expression);
    protRxn     = find(~cellfun(@isempty,protRxn));
    expression  = '(^draw_prot_)|(^prot_)|(_exchange$)';
    protId      = regexprep(ecModel.rxnNames(protRxn),expression,'');
    
    absUsage    = fluxes(protRxn);
    UB          = ecModel.ub(protRxn);
    capUsage    = absUsage./UB;
end

if true(iszero)
    nonzero     = absUsage>0;
    absUsage    = absUsage(nonzero);
    capUsage    = capUsage(nonzero);
    UB          = UB(nonzero);
    protId      = protId(nonzero);
end
end
