%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,count,sum] = measureAbundance(enzymes)
% 
% Benjamin J. Sanchez. Last edited: 2018-08-10
% Ivan Domenzain.      Last edited: 2018-11-26

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,count] = measureAbundance(enzymes,pIDs,abundance)
if nargin==1 % If 
    genes     = {};
    abundance = [];
    fileName = '../../databases/prot_abundance.txt';
    if exist(fileName,'file')~= 0
        fID         = fopen(fileName);
        data        = textscan(fID,'%s','delimiter','\n');
        headerLines = sum(startsWith(data{1},'#'));
        fclose(fID);
        %Read data file, excluding headerlines
        fID         = fopen('../../databases/prot_abundance.txt');
        data        = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',headerLines);
        genes       = data{2};
        %Remove internal geneIDs modifiers
        genes     = regexprep(genes,'(\d{4}).','');
        abundance = data{3};
        fclose(fID);
    else
        fileName = '../../databases/relative_proteomics.txt';
        if exist(fileName,'file')~= 0
            fID       = fopen(fileName);
            data      = textscan(fID,'%s %f','delimiter','\t','HeaderLines',1);
            genes     = data{1};
            abundance = data{2};
            fclose(fID);
        end
    end
end

if ~isempty(abundance)
    %Load swissprot data:
    data      = load('../../databases/ProtDatabase.mat');
    swissprot = data.swissprot;
    for i = 1:length(swissprot)
        swissprot{i,3} = strsplit(swissprot{i,3},' ');
    end
    %Main loop:
    MW_ave  = mean(cell2mat(swissprot(:,5)));
    concs   = zeros(size(abundance));
    counter = false(size(abundance));
    for i = 1:length(abundance)
        swIdx=0;
        MW = MW_ave;
        if exist('pIDs','var')
            swIdx = find(contains(swissprot(:,1),pIDs{i}));
        else
            %Find gene in swissprot database:
            for j = 1:length(swissprot)
                if sum(strcmp(swissprot{j,3},genes{i})) > 0
                    swIdx = j;
                end
            end
        end
        if swIdx~=0
            MW = swissprot{swIdx,5};	%g/mol
            %Check if uniprot is in model:
            if contains(swissprot(swIdx,1),enzymes)
                counter(i) = true;
            end
        end
        concs(i) = MW*abundance(i);     %g/mol(tot prot)
        if rem(i,100) == 0 || i == length(abundance)
            disp(['Calculating total abundance: Ready with ' num2str(i) '/' ...
                  num2str(length(abundance)) ' genes '])
        end
    end
    f     = sum(concs(counter),'omitnan')/sum(concs,'omitnan');
    count = [length(counter);sum(counter,'omitnan')];
else
    disp('Abundances not provided or files could not be found. A default value of f=0.5 is set instead')
    f     = 0.5;
    count = 0;
end
end

