%% aveGeneExpMatrix_subset_4viz
%% Subset and average particular conditions (or cells) from the gene
% expression data matrix for visualization. 
% * User specifies new names for conditions and experiments to average
% * User specifies the order of select conditions and identifies clusters
% of conditions (e.g., for visualization in heatmaps in vis_allTfTfOverlaps.m) 
%% OUTPUTS:
% 1. .txt file of gene expression values, with averaged conditions
%   ordered as specified. (e.g., can serve as input to jp_gene_viz)
% 2. .mat file of gene expresison values with cluster info (e.g., for 
%   visualization in heatmaps in vis_allTfTfOverlaps.m) 
%% NOTE: The output gene expression matrix is limited to TFA genes.

clear all
close all

%% INPUTS
matlabDir = '';
addpath(fullfile(matlabDir,'customMatlabFxns'))

geneExprTFAdir = 'pancreas_NEUROG3_induction/outputs/processedGeneExpTFA';
outDir = fullfile(geneExprTFAdir,'geneExpHeatmapInputs');
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

%% condGroups -- each row is 3D
% Col 0 = a field not used in this code
% Col 1 = name for group of conditions -- will form output file name
% Col 2 = table with three columns (1 --> add a group line on heatmap
%    2 --> a name for set of samples, 3--> cell containing individual
%    sample names (as found in conditionsc)
condGroups = {...
    {'pancreas_8h_8-72hpi';...
        {1,'100_8h(8hpi)',{'wt_100_8h_8hpi_1';'wt_100_8h_8hpi_2';'wt_100_8h_8hpi_3'};
        1,'100_8h(24hpi)',{'wt_100_8h_24hpi_1';'wt_100_8h_24hpi_2';'wt_100_8h_24hpi_3';};
        1,'100_8h(48hpi)',{'wt_100_8h_48hpi_1';'wt_100_8h_48hpi_2';'wt_100_8h_48hpi_3';};
        1,'100_8h(72hpi)',{'wt_100_8h_72hpi_1';'wt_100_8h_72hpi_2';'wt_100_8h_72hpi_3';};
        1,'0_8h(8hpi)',{'wt_0_8h_8hpi_1';'wt_0_8h_8hpi_2';'wt_0_8h_8hpi_3'};
        1,'0_8h(24hpi)',{'wt_0_8h_24hpi_1';'wt_0_8h_24hpi_2';'wt_0_8h_24hpi_3';};
        1,'0_8h(48hpi)',{'wt_0_8h_48hpi_1';'wt_0_8h_48hpi_2';'wt_0_8h_48hpi_3';};
        1,'0_8h(72hpi)',{'wt_0_8h_72hpi_1';'wt_0_8h_72hpi_2';'wt_0_8h_72hpi_3';};}};
    };

%% END INPUTS

mkdir(outDir)
disp(outDir)

totGroups = length(condGroups);
load(geneExprMat)

for groupInd = 1:totGroups
    condGroup = condGroups{groupInd};
    condGroupName = condGroup{1}; condInf = condGroup{2};
    startSpotsConds = [condInf{:,1}];
    condPrintNames = {condInf{:,2}}';
    indConds = {condInf{:,3}}';

    outBase = fullfile(outDir,[condGroupName]);

    %% Limit gene expression matrix to specific conditions
    aveCounts = [];      % averaged conditions (corresponding to condPrintNames)
    totBigConds = length(condPrintNames);
    aveGeneExprVals = [];
    for cind = 1:totBigConds
        currSamps = indConds{cind}';
        totIndConds = length(currSamps);    
        sampInds = zeros(totIndConds,1);
        for icind = 1:totIndConds
            sampInds(icind) = find(ismember(conditionsc,currSamps{icind}));
        end
        currVals = tfaGeneMat(:,sampInds);
        aveCounts = [aveCounts median(currVals,2)];
    end

    %% generate output text file
    outExprFile = [outBase '.txt'];
    fout = fopen(outExprFile,'w');
    fprintf(fout,['\t' strjoin(condPrintNames,'\t') '\n']);
    for gene = 1:length(tfaGenes)
        fprintf(fout,[tfaGenes{gene} '\t' ...
            strjoin(cellstr(num2str(aveCounts(gene,:)')),'\t') '\n']);
    end
    disp([outExprFile ' generated.'])
    disp([num2str(length(tfaGenes)) ' TFA genes included in output matrix.'])
    
    %% generate output .mat file
    genesc = tfaGenes;
    totConds = length(condPrintNames);
    zAveCounts = zscore(aveCounts')';
    meanC = mean(2.^aveCounts,2);
    l2AveCounts = log2((2.^aveCounts)./repmat(meanC,1,totConds));


    matOut = [outBase '.mat'];
    save(matOut,...
        'zAveCounts',...
        'l2AveCounts',...
        'meanC',...
        'aveCounts',...
...        'conds2average',...
        'condPrintNames',...
        'startSpotsConds',...
        'genesc')
    disp([matOut ' created.'])

end
