%% coreTF_analysis_sharedUnique_cellstr_heatMap_FDR1E5_listsSelTf
% listsSelTf -- some select TFs are automatically included
% get shared and unique core TFs
% Visualize enrichment of TRN with KO edges
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
% (2) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
% (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: April 24, 2018

clear all
close all
restoredefaultpath

currDir = '';

addpath(fullfile(currDir,'infLassoStARS'))
addpath(fullfile(currDir,'glmnet'))
addpath(fullfile(currDir,'customMatlabFxns'))


% plotting code
fontSize = 12;
load cbcolors.mat
load cbcolors
darkBlue = [0 0 170]./255;
mediumBlue = [0 85 255]./255;
lightBlue = [0.5843    0.8157    0.9882];
darkRed = [170 0 0]./255;
mediumRed = [228 26 28]./255;
pink = [ 0.9686    0.5059    0.7490];
blues = [darkBlue; mediumBlue; lightBlue];
reds = [darkRed; mediumRed;pink];
browns = [cbcolors(7,:)/2; cbcolors(7,:); mean([cbcolors(7,:);1 1 1])];
purples = [cbcolors(3,:)/2; cbcolors(3,:); mean([cbcolors(3,:);1 1 1])];
greens = [cbcolors(5,:)/2; cbcolors(5,:); mean([cbcolors(5,:);1 1 1])];
charcoals = [0 0 0; .6*[1 1 1]];


%% begin inputs

disp('See what happens with an overlap cutoff -- can get a sense with prior inf code.')

% GeneSetNameBase, potential Core TFs per subtype
coreInfs = {'8hpi','8hpiSet.txt';    
            '24hpi','24hpiSet.txt';    
            '48hpi','48hpiSet.txt';
            '72hpi','72hpiSet.txt';   };
geneSetDir = 'pancreas_NEUROG3_induction/inputs/geneSets';

padjMin = 1E-20;    % max for adjusted p-value

numTfs = '10';%'10'; % # of TFs per gene
netInf = 'prior_atac_Miraldi_q_ChIP_bias10_maxComb';
enrichmentDir = ['pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_' numTfs 'tfsPerGene/' netInf '/GSEA/' netInf '_cut01_8hpiSet_Praw0p1_dir_wCut0p0_minSet5'];
outDir0 = ['pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_' numTfs 'tfsPerGene/' netInf];

setName = 'Pancreas_8h_induction';
pTableSuffix = '_fdr100_SIGN_adjp.txt';
negFile = ['Pancreas_8h_induction' strrep(pTableSuffix,'SIGN','down')];
posFile = ['Pancreas_8h_induction' strrep(pTableSuffix,'SIGN','up')];

tfAnnotationInf = ''; % no annotations so far -- see figure 6 of Miraldi et al. Genome Research, where known TFs were marked, you can provide these lists here
% {'c','/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/th17/ciofani_tfList_TableS2.txt';
%     'y','/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors/mm10_merged/th17/yosef_tfList.txt'};

topNs = [Inf];%10 15 Inf];          % max # of TFs in bar graphs
FDRcutoff = .99;    % cutoff for inclusion of TF enrichment in bar graph

FDRcutoff2 = 1E-2;


totCores = size(coreInfs,1);

%% first get all p-values keeping track of all adjust p-values, then threshold below.
tmpTfNames = {};       % will keep track of TF names
% coreTfCellTypes = {};   % will also track what cell types they were significant in

tmpTfMatAct = zeros(0,totCores); % Activator core TFs
tmpTfMatRep = zeros(0,totCores); % Repressor TFs

indCoreNames = cell(totCores,1);

for cind = 1:totCores
    currCore = coreInfs{cind,1};
    coreShort = strrep(currCore,'_SI','');
    indCoreNames{cind} = coreShort;
    currPotTfs = coreInfs{cind,2};
    posPvalTables = {fullfile(enrichmentDir,posFile)};
    negPvalTables = {fullfile(enrichmentDir,posFile)};
    upGeneSet = [currCore '_up'];
    downGeneSet = [currCore '_down'];
    potCoreTfList= fullfile(geneSetDir,currPotTfs);    
    
    figure
    topN0 = Inf;
    [orderedCoreTfs, annoCoreTfs, posEnrichments, negEnrichments] ...
        = getCoreNetworks(posPvalTables,negPvalTables, upGeneSet,...
        downGeneSet, potCoreTfList, topN0, FDRcutoff, tfAnnotationInf, padjMin)        ;
    titleInf = strvcat([strrep(currCore,'_SI','') ', '  numTfs ' TFs/gene'], strrep(netInf,'_', ' '));
    title(titleInf,'FontSize',11)   
    
    disp(titleInf)
    
    % get a matrix of TFs
    [ov, coreInds, currInds] = intersect(tmpTfNames,orderedCoreTfs);
    if length(ov)
        tmpTfMatAct(coreInds,cind) = posEnrichments(currInds);
        tmpTfMatRep(coreInds,cind) = negEnrichments(currInds);
%         for tf = 1:length(coreInds)
%             coreTfCellTypes{coreInds(tf)} = [coreTfCellTypes{coreInds(tf)} ',' coreShort];
%         end
    end
    [newTfs, currInds] = setdiff(orderedCoreTfs,tmpTfNames);
    if length(newTfs)
        tmpTfNames = cellstr(strvcat(strvcat(tmpTfNames),strvcat(newTfs)));
        totNew = length(newTfs);        
%         coreTfCellTypes = cellstr(strvcat(strvcat(coreTfCellTypes),repmat(coreShort,totNew,1)));        
        newMat = zeros(totNew,totCores);
        newMat(:,cind) = posEnrichments(currInds);
        tmpTfMatAct = [tmpTfMatAct; newMat];
        newMat = zeros(totNew,totCores);
        newMat(:,cind) = negEnrichments(currInds);
        tmpTfMatRep = [tmpTfMatRep; newMat];
    end
        
%     strvcat(orderedCoreTfs)        
    % save figure
%     figOutBase = fullfile(outDir,[strrep(strrep(currCore,'+','p'),'-','m') '_' numTfs 'TFsPG']);
%     if figOutBase    
%         saveas(gcf,figOutBase,'fig')
%         fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [7 5.25]);
%         print('-painters','-dpdf','-r150',[figOutBase '.pdf'])
%         disp([figOutBase '.fig + .pdf created.'])
%     end

%     % output ordered TFs
%     currFig = [figOutBase '_notMerged.txt'];
%     fout = fopen(currFig,'w');
%     fprintf(fout,strjoin(upper(orderedCoreTfs),'\n'));
%     fclose(fout);
%     disp(['TFs out: ' currFig])
    
%     pause

end



%% now threshold and get core TFs and their enrichments
enrichmentsSum = tmpTfMatAct + tmpTfMatRep;
figure(1)
subplot(1,3,1)
imagesc(tmpTfMatAct)
set(gca,'YTick',1:size(tmpTfMatAct,1),'YTickLabel',tmpTfNames)
set(gca,'CLim',[-10 10])
set(gca,'XTick',1:totCores,'XTickLabel',indCoreNames,'XTickLabelRotation',90)
subplot(1,3,2)
imagesc(tmpTfMatRep)
set(gca,'YTick',1:size(tmpTfMatAct,1),'YTickLabel',tmpTfNames)
set(gca,'XTick',1:totCores,'XTickLabel',indCoreNames,'XTickLabelRotation',90)
set(gca,'CLim',[-10 10])
subplot(1,3,3)
imagesc(enrichmentsSum)
set(gca,'YTick',1:size(tmpTfMatAct,1),'YTickLabel',tmpTfNames)
set(gca,'XTick',1:totCores,'XTickLabel',indCoreNames,'XTickLabelRotation',90)
set(gca,'CLim',[-10 10])
colormap redblue

% gompers
%% just so happens that, for these networks nothing is both an activator and repressor in a cell type
maxes = max(abs(enrichmentsSum),[],2);
keepTfs = find(maxes > -log10(FDRcutoff2));

coreTfs = '';
coreTfCellTypes = '';
coreTfMat = [];
for tind = 1:length(tmpTfNames);
    currEnrich = enrichmentsSum(tind,:);
    testEnrich = max(abs([tmpTfMatAct(tind,:); tmpTfMatRep(tind,:)]));
    cellInds = find(abs(testEnrich) > -log10(FDRcutoff2));
    if cellInds
        coreTfs = strvcat(coreTfs,tmpTfNames{tind});
%         if length(cellInds) > 1
            coreTfCellTypes = strvcat(coreTfCellTypes,strjoin({indCoreNames{cellInds}},','));
%         else
%             coreTfCellTypes = strvcat(coreTfCellTypes,indCoreNames{cellInds});
%         end
        coreTfMat = [coreTfMat; currEnrich];
    end
end

coreTfCellTypes = cellstr(coreTfCellTypes);
coreTfs = cellstr(coreTfs);
totCoreTfs = length(coreTfs);

uniCores = unique(coreTfCellTypes);
totUniCores = length(uniCores);

tfIndsOrdered = zeros(totCoreTfs,1);
clusterSpots = zeros(totCoreTfs,1);
last = 0;


for topN = topNs    
    outBase = [setName '_2fdr' strrep(num2str(FDRcutoff2,'%1.1E'),'.','p') '_top' num2str(topN)]; 
    outDir = fullfile([outDir0 '/' outBase]);
    mkdir(outDir)
    alphDir = fullfile(outDir,'alph_upper');
    mkdir(alphDir)
end


for uind = 1:totUniCores
    currCore = uniCores{uind};
    currInds = find(ismember(coreTfCellTypes,currCore));
    coreInds = find(ismember(indCoreNames,strsplit(currCore,',')));
    currEnriches = abs(coreTfMat(currInds,coreInds));        
    % combine p-values using z-method
    zs = norminv(2*(10.^(-currEnriches)));
    meanZ = mean(zs,2);
    combEns = -log10(normcdf(meanZ)/2);
    [vals, inds] = sort(combEns,'descend');

    for topN = topNs    
        outBase = [setName '_2fdr' strrep(num2str(FDRcutoff2,'%1.1E'),'.','p') '_top' num2str(topN)]; 
        outDir = fullfile([outDir0 '/' outBase]);
        
        keepInds = currInds(inds(find(vals >= vals(min(topN,length(currInds)))))); % consider ties, and topN = Inf
        if length(keepInds) > 0
            tfsInCore = {coreTfs{keepInds}}';
            figOutBase = fullfile(outDir,[strrep(strrep(strrep(currCore,'+','p'),'-','m'),',','_') ]);
            figText = [figOutBase '.txt'];
            fout = fopen(figText,'w');
            fprintf(fout,strjoin(tfsInCore,'\n'));
            fclose(fout);
            disp(figText)
            
            
            figOutBase = fullfile(alphDir,[strrep(strrep(strrep(currCore,'+','p'),'-','m'),',','_')]);
            figText = [figOutBase '.txt'];
            fout = fopen(figText,'w');
            fprintf(fout,strjoin(sort(upper(tfsInCore)),', '));
            fclose(fout);
            disp(figText)

            
            
        end
    end
end
