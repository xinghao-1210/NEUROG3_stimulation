%% vis_th17_coreAnalysis
%% Goal: visualize inferred "core" regulators for Th17 cells, based on GSEA 
% of TF target genes (e.g., with tfTargets_GSEA.sh), where p-values of 
% enrichment were calculated. 
% To determine the core, we are interested in TFs that promote
% Th17 gene expression patterns. Specifically, we combine p-values of
% enrichment of (1) postive edges in Th17-up-regulated and (2) negative
% edges in Th17-down-regulated genes.
%% Reference: Miraldi et al. (2019) Genome Research.
%% Visualize most Top N most significant Th17 core regulators in a bar graph
%% We add metadata to gene symbol name IF the gene was a Th17
% regulator from Yosef et al. or Ciofani et al. 

clear all
close all

currDir = '';
addpath(fullfile(currDir,'customMatlabFxns'))

%% max-combined Th17 TF mRNA & TFA TRN
cond_time=72
outFileBase = sprintf('pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/GSEA/prior_atac_Miraldi_q_ChIP_bias10_maxComb_cut01_%dhpiSet_Praw0p1_dir_wCut0p0_minSet5',cond_time);
titleBit = sprintf('max-combined pancreas_NEUROG3_8h %dhpi TF mRNA & TFA TRN',cond_time);
outDir0 = sprintf('pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/%dhpi_Cores',cond_time);

potRegList = 'pancreas_NEUROG3_induction/inputs/targRegLists/potRegs_names.txt';
topN = 30;

padjMin = 1E-120; % max for adjusted p-value
setInfIn = sprintf('%dhpiset_fdr100',cond_time);   % for input set
FDR_cutoff = .01;   % cutoff for inclusion of enrichment in heatmap
setInf = [sprintf('%dhpiSets_fdr',cond_time) num2str(100*FDR_cutoff) '_top' num2str(topN)];     % for output set

prevTh17tfInfs = {'p','pancreas_NEUROG3_induction/inputs/geneSets/prevgenes.txt'};

outDir = fullfile([outDir0 '/' setInf]);
mkdir(outDir)

makeTextOut = 1; % output TF names as a text file

distance_metric = '';
data_scaling = '';
fc_cutoff = 0;
fontSize = 8;
lineWidth = .5;

th17promActive = {sprintf('%dhpi up',cond_time)};       % strip "_" (line100) for later intersect (line103), either not include "_" here or "setsTmp = C{1}" (line100)
th17promRepress = {sprintf('%dhpi down',cond_time)};    % strip "_" (line100) for later intersect (line103), either not include "_" here or "setsTmp = C{1}" (line100)

%% END parameters

titleInf = strrep([titleBit ' ' setInf],'_',' ');
disp(titleInf)

% get list of potential regulators
fin = fopen(potRegList,'r');
C = textscan(fin, '%s','HeaderLines',0);
fclose(fin);
potRegs = C{1};
totPotRegs = length(potRegs);

pvalStruct = [];
pvalStruct(1).type = 'up';
pvalStruct(2).type = 'down';
pvalStruct(1).factor = 1;
pvalStruct(2).factor = -1;
pvalStruct(1).Th17promGroups = th17promActive;    
pvalStruct(2).Th17promGroups = th17promRepress;    
% set up p-value structure for info
for sind = 1:2
    pvalStruct(sind).pvals = ones(length(pvalStruct(1).Th17promGroups),totPotRegs);
    % set to zero initially as z-scores
    pvalStruct(sind).tfs = potRegs;
end

%% load regulators from previous networks
totPrevs = size(prevTh17tfInfs,1);
prevSymbs = {prevTh17tfInfs{:,1}};
prevRegs = cell(totPrevs,1); % for storing TF names
for pr = 1:totPrevs
    currFile = prevTh17tfInfs{pr,2};
    fin = fopen(currFile,'r');
    C = textscan(fin,'%s','HeaderLines',0);
    fclose(fin);
    prevRegs{pr} = C{1};
end
    
%% load adjusted pvals from up and down
totTypes = length(pvalStruct);

for tind = 1:totTypes
    pvalFile = fullfile(outFileBase,[setInfIn '_' pvalStruct(tind).type '_adjp.txt']);
    fid = fopen(pvalFile);  % get first line, TFs
    tline=fgetl(fid);
    tfTmps = strvcat(strsplit(tline,'\t'));
    fclose(fid);
    tfsTmpc = cellstr(tfTmps);
    totTfs = length(tfsTmpc);
    fid = fopen(pvalFile);  % get pvals + sets
    C = textscan(fid,['%s' repmat('%f',1,totTfs)],'Delimiter','\t','Headerlines',1);
    setsTmp = strrep(C{1},'_',' '); % strip "_" (line100) for later intersect (line103), either not include "_" (line46/47) or "setsTmp = C{1}" here
    padjsTmp = [C{2:end}];
    %% 1. Get All p-values for sets of interests
    [sets,outInds,inInds] = intersect(pvalStruct(tind).Th17promGroups,setsTmp);
    [comTfs,pvalOrd,inOrd] = intersect(potRegs,tfTmps);
    pvalStruct(tind).pvals(outInds,pvalOrd) = padjsTmp(inInds,inOrd);
end

promTh17Tfs = {};
for tind = 1:totTypes
%% 2. Find Th17 promotor TFs
pmins = min(pvalStruct(tind).pvals,[],1);
[keep] = find(pmins<=FDR_cutoff);
promTh17Tfs = union(promTh17Tfs,{potRegs{keep}}');
end

%% combine the matrices including things significant in either
totTfs = length(promTh17Tfs);
tfInds = find(ismember(potRegs,promTh17Tfs));

allSets = '';
allEnriches = [];
for tind = 1:totTypes
    allEnriches = [allEnriches; pvalStruct(tind).factor * -log10(max(pvalStruct(tind).pvals(:,tfInds),padjMin))];
    allSets = strvcat(allSets,strvcat(pvalStruct(tind).Th17promGroups));    
    
end
allSets = cellstr(allSets);

% augment promTh17Tfs names to add metadata
totPrevs = size(prevTh17tfInfs,1);
prevSymbs = {prevTh17tfInfs{:,1}};
prevRegs = cell(totPrevs,1); % for storing TF names
annoTfs = cell(totTfs,1);
for tt = 1:totTfs
    currTf = promTh17Tfs{tt};
    baseExtra = '^{';
    prevOnce = 0;
    for pt = 1:totPrevs
        isPrev = length(find(ismember(prevRegs{pt},currTf)));
        if isPrev
            baseExtra = [baseExtra prevSymbs{pt} ','];
            prevOnce = 1;
        end
    end
    if prevOnce % remove last comma and add }'
        baseExtra = [baseExtra(1:end-1) '}'];
    else
        baseExtra = '';
    end
    annoTfs{tt} = [upper(currTf) baseExtra];
end

%% order based on significance
maxes = max(abs(allEnriches),[],1);
[vals, horderTfs] = sort(maxes,'ascend');

topN = min(totTfs,topN); % take care of infinity case
horderTfs = horderTfs(end-topN+1:end);
currTfs = {annoTfs{horderTfs}}';
[alphTfs, alphOrder] = sort(currTfs);
alphTfs = flipud(alphTfs);
alphOrder = flipud(alphOrder);

%% alphabetize TFs
outBit = sprintf('_%dhpipromTFs_padjsORD_alph',cond_time);
currFig = [outDir outBit];
if makeTextOut
    fout = fopen([currFig '.txt'],'w');
    fprintf(fout,strjoin({promTh17Tfs{horderTfs(alphOrder)}},'\n'));
    fclose(fout);
    disp(['TFs out: ' currFig '.txt'])
end

%% visualize enrichment
figure(2), clf
subplot(1,3,2:3)

disp('hard-coded bar graph plotting of enrichments')
barh(allEnriches(1,horderTfs(alphOrder)),'FaceColor','r')
hold on
barh(allEnriches(2,horderTfs(alphOrder)),'FaceColor','b')
set(gca,'YTick',1:topN,'YTickLabel',alphTfs,...
    'FontSize',fontSize,'FontWeight','Bold')
xlabel('(Edge Sign) * -log_{10}(P_{adj})','FontSize',fontSize)
title(titleInf,'FontSize',fontSize + 2)
grid on, grid minor

ax = axis();
axis([ax(1:2) 0 topN+1])

currFig = [outDir outBit];
saveas(gcf,currFig,'fig')
pause
save2pdf(currFig,gcf,300)
disp('Completed')
disp(currFig)

