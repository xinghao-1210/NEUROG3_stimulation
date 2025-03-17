%% example_workflow_Th17_r2Pred
% use out-of-sample prediction to select model size (e.g., average # of 
% TFs / gene) using mLASSO-StARS to build a TRN
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
%% Date: May 15, 2018

clear all
close all
restoredefaultpath

currDir = '';

addpath(fullfile(currDir,'infLassoStARS'))
addpath(fullfile(currDir,'glmnet'))
addpath(fullfile(currDir,'customMatlabFxns'))

geneExprTFAdir = 'HIOE_NEUROG3_induction/outputs/processedGeneExpTFA';

%% 1. If necessary, import gene expression data, list of regulators, list 
% of target genes into a Matlab .mat object
normGeneExprFile = 'HIOE_NEUROG3_induction/inputs/geneExpression/RNAseq_24_DESeq2_VSDcounts.txt';
targGeneFile = 'HIOE_NEUROG3_induction/inputs/targRegLists/targetGenes_names.txt';
potRegFile = 'HIOE_NEUROG3_induction/inputs/targRegLists/potRegs_names.txt';
tfaGeneFile = 'HIOE_NEUROG3_induction/inputs/targRegLists/genesForTFA.txt';
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

try
    ls(geneExprMat)
catch    
    disp('1. importGeneExpGeneLists.m')
    importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
        tfaGeneFile,geneExprMat)
end   

%% 2. Given a prior of TF-gene interactions, estimate transcription factor 
% activities (TFAs) using prior-based TFA and TF mRNA levels
priorName = 'prior_atac_Miraldi_q_ChIP';
priorFile = ['HIOE_NEUROG3_induction/inputs/priors/' priorName '.tsv']; % Th17 ATAC-seq prior
edgeSS = 50;
minTargets = 3;
[xx, priorName, ext] = fileparts(priorFile);
tfaMat = fullfile(geneExprTFAdir,[priorName '_ss' num2str(edgeSS) '.mat']);

try
    ls(tfaMat)
catch
    disp('2. integratePrior_estTFA.m')
    integratePrior_estTFA(geneExprMat,priorFile,edgeSS,...
         minTargets, tfaMat)
end

%% Evaluate models for a variety of leave-out gene expression datasets, 
%% lambda bias values, and both TFA methods

loInfo = {'HIOE_NEUROG3_induction/inputs/leaveOutLists/0-24hpi.txt','_0-24hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/0-48hpi.txt','_0-48hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/0-72hpi.txt','_0-72hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/0-96hpi.txt','_0-96hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/100-24hpi.txt','_100-24hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/100-48hpi.txt','_100-48hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/100-72hpi.txt','_100-72hpi';
    'HIOE_NEUROG3_induction/inputs/leaveOutLists/100-96hpi.txt','_100-96hpi';};

lambdaBiases = [1 .5 .25 .1]; % correspond to "no,moderate, and strong prior reinforcement
tfaOpts = {'','_TFmRNA'}; % the two TFA options

totLos = size(loInfo,1);
totTfas = length(tfaOpts);

%% parameters for Step 3 (calculating network instabilities w/ bStARS)
totSS = 50;
targetInstability = .05;
lambdaMin = 1e-5;
lambdaMax = 1;
extensionLimit = 1;
totLogLambdaSteps = 25; % will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5;
subsampleFrac = .9;
%% parameters for Step 4 (ranking TF-gene interactions)
meanEdgesPerGene = 25;
instabSource = 'Network';

for lind = 1:totLos  % you can use parfor loop here, if you have parallel computing toolbox  
    leaveOutSampleList = loInfo{lind,1};
    leaveOutInf = loInfo{lind,2};

for tind = 1:totTfas
    tfaOpt = tfaOpts{tind};
    
for lambdaBias = lambdaBiases

    instabilitiesDir = fullfile('HIOE_NEUROG3_induction/outputs',strrep(['instabilities_targ' ...
        num2str(targetInstability) '_SS' num2str(totSS) leaveOutInf '_bS' num2str(bStarsTotSS)],'.','p'));
    mkdir(instabilitiesDir)
    netSummary = [priorName '_bias' strrep(num2str(100*lambdaBias),'.','p') tfaOpt];
    instabOutMat = fullfile(instabilitiesDir,netSummary);

    %% 3. Calculate network instabilities using bStARS
    disp('3. estimateInstabilitiesTRNbStARS.m')
    estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
        totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,...
        subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)

    %% 4. For a given instability cutoff and model size, rank TF-gene
    % interactions, calculate stabilities and network file for jp_gene_viz
    % visualizations
    priorMergedTfsFile = ['HIOE_NEUROG3_induction/inputs/priors/' priorName '_mergedTfs.txt'];
    try % not all priors have merged TFs and merged TF files
        ls(priorMergedTfsFile) 
    catch
        priorMergedTfsFile = '';
    end
    networkDir = strrep(instabilitiesDir,'instabilities','networks');
    mkdir(networkDir);
    networkSubDir = fullfile(networkDir,[instabSource ...
        strrep(num2str(targetInstability),'.','p') '_' ...
        num2str(meanEdgesPerGene) 'tfsPerGene']);
    mkdir(networkSubDir)
    trnOutMat = fullfile(networkSubDir,netSummary);
    outNetFileSparse = fullfile(networkSubDir,[netSummary '_sp.tsv']);
    networkHistDir = fullfile(networkSubDir,'Histograms');
    mkdir(networkHistDir)
    subsampHistPdf = fullfile(networkHistDir,[netSummary '_ssHist']);

    disp('4. buildTRNs_mLassoStARS.m')
    buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile,...
        meanEdgesPerGene,targetInstability,instabSource,subsampHistPdf,trnOutMat,...
        outNetFileSparse)

    modSizes = 1:meanEdgesPerGene;
    r2OutMat = [instabOutMat '_r2pred'];

    disp('5. calcR2predFromStabilities')
    calcR2predFromStabilities(instabOutMat,trnOutMat,r2OutMat,modSizes)
    
end
end
end

%% 5. Calculate precision-recall relative to one or more gold standards
% to be added

