%% combine_Th17_TRNs
% max combine TFA and TF mRNA Th17 TRNs (15 TFs per gene) to get a 15 TFs per
% gene model

clear all
close all
restoredefaultpath

currDir = '';

addpath(fullfile(currDir,'infLassoStARS'))
addpath(fullfile(currDir,'glmnet'))
addpath(fullfile(currDir,'customMatlabFxns'))

%% parameters

meanEdgesPerGene = 15;

combine_list = {'max','mean'};
for combine_ind = 1:length(combine_list) 
    combineOpt = combine_list{combine_ind};

    combinedNetTsv = sprintf('pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_%dtfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_%sComb_sp.tsv',10,combineOpt);

    nets2combine = {sprintf('pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_%dtfsPerGene/prior_atac_Miraldi_q_ChIP_bias10.mat',10);
                    sprintf('pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_%dtfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_TFmRNA.mat',10)};

    combineTRNs(combinedNetTsv,combineOpt,meanEdgesPerGene,nets2combine)
end