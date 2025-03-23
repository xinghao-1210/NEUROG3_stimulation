#!/bin/bash

# core_priorsTable2Sparse.sh

# generate the  prior matrix (Integrated with timepoint ChIP data) from a sparse network (_sp) to enable
# downstream analyses, e.g. jp_gene visualization with ChIP integratation

# INPUTS
# script relies on prior-parsing code
scriptHome="../priorConstruction"

# condition/time
cond_list=(24)
for cond_time in ${cond_list[@]}
do
	inputDir="outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb/${cond_time}hpi_Cores"

	declare -a inFileNames=("Core_prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb_fdr5_HIOE_NEUROG3_${cond_time}hpiSet_All_ChIP.tsv" "Core_prior_atac_Miraldi_q_ChIP_x10_bias10_maxComb_fdr5_HIOE_NEUROG3_${cond_time}hpiSet_DE_ChIP.tsv")
	outputDirBase=${inputDir}
# END INPUTS

	for inFileName in ${inFileNames[@]}
	do
		inFile=${inputDir}/${inFileName}
		wc -l ${inFile}

		# convert filtered network to square form
		sqFull=${inFile/.tsv/_sp.tsv}
	 	python ${scriptHome}/priorsTable2Sparse.py ${inFile} ${sqFull} 0
		wc -l ${sqFull}
	done

done
echo "Finished!"
