#!/bin/bash

# core_priorsSparse2Table.sh

# use prior-parsing code to generate a prior matrix from the sparse network (to enable
# downstream analyses, e.g. TF-TF module analysis)

# INPUTS
# script relies on prior-parsing code
scriptHome="../priorConstruction"

# condition/time
cond_list=(24 48 72 96)
for cond_time in ${cond_list[@]}
do
	inputDir="outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/${cond_time}hpi_Cores"

	declare -a inFileNames=("Core_prior_atac_Miraldi_q_ChIP_bias10_maxComb_fdr5_HIOE_NEUROG3_${cond_time}hpiSet_All_sp.tsv" "Core_prior_atac_Miraldi_q_ChIP_bias10_maxComb_fdr5_HIOE_NEUROG3_${cond_time}hpiSet_DE_sp.tsv")
	outputDirBase=${inputDir}
# END INPUTS

	for inFileName in ${inFileNames[@]}
	do
		inFile=${inputDir}/${inFileName}
		wc -l ${inFile}

		# convert filtered network to square form
		sqFull=${inFile/_sp/}
	 	python ${scriptHome}/priorsSparse2Table.py ${inFile} ${sqFull} 0
		wc -l ${sqFull}
	done

done
echo "Finished!"
