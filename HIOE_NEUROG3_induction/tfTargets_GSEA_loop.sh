#!/bin/bash
# tfTarget_GSEA_loop.sh
# given databases of gene sets, test for enrichment of each set in the targets of each TF
# in the network. Positive and negative edges can be considered separately (see signOpts
# below). GSEA results are further organized into tables (gene sets X TFs) of raw and
# adjusted p-values, as well as other statistics, for downstream visualization (e.g., in
# MATLAB)

# Note that geneset databases must be formatted as .txt files with background files
# following naming convention below
setInDir="inputs/geneLists/gseaDatabases/human"
bckgrndBit="_bckgrnd.txt"
setBit=".txt"
declare -a setNames=('GO_all' 'CP_all')

# condition/time
cond_list=(24 48 72 96)
for cond_time in ${cond_list[@]}
do
	# Network files (NOTE: for signed-analysis of TF targets, we work with TRNs where edges
	# have been filtered based on partial correlation, e.g., as in filter_TRNs_by_pcorr.sh
	declare -a networkFiles=("outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/${cond_time}hpi_Cores/Core_prior_atac_Miraldi_q_ChIP_bias10_maxComb_fdr5_HIOE_NEUROG3_${cond_time}hpiSet_All_sp.tsv" \
													 "outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/${cond_time}hpi_Cores/Core_prior_atac_Miraldi_q_ChIP_bias10_maxComb_fdr5_HIOE_NEUROG3_${cond_time}hpiSet_DE_sp.tsv")

	scriptHome="../infLassoStARS"

	# Parameters for GSEA -- see tfTargetGSEA_sortPvals.py for definitions
	interactionWeightCut="0"
	# In initial GSEA, include an enrichment with a raw p-value <= to this cutoff
	rawPcut=".1"
	# Only test gene sets that overlap with our background gene set at this level:
	minSetSize=5
	# TF-gene interactions are signed, and so a TF's postive and negative targets can be
	# analyzed separately (signed or "dir" for directionally, signOpt = 1) or together
	# (unsigned, "abs" for absolute value of interaction, signOpt = 0). Signed analysis makes
	# sense in terms of finding pathways that are turned on or off (note this could be
	# imperfect if some genes in the pathway are negative regulators). Unsigned analysis makes
	# sense for GWAS-derived gene sets, where the genetic effect on gene function is often unknown.
	declare -a signOpts=("1" "0")
	# For putting results in table format for downstream visualization
	fdrCut="1"
	rawPcutTable=".1"

	function log {
		echo $(date +%F_%T) $$ $BASHPID $1
		}

	for geneSetAbbrev in ${setNames[@]}
	do
		geneSetFile=${setInDir}/${geneSetAbbrev}${setBit}
		backgroundGeneFile=${setInDir}/${geneSetAbbrev}${bckgrndBit}

		for networkFile in ${networkFiles[@]}
		do
			echo $networkFile
			fileBaseName=$(basename ${networkFile} _sp.tsv)
			# -- base for output file names
			outDir=$(dirname ${networkFile})/${fileBaseName}/GSEA
			mkdir -p ${outDir}
			echo ${outDir}
			# -- directory for output

			outDirBase=${outDir}/${fileBaseName}_${geneSetAbbrev}
			echo ${tfTargetOutDir}
			# output directory for motifEnrichemnts_TFtargets.py

			for signOpt in ${signOpts[@]}
			do
				if [ "${signOpt}" -eq "1" ]; then
					signInf="dir"
					echo "signed GSEA of TF targets"
				elif [ "${signOpt}" -eq "0" ]; then
					signInf="abs"
					echo "unsigned GSEA of TF targets"
				else
					echo "signOpt must be either 0 or 1 (for unsigned/signed GSEA)"
					exit 64
				fi

				log "GSEA of ${netBase} WRT ${geneSetAbbrev}"
				echo "python ${scriptHome}/tfTargetGSEA_sortPvals.py ${networkFile} ${geneSetFile} ${backgroundGeneFile} ${rawPcut} ${signOpt} ${interactionWeightCut} ${minSetSize} ${outDirBase}"
				python ${scriptHome}/tfTargetGSEA_sortPvals.py ${networkFile} ${geneSetFile} \
					${backgroundGeneFile} ${rawPcut} ${signOpt} ${interactionWeightCut} \
					${minSetSize} ${outDirBase}

				gseaFileName=${fileBaseName}_${geneSetAbbrev}_Praw${rawPcut/\./0p}_${signInf}_wCut0p0_minSet${minSetSize}
				ls ${outDir}/${gseaFileName}.txt

				log "Convert to table form, based on FDR cutoff -- for unsigned analysis"
				echo "python ${scriptHome}/makeGSEAmatrices_sepPosNeg_fxn.py ${outDir} ${fdrCut} ${outDir} ${gseaFileName} ${geneSetAbbrev}"
				python ${scriptHome}/makeGSEAmatrices_sepPosNeg_fxn.py \
					${outDir} \
					${fdrCut} \
					${outDir} \
					${gseaFileName} \
					${geneSetAbbrev}

				log "Convert to table form, based on raw p-value cutoff -- for unsigned analysis"
				gseaFileName=${fileBaseName}_${geneSetAbbrev}_Praw${rawPcutTable/\./0p}_${signInf}_wCut0p0_minSet${minSetSize}
				ls ${outDir}/${gseaFileName}.txt
				echo "python ${scriptHome}/makeGSEAmatrices_sepPosNeg_prawCut_fxn.py \
					${outDir} \
					${rawPcutTable} \
					${outDir} \
					${gseaFileName} \
					${geneSetAbbrev}"

				python ${scriptHome}/makeGSEAmatrices_sepPosNeg_prawCut_fxn.py \
					${outDir} \
					${rawPcutTable} \
					${outDir} \
					${gseaFileName} \
					${geneSetAbbrev}

			done
		done
	done
done

echo "Finished!"
