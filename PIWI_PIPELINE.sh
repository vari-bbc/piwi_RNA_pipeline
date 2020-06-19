#PBS -l walltime=10:00:00
#PBS -l mem=25gb
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -N PIWI_master

# CHANGE THESE AS NEEDED
dir=/my/working/dir/
get_length_distr_perlscript=${dir}AssignReadsToFeaturesWithLenDistr.pl
assign_counts_perlscript=${dir}AssignReadsToFeatures.pl
bowtie_index=/path/to/bowtie/index/ # e.g. DM6 UCSC
raw_data=${dir}raw_data/ # {$sample}_R1.fastq
mkdir -p ${dir}src/
trimmed=${dir}trimmed/
trimmed_suffix='_R1_trimmed.fq'
mkdir -p $trimmed
ppn=40 # processors per node for job submission
bowtie_out=${dir}bowtie/
mkdir -p ${bowtie_out}
#
# Reads will be assigned to these features in order
gtf_names=(rRNA tRNA snRNA snoRNA pre_miRNA repeatmasker piRNA_dm6 five_prime_utr three_prime_utr protein_coding_exon pseudogene ncRNA mitochondrion_genome)
gtf_path=/path/to/the/gtf/files/ # the gtf_names files are here and should be named e.g. rRNA.gtf, tRNA.gtf
strand=(sense antisense) # (1 2) argument to featureCounts
featureCounts_dir=${dir}featureCounts/
mkdir -p $featureCounts_dir

# this should be an array of all sample names which is passed to $assign_counts_perlscript
samples_argument=$(ls -1 ${raw_data}|grep _R1.fastq$|grep -v 'SV_PIWI_IP_5'| grep -v 'SV_Total_RNA_3'|grep -v 'SV_Total_RNA_4'|sed -e 's/_R1.fastq//g'|tr '\n' ','|sed -e 's/,$//') # updates 5/4/2020 to exclude 3 samples


#======================================================================= trim data:
echo "Generating trimming:"
i=0
for r1 in `ls -1 ${raw_data}|grep _R1.fastq$`; do # may need to change this
	sample=${r1/_R1.fastq/}

cat > ${dir}src/pbs.trimgalore.$i.sh << EOF

#PBS -l walltime=400:00:00
#PBS -l mem=10gb
#PBS -l nodes=1:ppn=1
#PBS -M ian.beddows@vai.org
#PBS -m a
#PBS -N ${sample}_trimgalore


cd ${trimmed}

trim_galore \
${raw_data}${r1} \
-q 20 \
--fastqc


EOF

	echo "	(${i})	${sample}	pbs.trimgalore.$i.sh"
	((i++))
done

#======================================================================= Generate feature counts table
echo "Generating bowtie scripts:"
echo "featureCounts with the following order:" && for i in "${gtf_names[@]}"; do echo "	$i"; done
i=0
for r1 in `ls ${trimmed}|grep _R1_trimmed.fq$`; do
	sample=${r1/trimmed-/}
	sample=${sample/_R1_trimmed.fq/}
	bam=${sample}.sorted.F4.bam

cat > ${dir}src/pbs.bowtie.$i.sh << EOF

#PBS -l walltime=400:00:00
#PBS -l mem=10gb
#PBS -l nodes=1:ppn=${ppn}
#PBS -M ian.beddows@vai.org
#PBS -m a
#PBS -N ${sample}_bowtie

module load bbc/bowtie/bowtie-1.2.3
module load bbc/subread/subread-2.0.0
cd ${bowtie_out}

if [ ! -f ${bam} ]; then
	echo "Mapping ..."
	bowtie \
	-v 1 \
	-p ${ppn} \
	-M 1 \
	--best \
	--strata \
	${bowtie_index} \
	${trimmed}${r1} \
	--seed 666 \
	--sam | samtools view -h -F 4 - |samtools sort -@ ${ppn} -O BAM -o ${bam} -
	#-F 4 filters out only mapped reads
fi
if [ ! -f ${bam}.bai ]; then
	echo "Indexing ..."
	samtools index ${bam}
fi

# now get the featureCounts for sense and antisense
gtf_names=(${gtf_names[@]})
strand=(${strand[@]})

for gtf in \${gtf_names[@]}; do
	for (( strand_arg=1; strand_arg<=${#strand[@]}; strand_arg++ )); do
		gtf_arg=${gtf_path}\${gtf}.gtf
		strand_name=\${strand[\${strand_arg}-1]}
		outfile=${featureCounts_dir}${sample}.\${strand_name}.\${gtf}.featureCounts
		if [ ! -f \$outfile ]; then
			echo "		Featurecounts - $sample"
			echo "			bam: ${bam}"
			echo "			gtf: \$gtf_arg"
			echo "			strand: \${strand_arg} (\${strand_name})"
			echo "			outfile: \${outfile}"

featureCounts \
-s \${strand_arg} \
-M \
-a \${gtf_arg} \
-F GTF \
-t gene \
-T ${ppn} \
-g gene_name \
-O \
-R CORE \
--Rpath ${featureCounts_dir} \
--fraction \
--minOverlap 15 \
-o \${outfile} \
${bam}

# this results in the file ${featureCounts_dir}${bam}.featureCounts
# this file will be overwritten for the different gtf files,
# so it is critical to rename it ...
OUTFILE=${featureCounts_dir}${sample}.\${strand_name}.\${gtf}.READINFO.featureCounts

mv ${featureCounts_dir}${bam}.featureCounts \${OUTFILE}

shuf \${OUTFILE} > \${OUTFILE}.rand

rm \${OUTFILE}

		else
			echo "Featurecounts - $sample - \${gtf_arg} - DONE!"
			echo "	featureCounts \
-s \${strand_arg} \
-M \
-a \${gtf_arg} \
-F GTF \
-t gene \
-T ${ppn} \
-g gene_name \
-O \
-R CORE \
--fraction \
--minOverlap 15 \
-o \${outfile} \
${bam}"
		fi




	done
done

# now all featureCounts have been generated.
# the next step is to iteratively assign reads to features
# based on the order of the gtf_names variable initialized earlier



EOF

	echo "	(${i})	$sample	${dir}src/pbs.bowtie.$i.sh"
	((i++))
done

#======================================================================= Assign featureCounts
echo "Assign Feature Counts"
gtf_argument=$(echo "${gtf_names[@]}"|sed -e 's/ /,/g')
file_extension='READINFO.featureCounts.rand'

for STRAND in "${strand[@]}"; do
	cmd="perl ${assign_counts_perlscript} $gtf_argument $samples_argument $STRAND ${featureCounts_dir} ${file_extension}"
	echo "	${dir}src/pbs.assign.${strand}.sh"
	cat > ${dir}src/pbs.assign.${STRAND}.sh << EOF
	#PBS -l walltime=400:00:00
	#PBS -l mem=100gb
	#PBS -l nodes=1:ppn=1
	#PBS -M ian.beddows@vai.org
	#PBS -m a
	#PBS -N ${STRAND}_assign_features_piwi

	echo "#======================================================================="
	echo "${STRAND}"
	$cmd
EOF

done
#======================================================================= Get length distribution
echo "Get Length Distribution"
for STRAND in "${strand[@]}"; do
	cmd="perl ${get_length_distr_perlscript} $gtf_argument $samples_argument $STRAND ${featureCounts_dir} ${file_extension} ${trimmed} ${trimmed_suffix}"
	echo "	${dir}src/pbs.length_distr.${STRAND}.sh"
	cat > ${dir}src/pbs.length_distr.${STRAND}.sh << EOF
	#PBS -l walltime=400:00:00
	#PBS -l mem=100gb
	#PBS -l nodes=1:ppn=1
	#PBS -M ian.beddows@vai.org
	#PBS -m a
	#PBS -N ${STRAND}_feature_length_piwi

	echo "#======================================================================="
	echo "${STRAND}"
	$cmd
EOF

done
exit
#======================================================================= Create stranded bigwig files
source ~/.bashrc
conda activate ucsc_bigwig
bigwig=${dir}bigwig/
mkdir -p $bigwig
for bam in `ls ${bowtie_out}|grep .bam$`; do
	outfile=${bam/.sorted.F4.bam/.forward.strand.cpm.bw}
	if [ ! -f ${outfile} ]; then
		echo "generating bigwig for $bam -> ${outfile}"
		bamCoverage -b ${bowtie_out}${bam} -o ${bigwig}${outfile} -p ${ppn} -of bigwig --normalizeUsing CPM --binSize 10 --filterRNAstrand forward
	fi
	outfile=${bam/.sorted.F4.bam/.reverse.strand.cpm.bw}
	if [ ! -f ${outfile} ]; then
		echo "generating bigwig for $bam -> ${outfile}"
		bamCoverage -b ${bowtie_out}${bam} -o ${bigwig}${outfile} -p ${ppn} -of bigwig --normalizeUsing CPM --binSize 10 --filterRNAstrand reverse
	fi
done
