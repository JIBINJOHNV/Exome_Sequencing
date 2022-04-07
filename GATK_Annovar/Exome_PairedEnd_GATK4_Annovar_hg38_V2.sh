Fasta="/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/Homo_sapiens_assembly38.fasta"
dbsnp='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/Homo_sapiens_assembly38.dbsnp138.vcf'
GoldIndels='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
hapmap='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/hapmap_3.3.hg38.vcf.gz'
Omni25='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/1000G_omni2.5.hg38.vcf.gz'
OneKgvcf='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/1000G_omni2.5.hg38.vcf.gz'
Annovar='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/annovar/'
AnnovarDb="/mnt/Data/NGS_Analysis_Database_Tools/Hg38/annovar/humandb_hg38/"




##Bed file; Depending on the target capture kit bed file need to change
Exome_bed='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/agilent_v8_targets_hg38.bed'
Exome_Vcfbed='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/agilent_v8_targets_hg38_10.bed'

#Exome_bed='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/agilent_v8_targets_hg38.bed'
#Exome_Vcfbed='/mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources/agilent_v8_targets_hg38_10.bed'


gvcfFolder="/mnt/Data/gvcf_all/" #Folder to save all the gvcf files generated from single sample vari	#ant calling
BamFolder="/mnt/Data/bam_all/" #Folder to save the final bam generated from single sample variant call	#ing


#sOFTWARES
picard='java -jar /mnt/Data/NGS_Analysis_Database_Tools/picard-tools-1.119/MarkDuplicates.jar'
BuildIndex='java -jar /mnt/Data/NGS_Analysis_Database_Tools/picard-tools-1.119/BuildBamIndex.jar'
Mosedepth_plots='python /mnt/Data/NGS_Analysis_Database_Tools/Mosedepth_plots/plot-dist.py'
MosedepthTargetComparison="/mnt/Data/NGS_Analysis_Database_Tools/MosedepthTargetComparison.py"
MosedepthTargetRegionBed_toPrecentage="/mnt/Data/NGS_Analysis_Database_Tools/MosedepthTargetRegionBed_toPrecentage.py"
RefinedGenotype="/mnt/Data/NGS_Analysis_Database_Tools/RefinedGenotype.py"
AnnovarFiltering="/mnt/Data/NGS_Analysis_Database_Tools/Hg38/annovar/AnnovarFiltering.py"



N=8 ; Number of jobs can be run parallel 
threads=5





#Creating Sample names based on fastq File
ls *.fastq.gz >fasqfile.txt && rm Sample_Names.txt temp.txt
cat fasqfile.txt
 
while read Fastq; do echo $Fastq | cut -d "_" -f 1 >> temp.txt ;done < fasqfile.txt

cat temp.txt |sort| uniq |sed '/^$/d' > Sample_Names.txt 
cat -n Sample_Names.txt 

#####################----------------------------------Analysis started---------------------------------------------------------------------

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

unset array && array=($(ls ${Samplename}_*fastq.gz)) 

echo "Step1: Prealignment QC" 

   echo "1.A) QC checking of Raw FastQ data using Fatqc"
mkdir -p QC/FastQC/Raw_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/ ${array[0]} & ##If fasq is not ending like '_R1.fq.gz'  it will not work
fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/ ${array[1]} & ##If fasq is not ending like '_R2.fq.gz'  it will not work


done < Sample_Names.txt 
wait 

mkdir -p MultiQC/FastQC/RawFastq_QC
multiqc -o MultiQC/FastQC/RawFastq_QC/ QC/FastQC/Raw_Fastq_QC/




echo "1.B) Remove low quality reads, adaptor and bases using TrimGalore OR fASTP (https://github.com/OpenGene/fastp#adapters)"

mkdir -p QC/TrimGalore.stdout_stats
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

unset array && array=($(ls ${Samplename}_*fastq.gz)) 

trim_galore --paired ${array[0]} ${array[1]}   \
--length 14 --length_1 15 --keep --length_2 15 \
--retain_unpaired -basename ${Samplename} --cores 8 >& QC/TrimGalore.stdout_stats/${Samplename}_TrimGalore.stdout_stats.txt &

mkdir -p QC/FastQC/TrimGalore_Fastq_QC/

done < Sample_Names.txt	;wait 







echo "1.C) FastQC of cleaned fastq file  "

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p QC/FastQC/TrimGalore_Fastq_QC/
fastqc -t ${threads} -o QC/FastQC/TrimGalore_Fastq_QC/ ${Samplename}_val_1.fq.gz &
fastqc -t ${threads} -o QC/FastQC/TrimGalore_Fastq_QC/ ${Samplename}_val_2.fq.gz &	
done < Sample_Names.txt ;wait 



mkdir -p MultiQC/FastQC/TrimGalore_Fastq_QC
multiqc -o MultiQC/FastQC/TrimGalore_Fastq_QC/ QC/FastQC/TrimGalore_Fastq_QC/




echo "1.D) Moving the TrimGalore Filtered/Trimed fastq files "

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p CleanedFastQ_Files/TrimGalore_QC/  
mv  ${Samplename}_val_1.fq.gz  CleanedFastQ_Files/TrimGalore_QC/ 
mv  ${Samplename}_val_2.fq.gz   CleanedFastQ_Files/TrimGalore_QC/ 
 	
mkdir -p QC_FailedFastQ

mv ${Samplename}*{1,2}.fastq.gz_trimming_report.txt QC/TrimGalore.stdout_stats/
mv ${Samplename}*{1,2}_unpaired_{1,2}.fq.gz QC_FailedFastQ

done < Sample_Names.txt ;wait # for all the something with stuff






echo "Step2: Alignment and Post Alignment Processing"

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
   echo "2.A) Alignment ;" 
ID="HTYYCBBXX.1.${Samplename}"
SM="${Samplename}"
LB=Lib-1
PU=HTYYCBBXX.1.TCAATCCG+TTCGCAGT
PL=ILLUMINA

#SM : give the name of the sample, CH_plate1_A01, the read is from,
#LB : gives the name of the library prep, Lib-1, the read is from,
#PU : gives the name of the flowcell, HTYYCBBXX, and the lane, 1, the read is from.

mkdir -p BamFiles
bwa mem  -R "@RG\tID:${Samplename}\tSM:${Samplename}\tLB:Lib\tPU:PU\tPL:Illumina" \
-M -t ${threads} \
${Fasta} \
CleanedFastQ_Files/TrimGalore_QC/${Samplename}_val_1.fq.gz \
CleanedFastQ_Files/TrimGalore_QC/${Samplename}_val_2.fq.gz | samtools sort -@ ${threads} -o BamFiles/${Samplename}.bam  \
 & done < Sample_Names.txt ; wait



#2.B) 

echo "mark duplicates ; https://gatk.broadinstitute.org/hc/en-us/articles/360057438771-MarkDuplicatesSpark; https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-"

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
   #gatk MarkDuplicatesSpark \
   #-I BamFiles/${Samplename}.bam \
   #-O BamFiles/${Samplename}_Mduplicates.bam \
   #-M BamFiles/${Samplename}_marked_dup_metrics.txt \
   #--create-output-bam-index true & done < Sample_Names.txt ; wait

${picard} I=BamFiles/${Samplename}.bam \
O=BamFiles/${Samplename}_Mduplicates.bam \
M=BamFiles/${Samplename}_marked_dup_metrics.txt  & done < Sample_Names.txt ; wait


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
${BuildIndex} I=BamFiles/${Samplename}_Mduplicates.bam O=BamFiles/${Samplename}_Mduplicates.bai & done < Sample_Names.txt ; wait





#echo "2.C) Base (Quality Score) Recalibration ; https://gatk.broadinstitute.org/hc/en-us/articles/360056969412-BaseRecalibrator"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
mkdir -p Recal_Files
gatk BaseRecalibrator \
-I BamFiles/${Samplename}_Mduplicates.bam \
-R ${Fasta} \
--known-sites ${dbsnp} \
--known-sites ${GoldIndels} \
-O Recal_Files/${Samplename}_recal_data.table & done < Sample_Names.txt ; wait


echo "2.D) ApplyBQSR ; https://gatk.broadinstitute.org/hc/en-us/articles/360056968652-ApplyBQSR"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk ApplyBQSR \
-R ${Fasta} \
-I BamFiles/${Samplename}_Mduplicates.bam \
--bqsr-recal-file Recal_Files/${Samplename}_recal_data.table \
-O BamFiles/${Samplename}_recal.bam \
--create-output-bam-index true & done < Sample_Names.txt ; wait

echo "2.E) Create Post recal table #https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk BaseRecalibrator \
-I BamFiles/${Samplename}_recal.bam \
-R ${Fasta} \
--known-sites ${dbsnp} \
--known-sites ${GoldIndels} \
--output  Recal_Files/${Samplename}_post_recal_data.table & done < Sample_Names.txt ; wait


echo "2.F) AnalyzeCovariates ; https://gatk.broadinstitute.org/hc/en-us/articles/360056967752-AnalyzeCovariates"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk AnalyzeCovariates \
-before Recal_Files/${Samplename}_recal_data.table \
-after Recal_Files/${Samplename}_post_recal_data.table \
-plots Recal_Files/${Samplename}_AnalyzeCovariates.pdf & done < Sample_Names.txt ; wait







###------------------------------------------------------------If required # shoud remove---------------------------------------------

#echo "Step3: Post Alignment Summary and Qc"

#while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
#mkdir -p QC/AlignmentQC/picard
#echo "3.A) CollectAlignmentSummaryMetrics ; https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/"
#gatk \
#      CollectAlignmentSummaryMetrics \
#      R=${Fasta} \
#      I=BamFiles/dups/${Samplename}_Mduplicates.bam \
#       O=QC/AlignmentQC/picard/${Samplename}_alignment_metrics.txt &

#echo "3.B) InsertSizeMetrics ; https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/"
#gatk \
 #   CollectInsertSizeMetrics \
  #  INPUT=BamFiles/dups/${Samplename}_Mduplicates.bam \
   # OUTPUT=QC/AlignmentQC/picard/${Samplename}_insert_metrics.txt \
   #HISTOGRAM_FILE=QC/AlignmentQC/picard/${Samplename}_insert_size_histogram.pdf


#mkdir -p QC/AlignmentQC/Mosdepth/GeneExone/
#mosdepth --by ${GeneExone_bed} QC/AlignmentQC/Mosdepth/GeneExone/${Samplename}_GeneExone BamFiles/dups/${Samplename}_Mduplicates.bam \
	#--thresholds 1,10,20,30,50 &

#mkdir -p QC/AlignmentQC/QualimapOutput/
#qualimap bamqc \
	 #  -bam BamFiles/dups/${Samplename}_Mduplicates.bam \
	  #  -gff ${QualiExomebed} \
	   #  -outdir  QC/AlignmentQC/QualimapOutput/ \
	    #  -outfile ${Samplename} \
	     #  --java-mem-size=16G &

#done < Sample_Names.txt ; wait

#mkdir -p MultiQC/AlignmentQC/picard
#multiqc -o MultiQC/AlignmentQC/picard/ QC/AlignmentQC/picard
#mkdir -p MultiQC/AlignmentQC/Mosdepth/GeneExone
#multiqc -o MultiQC/QC/AlignmentQC/Mosdepth/GeneExone/ QC/AlignmentQC/Mosdepth/GeneExone

#mkdir -p MultiQC/AlignmentQC/QualimapOutput/
#multiqc -o MultiQC/AlignmentQC/QualimapOutput/ QC/AlignmentQC/QualimapOutput/

#-------------------------------------------------------------------------------------------------------------------------------------------------------------





while read Samplename ; do ((i=i%N)); ((i++==0)) && wait

echo "3.C) Mos depth global and target region and gene exon"

mkdir -p QC/AlignmentQC/Mosdepth/Exome/


mosdepth -x --by ${Exome_bed} QC/AlignmentQC/Mosdepth/Exome/${Samplename}_Exome BamFiles/${Samplename}_Mduplicates.bam \
--thresholds 1,10,20,30,50,100 &


mkdir -p QC/AlignmentQC/flagstat
samtools flagstat BamFiles/${Samplename}_Mduplicates.bam > QC/AlignmentQC/flagstat/${Samplename}_flagstat.txt &

done < Sample_Names.txt ; wait


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
${Mosedepth_plots}  QC/AlignmentQC/Mosdepth/Exome/*.dist.txt --output QC/AlignmentQC/Mosdepth/Exome/${Samplename}_MosdepthCoverage.html 
done < Sample_Names.txt ;wait



mkdir -p MultiQC/AlignmentQC/Mosdepth/Exome
multiqc -o MultiQC/QC/AlignmentQC/Mosdepth/Exome/ QC/AlignmentQC/Mosdepth/Exome


mkdir -p MultiQC/CompleteQC
multiqc -o MultiQC/CompleteQC/ QC



##---------------------------------------Variant calling-----------------------------------------------------------------------"

echo "#Step4 ;  Variant calling ; "
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller ; 
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360057439091-StrandBiasBySample
                            #When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).
                            #https://gatk.broadinstitute.org/hc/en-us/articles/360035890551-Allele-specific-annotation-and-filtering-of-germline-short-variants


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait  
    echo "4.A) HaplotypeCaller; gvcf creation"
 gatk --java-options "-Xmx16g -XX:ParallelGCThreads=12"  HaplotypeCaller  \
   -R ${Fasta} \
   -I BamFiles/${Samplename}_recal.bam \
   -O ${Samplename}.g.vcf.gz \
   -ERC GVCF \
   -A StrandBiasBySample \
   --native-pair-hmm-threads 8 \
   -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation &

done < Sample_Names.txt ; wait 



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait  
echo "4.B) GenotypeGVCFs ; create single sample vcf file ; https://gatk.broadinstitute.org/hc/en-us/articles/360036899732-GenotypeGVCFs"
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=12"  GenotypeGVCFs \
 -R ${Fasta} \
 -V ${Samplename}.g.vcf.gz \
 -O ${Samplename}_raw.vcf.gz \
 -A StrandBiasBySample \
 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation & 
 done < Sample_Names.txt ; wait

##---------------------------------------Variant recalibration----------------------------------------------------------------------
#Step5 ;VariantRecalibrator  ; Not performing ; insted hard filter will do ; https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering




#normalize variants ;
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	bcftools norm \
	-f ${Fasta} \
	${Samplename}_raw.vcf.gz \
	--multiallelics -any \
	-o ${Samplename}_Laligned.vcf -O v  &
done < Sample_Names.txt ; wait



Subset to SNPs-only callset with SelectVariants
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	gatk SelectVariants \
	-V ${Samplename}_Laligned.vcf \
	-select-type SNP \
	-O ${Samplename}_snps.vcf.gz  & done < Sample_Names.txt ; wait


Subset to indels-only callset with SelectVariants
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	gatk SelectVariants \
	-V ${Samplename}_Laligned.vcf \
	-select-type INDEL \
	-O ${Samplename}_INDEL.vcf.gz  &  done < Sample_Names.txt ; wait




i#Hard-filter SNPs on multiple expressions using VariantFiltration
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	gatk VariantFiltration \
	    -V ${Samplename}_snps.vcf.gz \
	   -filter "QD < 2.0" --filter-name "QD2" \
           -filter "QUAL < 30.0" --filter-name "QUAL30" \
           -filter "SOR > 3.0" --filter-name "SOR3" \
           -filter "FS > 60.0" --filter-name "FS60" \
           -filter "MQ < 40.0" --filter-name "MQ40" \
           -filter "DP < 10.0" --filter-name "DP10" \
           -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
           -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O ${Samplename}_snps_filtered.vcf.gz & done < Sample_Names.txt ; wait


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	gatk VariantFiltration \
	    -V ${Samplename}_INDEL.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2" \
	    -filter "QUAL < 30.0" --filter-name "QUAL30" \
	    -filter "FS > 200.0" --filter-name "FS200" \
	    -filter "DP < 10.0" --filter-name "DP10" \
	    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	   -O ${Samplename}_INDEL_filtered.vcf.gz & done < Sample_Names.txt ; wait



#STEP16: Combine SNP and INDEL
	while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
 		gatk MergeVcfs  \
	        -R ${Fasta} \
	        -I "${Samplename}_"snps_filtered.vcf.gz \
	        -I "${Samplename}_"INDEL_filtered.vcf.gz \
	         --CREATE_INDEX true \
	        -O "${Samplename}_"hardfiltered.vcf.gz & done < Sample_Names.txt ; wait




#STEP16: Combine SNP and INDEL
	while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
		gatk MergeVcfs  \
		-R ${Fasta} \
		-I "${Samplename}_"snps_filtered.vcf.gz \
		-I "${Samplename}_"INDEL_filtered.vcf.gz \
		--CREATE_INDEX true \
		-O "${Samplename}_"hardfiltered.vcf.gz & done < Sample_Names.txt ; wait



#Remove the intermediate files
	while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	rm ${Samplename}_Laligned.vcf ${Samplename}_snps.vcf.gz*  ${Samplename}_INDEL.vcf.gz*
	rm ${Samplename}_snps_filtered.vcf.gz* ${Samplename}_INDEL_filtered.vcf.gz* ${Samplename}_raw.vcf.gz* & done < Sample_Names.txt ; wait




#Select variants from target region
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
	gatk SelectVariants \
       -V "${Samplename}_"hardfiltered.vcf.gz \
        -L ${Exome_Vcfbed} \
        -O "${Samplename}_"hardfiltered_Target.vcf
        done < Sample_Names.txt ; wait


	
#----------------------------------------------------Not performing this steps if required uncomment the commands-------------------------------------------
##EXTRACT ONLY PASSED VARIANTS
#while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
#	        vcffilter -f "FILTER = PASS" "${Samplename}_"hardfiltered_Target.vcf  > "${Samplename}_"hardfiltered_Target_Pass.vcf ; done < Sample_Names.txt ; wait

#------------------------------------------------------------------------------------------------------------------------------------------------------------------


###ANNOTATION OF VCF FILE USING aNNOVAR
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 ${SAMPLENAME}_hardfiltered_Target.vcf \
		-outfile ${SAMPLENAME}.avinput \
		-includeinfo -withzyg &
done < Sample_Names.txt ; wait


while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait

perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME} \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done < Sample_Names.txt ; wait


##Filter annovar output file based on the defined parameters
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait 
	python ${AnnovarFiltering} --AnnovarFile ${SAMPLENAME}.hg38_multianno.txt &
done < Sample_Names.txt ; wait



#Edit column headers and Incorporate additional columns
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait 
	find="Otherinfo4"
	replace=$(grep "#CHROM"  ${SAMPLENAME}_hardfiltered_Target.vcf)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv 


	find="Otherinfo4"
	replace=$(grep "#CHROM" ${SAMPLENAME}_hardfiltered_Target.vcf)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv 
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv 

done < Sample_Names.txt ; wait

rm *.avinput *.hg38_multianno.txt *raw.vcf.gz*



##Additional Analysis
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait
	
#Exome region based coverage %
python ${MosedepthTargetComparison} -Mfolder QC/AlignmentQC/Mosdepth/Exome/ -Outputfolder  QC/AlignmentQC/Mosdepth/Exome/

##dEFINE gENOTYPE(rEFINED GENOTYPE ) BASED ON THE DEFINED PARAMETERS
python ${MosedepthTargetRegionBed_toPrecentage} -Mfolder QC/AlignmentQC/Mosdepth/Exome/

python ${RefinedGenotype} --AnnovarOutputFile ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv

python ${RefinedGenotype} --AnnovarOutputFile ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv
done < Sample_Names.txt ; wait




#Remove Intermediate files and move gvcf and bam file to respective folders and create symbolic link
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait

mv BamFiles/${SAMPLENAME}_recal.ba* ${BamFolder}
mv ${SAMPLENAME}.g.vcf.gz* ${gvcfFolder}

ln -s ${BamFolder}${SAMPLENAME}_recal.bam ${SAMPLENAME}_recal.bam
ln -s ${BamFolder}${SAMPLENAME}_recal.bai ${SAMPLENAME}_recal.bai

done < Sample_Names.txt ; wait


rm -rf fasqfile.txt temp.txt *_hardfiltered.vcf.gz* BamFiles *idx
rm -rf CleanedFastQ_Files QC_FailedFastQ






