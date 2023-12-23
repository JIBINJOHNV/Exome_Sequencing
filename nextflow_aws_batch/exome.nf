params.publish_dir="s3://rawgdata/Results/"
params.samples_tsv ='/home/ec2-user/environment/single_end_fastq_files_locations.csv'

params.genome="/home/ec2-user/environment/human_index_file.csv"
params.dbsnp="s3://g-reference/Homo_sapiens_assembly38.dbsnp138.vcf"
params.goldindels="s3://g-reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.hapmap="s3://g-reference/hapmap_3.3.hg38.vcf.gz"
params.omni25="s3://g-reference/1000G_omni2.5.hg38.vcf.gz"
params.oneKgvcf="s3://g-reference/1000G_omni2.5.hg38.vcf.gz"


log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    genome        : ${params.genome}
    dbsnp         : ${params.dbsnp}
    goldindels    : ${params.goldindels}
    hapmap        : ${params.hapmap}
    omni25        : ${params.omni25}
    oneKgvcf      : ${params.oneKgvcf}
    """
    .stripIndent(true)



process FastQC {
    tag { FASTQ: meta.id }
    publishDir "${params.publish_dir}/QC/FastQC/${meta.id}", pattern: "*html", mode: 'copy'
    publishDir "${params.publish_dir}/QC/FastQC/${meta.id}", pattern: "*zip", mode: 'copy'

    input:
    tuple  val(meta), path(reads)

    output:
    path '*'

   script:
        """
        fastqc -t 8 -q  ${reads[0]} ${reads[1]}
    
        """ 
}


process FASTP {
 
    tag { FASTQ: meta.id }
    publishDir "${params.publish_dir}/CleanedFastQ_Files/Fastp", pattern: "*fastp_R*.gz", mode: 'copy'
    publishDir "${params.publish_dir}/QC/Fastp_stdout_stats/", pattern: "*_fastp.html" , mode: 'copy'
    publishDir "${params.publish_dir}/QC/Fastp_stdout_stats/", pattern: "*fastp.json" , mode: 'copy'

    input:
    tuple  val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastp_R*.gz"), emit: trimmed_fastqs
    path("${meta.id}_fastp.json"), emit: json
    path("${meta.id}_fastp.html"), emit: html

   script:
        """
	       fastp --thread 4 \
            --in1 ${reads[0]} --in2 ${reads[1]} --out1 ${meta.id}_fastp_R1.gz --out2 ${meta.id}_fastp_R2.gz \
            --cut_front cut_front_window_size 3 --cut_front_mean_quality 20 --cut_tail cut_tail_window_size 3 --cut_tail_mean_quality 20 \
            --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 30 --compression 4 --trim_poly_x  --correction   \
            --json ${meta.id}_fastp.json --html ${meta.id}_fastp.html 
        
        """ 
}


process BWA {

    debug true

    // Publish directories for output files
    publishDir "$params.publish_dir/Raw_Bam/", pattern: "*.bam", mode: 'copy'


    input:
    tuple val(meta), path(reads)
    path(bwa_index_ch)


    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam



    script:
    """
    bwa mem -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:Lib-1\\tPU:HTYYCBBXX.1\\tPL:ILLUMINA" -M -t 30 ${bwa_index_ch[0]} \
    ${reads[0]} ${reads[1]} | samtools sort -o ${meta.id}.bam
    """
}


workflow { 

    ch_samplesheet = Channel.fromPath(params.samples_tsv)
    fastq_channel=ch_samplesheet.splitCsv( header: true, sep: ',' )
        .map {
        R1 = it['R1']
        R2 = it['R2']
        is_singleEnd = R2.toString()=='' ? true : false
        meta = [id: it['samplename'], single_end: is_singleEnd]
        R2.toString()=='' ? [meta, [R1]] : [meta, [R1, R2]]
        }

    ch_genome_sheet = Channel.fromPath(params.genome)
    ch_genome=Channel.fromPath(params.genome).splitCsv( header: true, sep: ',' )
                | map { row ->tuple(row.location)}
    ch_genome=ch_genome.collect()
    
    
    //FastQC(fastq_channel)
    Fastp_ch=FASTP(fastq_channel)
    bwa_ch=BWA(Fastp_ch.trimmed_fastqs,ch_genome)

}
