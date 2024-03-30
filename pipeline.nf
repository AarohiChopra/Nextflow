#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
 
// input params
params.inputFile = "/home/achopra/BPA_Alt_Human/BPA/Fastq_dump/input/SRR_Acc_List.txt"
params.fd_outDir = "/home/achopra/BPA_Alt_Human/BPA/Fastq_dump/"
params.qc_outdir = "/home/achopra/BPA_Alt_Human/BPA/Fast_qc/"
params.qc_trimming = "/home/achopra/BPA_Alt_Human/BPA/Fast_p/"
params.multiqc = "/home/achopra/BPA_Alt_Human/BPA/Multi_qc/"
params.base_dir = "/home/achopra/BPA_Alt_Human/BPA/"
params.genome_dir = "/home/achopra/BPA_Alt_Human/BPA/Mapping/"
params.star_output = "/home/achopra/BPA_Alt_Human/BPA/STAR_output/"
params.sort_by_coordinate = "/home/achopra/BPA_Alt_Human/BPA/Picard_coordinate/"
params.mark_duplicate = "/home/achopra/BPA_Alt_Human/BPA/Mark_duplicates/"
params.picard_jar_path = "/opt/picard/picard.jar"
params.feature_counts = "/home/achopra/BPA_Alt_Human/BPA/Feature_counts/"
params.annotation_path="/home/achopra/BPA_Alt_Human/BPA/Human_ref_files/hg38.ncbiRefSeq.gtf"


log.info """\
    =====================================================
       R N A S E Q - T E M P O S E Q   P I P E L I N E
    =====================================================
	 READS : ${params.inputFile}
	 """.stripIndent()
// Create a channel from the file at params.inputFile
// Split the contents of the file into separate lines
// trim whitespace from each line  
// Assigns name to the channel 
Channel
    .fromPath(params.inputFile)
    .splitText()
    .map { it.trim() }
    .set { FastqDump_ch }
    
process FastqDump {
    label 'Download_Raw_Reads' // assigns name to the process
    cpus 25
    clusterOptions "--nodes=1" // each job takes on node 
    tag "on ${sra_id}" // Tags each process execution with a specific identifier
    publishDir "${params.fd_outDir}${sra_id}", mode: 'move', pattern: "*.fastq" // Specifies the dir where the results should be saved                      
    
    input:
    val sra_id // no from FastqDump_ch not allowed in dsl 2
    
    output:
    tuple val(sra_id), path("${sra_id}.fastq")
    
    script:
    """
    echo "Downloading: ${sra_id}" 
    mkdir -p ${params.fd_outDir}${sra_id}
    module load sratoolkit-3.0.2
    prefetch ${sra_id}
    fasterq-dump ${sra_id} -e 4 --skip-technical
    """
}
process QualityCheck {
    label "Quality_Checks" 
    cpus 25
    clusterOptions "--nodes=1" 
    tag "on ${sra_id}"
    publishDir "${params.qc_outdir}${sra_id}", mode:'copy'
    
    input:
    tuple val(sra_id). path(fastqd)
    
    output:
    tuple val(sra_id), path("${sra_id}_fastqc.html"), path("${sra_id}_fastqc.zip")
    
    script:
    """
    echo "Processing fastqc file: ${sra_id}" 
    mkdir -p ${params.qc_outdir}${sra_id}
    module load fastqc
    fastqc ${params.fd_outDir}${sra_id}/${sra_id}.fastq --outdir .
    """
}
process FastpTrimming {
    label "Trimming_Adapter" 
    cpus 25
    clusterOptions "--nodes=1" 
    tag "on ${sra_id}"
    publishDir "${params.qc_trimming}${sra_id}", mode:'move'
    conda "/home/achopra/miniconda3/envs/fastp/fastp.yaml"
    
    input:
    tuple val(sra_id), path(fastqc), path(fastqczip)
    
    output:
    tuple  val(sra_id), path("${sra_id}_fastp_trimmed.fastq"), path("${sra_id}_fastp.html"), path("${sra_id}_fastp.json")
    
    script:
    """
    mkdir -p ${params.qc_trimming}${sra_id}
    fastp \
        -i ${params.fd_outDir}${sra_id}/${sra_id}.fastq \
        -o ${sra_id}_fastp_trimmed.fastq \
        -h ${sra_id}_fastp.html \
        -j ${sra_id}_fastp.json
    """
}
process MultiQCFastqcData {
    label "MultiQC_On_QCfiles" 
    cpus 25
    clusterOptions "--nodes=1" 
    publishDir "${params.multiqc}", mode:'move'
    conda "/home/achopra/miniconda3/envs/multiqc/multiqc.yaml"
    
    input:
    path(filePaths)
    
    output:
    path("multiqc_FastQC_report.html")
   
    script:
    """
    multiqc ${filePaths} -n multiqc_FastQC_report.html -o . 
    """
}
process RNAStar {
    label "RNA_Star_Aligning"
    cpus 25
    clusterOptions "--nodes=1"
    publishDir "${params.star_output}${sra_id}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/bpaenv/bpaenv.yaml"
    memory '60 GB'
    
    input:
    tuple val(sra_id), path(trimmed_fastq), path(ignore1), path(ignore2)
    
    output:
    tuple val(sra_id), path("${sra_id}Log.out"), path("${sra_id}Log.progress.out"), path("${sra_id}Log.final.out"), path("${sra_id}SJ.out.tab"), path("${sra_id}Aligned.out.bam")
    
    script:
    """
    mkdir ${params.star_output}${sra_id}
    STAR --runMode alignReads \
         --genomeDir ${params.genome_dir} \
         --runThreadN 4 \
         --readFilesIn "${params.qc_trimming}${sra_id}/${sra_id}_fastp_trimmed.fastq" \
         --outFileNamePrefix "${sra_id}" \
         --outTmpDir ${params.base_dir}RNAscratch${sra_id} \
         --outSAMtype BAM Unsorted
    rm -rf "${params.base_dir}RNAscratch${sra_id}"
    """
}
process MultiQCAlignedData {
    label "MultiQC_On_Alignedfiles" 
    cpus 25
    clusterOptions "--nodes=1" 
    publishDir "${params.multiqc}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/multiqc/multiqc.yaml"
    
    input:
    path(filePaths)
    
    output:
    path("multiqc_Aligned_report.html")
   
    script:
    """
    multiqc ${filePaths} -n multiqc_Aligned_report.html -o . 
    """
}
process SortByCoordinate {
    label "Sort_By_Coordinate"
    cpus 25
    clusterOptions "--nodes=1"
    tag "on ${sra_id}"
    publishDir "${params.sort_by_coordinate}${sra_id}", mode:'move'
    
    input:
    tuple val(sra_id), path(ignore1), path(ignore2), path(ignore3), path(ignore4), path(ignore5)
    
    output:
    tuple val(sra_id), path("${sra_id}SC.bam")
    
    script:
    """
    module load java
    module load picard
    mkdir ${params.sort_by_coordinate}${sra_id}
    java -jar ${params.picard_jar_path} SortSam \
    INPUT="${params.star_output}${sra_id}/${sra_id}Aligned.out.bam" \
    OUTPUT="${sra_id}SC.bam" \
    SORT_ORDER="coordinate"
    """
}
process MarkDuplicates {
    label "Mark_Duplicates"
    cpus 25
    clusterOptions "--nodes=1"
    tag "on ${sra_id}"
    publishDir "${params.mark_duplicate}${sra_id}", mode:'copy'
    
    input:
    tuple val(sra_id), path(ignore1)
    
    output:
    tuple val(sra_id), path("${sra_id}MD.bam"), path("${sra_id}MD.txt")
    
    script:
    """
    module load java
    module load picard
    mkdir ${params.mark_duplicate}${sra_id}
    java -jar ${params.picard_jar_path} MarkDuplicates \
    INPUT="${params.sort_by_coordinate}${sra_id}/${sra_id}SC.bam" \
    OUTPUT="${sra_id}MD.bam" \
    M="${sra_id}MD.txt"
    """
}
process MultiQCBAMData {
    label "MultiQC_On_BAMfiles" 
    cpus 25
    clusterOptions "--nodes=1" 
    publishDir "${params.multiqc}", mode:'move'
    conda "/home/achopra/miniconda3/envs/multiqc/multiqc.yaml"
    
    input:
    path(filePaths)
    
    output:
    path("multiqc_BAM_report.html")
   
    script:
    """
    multiqc ${filePaths} -n multiqc_BAM_report.html -o . 
    """
}
process FeatureCounts {
    label "Feature_Counts"
    cpus 25
    clusterOptions "--nodes=1"
    tag "on ${sra_id}"
    publishDir "${params.feature_counts}", mode:'move'
    conda "/home/achopra/miniconda3/envs/featurecounts/featurecounts.yaml"
    
    input:
    path(bamFiles)
    
    output:
    path("combined_feature_counts.txt")
    
    script:
    """
    featureCounts -T 5 \
              -t exon \
              -g gene_id \
              -a "${params.annotation_path}" \
              -o "combined_feature_counts.txt" ${bamFiles}
    """
}
workflow {
    Fastq_Dump_Files = FastqDump(FastqDump_ch)
    Fastq_QC_Files = QualityCheck(Fastq_Dump_Files)
    
    Trimmed_files = FastpTrimming(Fastq_QC_Files)
    
    QC_Zip_Files = QualityCheck.out.collect(){ items -> items.findAll { it != null && it.toString().endsWith(".zip") }}
    MultiQC_file_qc_data = MultiQCFastqcData(QC_Zip_Files)
    
    Aligned_files = RNAStar(Trimmed_files)
    
    Aligned_Zip_Files = RNAStar.out.collect(){ items -> items.findAll { it != null && it.toString().endsWith("Log.final.out") }}
    MultiQC_file_aligned_data = MultiQCAlignedData(Aligned_Zip_Files)
    
    Sort_by_coordinate_files = SortByCoordinate(Aligned_files)
    
    Mark_duplicate_files = MarkDuplicates(Sort_by_coordinate_files)
    
    Mark_duplicate_BAM_files = MarkDuplicates.out.collect(){ items -> items.findAll { it != null && it.toString().endsWith(".txt") }}
    MultiQC_file_bam_data = MultiQCBAMData(Mark_duplicate_BAM_files)
    
    feature_count_File = MarkDuplicates.out.collect(){ items -> items.findAll { it != null && it.toString().endsWith(".bam") }}
    FeatureCounts(feature_count_File)
    
    workflow.onComplete {
    println("Pipeline completed successfully.")
}
}
