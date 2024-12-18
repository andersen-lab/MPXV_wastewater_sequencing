nextflow.enable.dsl=2

params.reads = '/path-to-fastq-folder/*.fastq'
params.reference = './mpxv_ref.fasta'
params.bed_file = './path-to-bed-file/MPXV_v1.1.bed'
params.outdir = '/demix_outdir/'
params.bin_path = '/path-to-dorado-bin/'

workflow {
    // Create a channel from the FASTQ files
    reads = Channel.fromPath(params.reads, checkIfExists: true)
    // Transform the channel to pad the sample name
    reads_tuple = reads.map { file ->
        def sample_id = file.baseName
        tuple(sample_id, file)
    }.filter { tuple ->
        def sample_id = tuple[0]
        // Exclude samples with 'unclassified' or 'CustomBarcodes' in the sample_id
        !sample_id.contains('unclassified') && !sample_id.contains('CustomBarcodes')
    }
    trimmed_fastq = processTrim(reads_tuple)
    reads_bed = trimmed_fastq.map { sample_id, reads -> tuple(sample_id, reads, file(params.bed_file)) }
    bam = processAlign(reads_bed)
    bam_bed = bam.map { sample_id, bam_file, bai_file -> tuple(sample_id, bam_file, bai_file, file(params.bed_file)) }
    read_depth = processCalculateReadDepth(bam_bed)
    processFreyja(bam_bed)
    processCombineCSV(read_depth.collect())
}

process processTrim {
    conda "environment.yml"

    publishDir "${params.outdir}/trimmed_fastq/", mode: 'copy', pattern: "${sample_id}_trimmed.fastq"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq")

    script:
    """
    seqkit seq -Q 10 ${reads} -o ${sample_id}_trimmed.fastq
    """
}

process processAlign {
    conda "environment.yml"

    publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "${sample_id}_trimmed_sorted.bam"    
    publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "${sample_id}_trimmed_sorted.bam.bai"

    input:
    tuple val(sample_id), path(reads), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_sorted.bam"), path("${sample_id}_trimmed_sorted.bam.bai")

    script:
    """
    ${params.bin_path}/dorado aligner ${params.reference} ${reads} | ivar trim -b ${bed} -q 0 -m 200 | samtools sort -o ${sample_id}_trimmed_sorted.bam
    samtools index ${sample_id}_trimmed_sorted.bam
    """
}

process processCalculateReadDepth {
    conda "environment.yml"
    
    publishDir "${params.outdir}/amplicon_depths/", mode: 'copy', pattern: "${sample_id}_read_depth.csv"

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(bed_file)

    output:
    path "${sample_id}_read_depth.csv"

    script:
    """
    calculate_amplicon_depth.py --bam ${bam_file} --bed ${bed_file} --out ${sample_id}_read_depth.csv
    """
}

process processCombineCSV {
    conda "environment.yml"

    publishDir "${params.outdir}/amplicon_depths/", mode: 'copy'

    input:
    path csv_files

    output:
    path "combined_read_depth.csv"

    script:
    """
    combine_csvs.py ${csv_files.join(" ")} combined_read_depth.csv
    """
}

process processFreyja {
    conda "environment.yml"
    
    publishDir "${params.outdir}/freyja/", mode: 'copy', pattern: "${sample_id}.demixed.tsv"

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed_file)

    output:
    tuple val(sample_id), path("${sample_id}.demixed.tsv")

    script:
    """
    freyja update --pathogen MPXV
    freyja variants ${bam} --variants ${sample_id}.tsv --depths ${sample_id}.depth --minq 0 --ref ${params.reference}
    freyja demix ${sample_id}.tsv ${sample_id}.depth --output "${sample_id}.demixed.tsv" --pathogen MPXV --depthcutoff 1 --adapt 0.035
    """
}
