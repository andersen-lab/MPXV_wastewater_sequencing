nextflow.enable.dsl=2

params.reads = '/path-to-pod5-folder/*.pod5'
params.samplesheet = '/path-to-sample-sheet/SampleSheet.csv'
params.model = './dna_r10.4.1_e8.2_400bps_hac@v4.3.0'
params.reference = './mpxv_ref.fasta'
params.kitname = 'CustomBarcodes'
params.outdir = '/outdir/'
params.bin_path = '/home/alab/dorado-0.8.1-linux-x64/bin/'
params.toml = './CustomBarcodes.toml'
params.barcode_seqs = './CustomBarcodes.fasta'


workflow {
    pod5_files = Channel.fromPath(params.reads, checkIfExists: true)

    // Collect basecalled BAM files before merging
    reads = processBasecalling(pod5_files)
    merged_bam = processMergeBAM(reads.collect())
    demux_results = processDemux(merged_bam)
}

process processBasecalling {
    maxForks 1

    input:
    path pod5_file

    output:
    path "*.bam"

    shell:
    '''
    export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:30
    !{params.bin_path}/dorado basecaller -x cuda:0 --no-trim !{params.model} !{pod5_file} > !{pod5_file.baseName}.bam
    '''
}

process processMergeBAM {
    input:
    path(bam_files)

    output:
    path "merged_sorted.bam"

    shell:
    '''
    # Merge and sort BAM files in one command
    samtools merge -O bam - !{bam_files.join(' ')} | samtools sort -o merged_sorted.bam
    '''
}

process processDemux {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path merged_bam_file

    output:
    path "fastq"

    shell:
    '''
    !{params.bin_path}/dorado demux --sample-sheet !{params.samplesheet} --kit-name !{params.kitname} --no-trim  --output-dir fastq  --barcode-arrangement !{params.toml} --barcode-sequences !{params.barcode_seqs} --emit-fastq !{merged_bam_file}
    
    # Rename FASTQ files: remove everything before the first '_'
    for file in fastq/*.fastq; do
        new_name=${file#*/} # Remove the directory part
        new_name=${new_name#*_} # Remove everything before the first '_'
        mv "$file" "fastq/$new_name"
    done
    '''
}