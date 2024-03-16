process CHORD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-chord:2.03--r43hdfd78af_0' :
        'quay.io/biocontainers/r-chord:2.03--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(smlv_vcf), path(sv_vcf)
    val genome_ver

    output:
    tuple val(meta), path('chord/'), emit: chord_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env Rscript
    library('CHORD')

    sampleName <- '${meta.sample_id}'
    snvIndVcf <- '${smlv_vcf}'
    svVcf <- '${sv_vcf}'
    refGenomeVsn <- '${genome_ver}'

    sigOutTxt <- 'chord/${meta.sample_id}_chord_signatures.txt'
    prdOutTxt <- 'chord/${meta.sample_id}_chord_prediction.txt'

    dir.create('chord/')

    if (refGenomeVsn == '37') {
        library(BSgenome.Hsapiens.UCSC.hg19)
        refGenome <- BSgenome.Hsapiens.UCSC.hg19
    } else if (refGenomeVsn == '38') {
        library(BSgenome.Hsapiens.UCSC.hg38)
        refGenome <- BSgenome.Hsapiens.UCSC.hg38
    } else {
        stop('Unsupported ref genome version: ', refGenomeVsn, ' (should be 37 or 38)\\n')
    }

    cat('[INFO] Performing chord signature extraction\\n')
    signatures <- CHORD::extractSigsChord(
        vcf.snv=snvIndVcf,
        vcf.indel=snvIndVcf,
        vcf.sv=svVcf,
        sample.name=sampleName,
        sv.caller='gridss',
        vcf.filters=list(snv='PASS', indel='PASS', sv='PASS'),
        ref.genome=refGenome
    )

    cat('[INFO] Performing chord HRD prediction\\n')
    prediction <- chordPredict(
        signatures,
        hrd.cutoff=0.5
    )

    cat('[INFO] Writing output file:', sigOutTxt,'\\n')
    write.table(signatures, file=sigOutTxt, sep='\\t')

    cat('[INFO] Writing output file:', prdOutTxt,'\\n')
    write.table(prediction, file=prdOutTxt, sep='\\t', quote=FALSE, row.names=FALSE)

    cat('[INFO] FINISHED CHORD signature extraction and HRD prediction\\n')

    sink('versions.yml')
    writeLines('"${task.process}":')
    writeLines(paste('    CHORD:', packageVersion('CHORD')))
    writeLines(paste('    mutSigExtractor:', packageVersion('mutSigExtractor')))
    sink()
    """

    stub:
    """
    mkdir -p chord/
    touch chord/${meta.sample_id}_chord_signatures.txt
    touch chord/${meta.sample_id}_chord_prediction.txt
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
