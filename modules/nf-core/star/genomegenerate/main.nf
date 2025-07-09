process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.3a--0' :
        'quay.io/biocontainers/star:2.7.3a--0' }"

    input:
    path fasta
    path gtf

    output:
    path "star_index"  , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args_list   = args.tokenize()
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''

    """
    mkdir -p star_index/

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
        """

    stub:
    if (gtf) {
        """
        mkdir -p star_index/
        touch star_index/Genome
        touch star_index/Log.out
        touch star_index/SA
        touch star_index/SAindex
        touch star_index/chrLength.txt
        touch star_index/chrName.txt
        touch star_index/chrNameLength.txt
        touch star_index/chrStart.txt
        touch star_index/exonGeTrInfo.tab
        touch star_index/exonInfo.tab
        touch star_index/geneInfo.tab
        touch star_index/genomeParameters.txt
        touch star_index/sjdbInfo.txt
        touch star_index/sjdbList.fromGTF.out.tab
        touch star_index/sjdbList.out.tab
        touch star_index/transcriptInfo.tab

        echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
        """
    } else {
        """
        mkdir -p star_index/
        touch star_index/Genome
        touch star_index/Log.out
        touch star_index/SA
        touch star_index/SAindex
        touch star_index/chrLength.txt
        touch star_index/chrName.txt
        touch star_index/chrNameLength.txt
        touch star_index/chrStart.txt
        touch star_index/genomeParameters.txt

        echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
        """
    }
}
