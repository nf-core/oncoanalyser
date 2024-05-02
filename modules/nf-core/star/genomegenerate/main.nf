process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

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
    if (args_list.contains('--genomeSAindexNbases')) {
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
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    } else {
        """
        samtools faidx $fasta
        NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

        mkdir -p star_index/
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star_index/ \\
            --genomeFastaFiles $fasta \\
            $include_gtf \\
            --runThreadN $task.cpus \\
            --genomeSAindexNbases \$NUM_BASES \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }
}
