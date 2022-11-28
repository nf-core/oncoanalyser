process INDEX {
    container 'docker.io/scwatts/gridss:2.13.2--3'

    input:
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_bwa_index_dir, stageAs: 'bwa_index'
    path genome_bwa_index_image
    val indices

    output:
    path '*.dict'       , emit: dict, optional: true
    path '*.img'        , emit: img, optional: true
    path '*.gridsscache', emit: index, optional: true
    path 'versions.yml' , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sequence_dict_arg = indices.contains('samtools_dict') ? 'true' : 'false'
    def bwa_index_image_arg = indices.contains('bwa_index_image') ? 'true' : 'false'
    def gridss_index_arg = indices.contains('gridss_index') ? 'true' : 'false'

    """
    # Symlink BWA indices next to assembly FASTA
    ln -s \$(find -L ${genome_bwa_index_dir} -type f) ./

    # Run
    java -Xmx${task.memory.giga}g \\
      -XX:ParallelGCThreads=${task.cpus} \\
      -Dsamjdk.reference_fasta=${genome_fasta} \\
      -Dsamjdk.use_async_io_read_samtools=true \\
      -Dsamjdk.use_async_io_write_samtools=true \\
      -Dsamjdk.use_async_io_write_tribble=true \\
      -Dsamjdk.buffer_size=4194304 \\
      -Dsamjdk.async_io_read_threads=${task.cpus} \\
      -cp ${task.ext.jarPath} \\
        REFERENCE_SEQUENCE=${genome_fasta} \\
        CREATE_SEQUENCE_DICTIONARY=${sequence_dict_arg} \\
        CREATE_BWA_INDEX_IMAGE=${bwa_index_image_arg} \\
        CREATE_GRIDSS_REFERENCE_CACHE=${gridss_index_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(java -cp "${task.ext.jarPath}" gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
    END_VERSIONS
    """

    stub:
    """
    touch ${genome_fasta.name}.dict
    touch ${genome_fasta.name}.img
    touch ${genome_fasta.name}.gridsscache
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

