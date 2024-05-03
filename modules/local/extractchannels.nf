process EXTRACTCHANNELS {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/kbestak/unstitch_restitch:0.0.1'

    input:
    tuple val(meta), path(image)
    path(marker_sheet)

    output:
    tuple val(meta), path("*.tif"), emit: extracted_channel
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_channels.py \\
        --input ${image} \\
        --output . \\
        --markers ${marker_sheet} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extractchannels: \$(extract_channels.py --version)
    END_VERSIONS
    """
}
