process UNSTITCH {
    tag "$meta.id"
    label 'process_single'

    container  'docker.io/kbestak/unstitch_restitch:0.0.1'

    input:
    tuple val(meta), path(image)
    path(marker_sheet)

    output:
    tuple val(meta), path("*.tif"), emit: tiles
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    unstitch_2x2.py \\
        --input ${image} \\
        --output . \\
        --markers ${marker_sheet} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unstitch: \$(unstitch_2x2.py --version)
    END_VERSIONS
    """
}
