process SPOT2CELL {
    tag   "$meta.id"
    label 'process_single'

    container 'docker.io/kbestak/unstitch_restitch:0.0.1'

    input:
    tuple val(meta), path(spots), path(seg_mask)

    output:
    tuple val(meta), path("*.csv"), emit: matched_spots
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION='0.0.1'
    """
    spot2cell.py \\
        --spots ${spots} \\
        --mask ${seg_mask} \\
        --output ${prefix}.csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spot2cell: $VERSION
    END_VERSIONS
    """
}
