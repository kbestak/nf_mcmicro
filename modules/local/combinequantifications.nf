process COMBINEQUANTIFICATIONS {
    tag   "$meta.id"
    label 'process_single'

    container 'docker.io/kbestak/unstitch_restitch:0.0.1'

    input:
    tuple val(meta), path(spot2cell), path(mcquant)

    output:
    tuple val(meta), path("*.csv"), emit: combined_quants
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION='0.0.1'
    """
    combine_quantifications.py \\
        --spot2cell ${spot2cell} \\
        --mcquant ${mcquant} \\
        --output ${prefix}.csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spot2cell: $VERSION
    END_VERSIONS
    """
}
