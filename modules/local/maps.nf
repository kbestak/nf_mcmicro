process MAPS {
    tag   "$meta.id"
    label 'process_medium'

    container "docker.io/kbestak/maps:1.0.0_cv1"

    input:
    tuple val(meta), path(csv)
    path(model)
    path(labelID)
    path(marker_sheet)

    output:
    tuple val(meta), path("*_pred.csv"), emit: phenotypes
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0.0"
    """
    maps_cli.py \\
        --input_csv ${csv} \\
        --output ./${prefix}_pred.csv \\
        --pretrained_model ${model} \\
        --labelID ${labelID}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maps: $VERSION
    END_VERSIONS
    """
}
