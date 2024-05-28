process MCQUANT {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "docker.io/kbestak/mcquant:1.5.4_cv6"

    input:
    tuple val(meta), path(image)
    tuple val(meta2), path(mask)
    tuple val(meta3), path(markerfile)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.5.4_cv6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    python /app/CommandSingleCellExtraction.py \\
        --masks $mask \\
        --image $image \\
        --channel_names $markerfile \\
        --output . \\
        --output_filename ${prefix}.csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcquant: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.5.4_cv4'
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcquant: $VERSION
    END_VERSIONS
    """
}
