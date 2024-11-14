process ROADIE_RECYZE {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/labsyspharm/mcmicro:roadie-2023-10-25"

    input:
    tuple val(meta), path(image), path(markers)

    output:
    tuple val(meta), path("*_segprep.tif")           , emit: seg_prep
    tuple val(meta), path("*_spotdetection.tif")     , emit: spot_detection, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    recyze.py \\
        --in ${image} \\
        --out ${prefix} \\
        --markers ${markers} \\
        --num-threads $task.cpus \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        recyze: \$(recyze.py --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_segmentation.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        recyze: \$(recyze.py --version)
    END_VERSIONS
    """
}
