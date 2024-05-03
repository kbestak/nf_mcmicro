process RESTITCH {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/kbestak/unstitch_restitch:0.0.1'

    input:
    tuple val(meta), path(tile0), path(tile1), path(tile2), path(tile3)

    output:
    tuple val(meta), path("*.tif"), emit: mask
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stitch_maskgrid_2x2.py \\
        --tile_0 ${tile0} \\
        --tile_1 ${tile1} \\
        --tile_2 ${tile2} \\
        --tile_3 ${tile3} \\
        --output ${prefix}.tif \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        restitch: \$(stitch_maskgrid_2x2 --version)
    END_VERSIONS
    """
}
