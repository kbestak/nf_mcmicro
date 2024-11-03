process MCMICRO_ROADIE {
    tag "$meta.id"
    label 'process_single'

    conda "numpy,ome_types,pandas,tifffile,zarr"
    container "community.wave.seqera.io/library/numpy_ome-types_pandas_tifffile_zarr:e51681255ba835b1"
    //container "ghcr.io/labsyspharm/mcmicro:roadie-2023-10-25"

    input:
    tuple val(meta), path(image), path(markers)

    output:
    tuple val(meta), path("*_extracted_segmentation_data.tif") , emit: segmentation
    tuple val(meta), path("*_extracted_spot_channel.tif")      , emit: spots
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    roadie.py \\
        --in ${image} \\
        --markers ${markers} \\
        -o ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcmicro_roadie: \$(roadie.py --version)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_extracted_segmentation_data.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcmicro_roadie: \$(roadie.py --version)
    END_VERSIONS
    """
}
