/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
import nextflow.Nextflow

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcmicro_pipeline'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

include { EXTRACTCHANNELS } from '../modules/local/extractchannels'
include { UNSTITCH        } from '../modules/local/unstitch'
include { RESTITCH        } from '../modules/local/restitch'
include { SPOTIFLOW       } from '../modules/local/spotiflow'
include { MAPS            } from '../modules/local/maps'

include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { BASICPY } from '../modules/nf-core/basicpy/main'
include { ASHLAR } from '../modules/nf-core/ashlar/main'
include { BACKSUB } from '../modules/nf-core/backsub/main'
include { CELLPOSE } from '../modules/nf-core/cellpose/main'
include { DEEPCELL_MESMER } from '../modules/nf-core/deepcell/mesmer/main'
include { MCQUANT } from '../modules/nf-core/mcquant/main'
include { SCIMAP_MCMICRO } from '../modules/nf-core/scimap/mcmicro/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    input_type = params.input_cycle ? "cycle" : "sample"

    ch_from_marker_sheet = Channel.fromSamplesheet(
        "marker_sheet",
        skip_duplicate_check: false
        )
    ch_from_marker_sheet.view()
    //
    // MODULE: BASICPY
    //

    if ( params.illumination ) {

        if (params.illumination == 'basicpy') {

            ch_samplesheet
                .transpose()
                .map { [[it[1].split('/')[-1][0..-5],it[0]], it[1]] }
                .set { ashlar_input_keyed }

                ch_samplesheet
                    .transpose()
                    .set { ch_basicpy_input }

            BASICPY(ch_basicpy_input)
            ch_versions = ch_versions.mix(BASICPY.out.versions)

            BASICPY.out.fields
                .transpose()
                .map { [[it[1].getBaseName()[0..-5],it[0]], it[1]]}
                .groupTuple()
                .set { correction_files_keyed }

            ashlar_input_keyed
                .concat(correction_files_keyed)
                .groupTuple()
                .map { [it[0][1], it[1][1]] }
                .transpose()
                .branch {
                    dfp: it =~ /-dfp.tiff/
                    ffp: it =~ /-ffp.tiff/
                }
                .set { ordered_correction_files }
            ch_dfp = ordered_correction_files.dfp
                .groupTuple()
                .map { it[1] }
            ch_ffp = ordered_correction_files.ffp
                .groupTuple()
                .map { it[1] }

        } else if(params.illumination == 'manual') {

            if (input_type == "cycle") {
                samplesheet = "input_cycle"
                dfp_index = 4
                ffp_index = 5
            } else if (input_type == "sample") {
                samplesheet = "input_sample"
                dfp_index = 3
                ffp_index = 4
            }
            ch_manual_illumination_correction = Channel.fromSamplesheet(
                samplesheet,
                skip_duplicate_check: false
            )
            .multiMap
                { it ->
                    dfp: it[dfp_index]
                    ffp: it[ffp_index]
                }

            ch_dfp = ch_manual_illumination_correction.dfp
            ch_ffp = ch_manual_illumination_correction.ffp
        }

    } else {
        ch_dfp = []
        ch_ffp = []
    }

    INPUT_CHECK( input_type, params.input_sample, params.input_cycle, params.marker_sheet )
    // ASHLAR(ch_samplesheet.ashlar_input, [], [])
    if ( !params.skip_registration ) {
        ASHLAR(ch_samplesheet, ch_dfp, ch_ffp)
        ch_versions = ch_versions.mix(ASHLAR.out.versions)
    }

    if ( params.extract_channel ){
        EXTRACTCHANNELS(ch_samplesheet, params.marker_sheet)
        ch_versions = ch_versions.mix(EXTRACTCHANNELS.out.versions)

        SPOTIFLOW(EXTRACTCHANNELS.out.extracted_channel)
        ch_versions =ch_versions.mix(SPOTIFLOW.out.versions)
    }

    if ( params.unstitch_restitch ){
        UNSTITCH(ch_samplesheet, params.marker_sheet)
        ch_versions = ch_versions.mix(EXTRACTCHANNELS.out.versions)

        UNSTITCH.out.tiles.flatMap { meta, fileList ->
            fileList.collect { file ->
                def tileId = file.toString().tokenize('_').last().tokenize('.').first()
                def newMeta = meta.clone()
                newMeta.put("tileID", tileId)
                tuple(newMeta, file)
            }
        }.set { segmentation_input }
    }

    params.unstitch_restitch ? segmentation_input.set {seg_input} : params.skip_registration ? ch_samplesheet.set { seg_input } : ASHLAR.out.tif.set { seg_input }

    // // Run Background Correction
    // BACKSUB(ASHLAR.out.tif, ch_markers)
    //BACKSUB(ASHLAR.out.tif, [[id: "backsub"], params.marker_sheet])
    //ch_versions = ch_versions.mix(BACKSUB.out.versions)

    // Run Segmentation

    DEEPCELL_MESMER(seg_input, [[:],[]])
    ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)

    if ( params.unstitch_restitch ){
        DEEPCELL_MESMER.out.mask
            .map {
                meta, tiff -> [meta.subMap("id"), tiff, meta.tileID]
            }.groupTuple()
            .map{
                meta, paths, tileids -> [meta, [paths[0], tileids[0]], [paths[1], tileids[1]], [paths[2], tileids[2]], [paths[3], tileids[3]]]
            }.map{
                meta, tile1, tile2, tile3, tile4 -> [meta, [tile1, tile2, tile3, tile4].sort{ it[1] }]
            }.map{
                meta, list -> [meta, list[0], list[1], list[2], list[3]] // sorted will have null as first list
            }.map{
                [it[0],it[1][0],it[2][0],it[3][0],it[4][0]]
            }.set { restitch_in }
        RESTITCH(restitch_in)
        ch_versions = ch_versions.mix(RESTITCH.out.versions)
    }

    params.unstitch_restitch ? RESTITCH.out.mask.set { segmentation_out } : DEEPCELL_MESMER.out.mask.set { segmentation_out }

    // Run Quantification
    mcquant_in = ch_samplesheet.join(segmentation_out).multiMap { it ->
        image: [it[0], it[1]]
        mask: [it[0], it[2]]
    }
    MCQUANT(mcquant_in.image,
            mcquant_in.mask,
            [[:], file(params.marker_sheet)])
    ch_versions = ch_versions.mix(MCQUANT.out.versions)

    MAPS(MCQUANT.out.csv,
        Channel.fromPath('/workspace/nf_mcmicro_files/best_model_20240527.pt'),
        Channel.fromPath('/workspace/nf_mcmicro_files/label_id.csv'),
        Channel.fromPath(params.marker_sheet))
    ch_versions = ch_versions.mix(MAPS.out.versions)
    /*
    // // Run Reporting
    SCIMAP_MCMICRO(MCQUANT.out.csv)
    ch_versions = ch_versions.mix(SCIMAP_MCMICRO.out.versions)
    */

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    //MULTIQC (
    //    ch_multiqc_files.collect(),
    //    ch_multiqc_config.toList(),
    //    ch_multiqc_custom_config.toList(),
    //   ch_multiqc_logo.toList()
    //)

    emit:
    multiqc_report = MCQUANT.out.csv // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
