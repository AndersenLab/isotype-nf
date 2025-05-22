#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Isotype Calling Pipeline
    Authors:
    - Mike Sauria <mike.sauria@jhu.edu>
*/

include { FIND_COVERAGES    } from "./modules/local/find_coverages/main"
include { ENCODE_VCF        } from "./modules/python/encode_vcf/main"
include { GTCHECK           } from "./modules/python/gtcheck/main"
include { CALL_ISOTYPES     } from "./modules/python/call_isotypes/main"
include { ADD_UNSEQ_STRAINS } from "./modules/local/add_unseq_strains/main"

// Needed to publish results
nextflow.preview.output = true

date = new Date().format( 'yyyyMMdd' )

// Debug
if (params.debug) {
    bam_folder = "${workflow.projectDir}/test_data/bams"
    vcf_file = "${workflow.projectDir}/test_data/hard-filtered.vcf.gz"
    previous_isotypes = "${workflow.projectDir}/test_data/old_isotype_groups.tsv"
    species = "c_elegans"
    cutoff = 0.9997
}

if (params.help == false & params.debug == false) {
    if (params.vcf_file == null) {
        println """
        Please specify a vcf_file with the option --vcf_file
        """
        exit 1
    } else {
        vcf_file = params.vcf_file
    }
    if (params.species == null) {
        if (params.bam_folder == null | params.previous_isotypes == null | params.cutoff) {
            println """
            Please specify a species with option --species or bam folder, previous isotype file, and concordance cutoff with --bam_folder, --previous_isotypes, and --cutoff
            """
            exit 1
        } else {
            bam_folder = params.bam_folder
            previous_isotypes = params.previous_isotypes
            cutoff = params.cutoff
        }
        species = params.species
    } else {
        species = params.species
        if (species == "c_elegans") {
            unsequenced_strains = "${workflow.projectDir}/data/c_elegans_unsequenced_strains.tsv"
        }
        if (params.bam_folder != null) {
                bam_folder = params.bam_folder
        } else if (species == "c_elegans" | species == "c_briggsae" | species == "c_tropicalis") {
            bam_folder = "${params.dataDir}/${species}/WI/alignments/"
        } else {
            println """
            When using a species other than c_elegans, c_briggsae, or c_tropicalis,
            a bam folder must be specified with the option --bam_folder
            """
            exit 1
        }
        if (params.previous_isotypes != null) {
                previous_isotypes = params.previous_isotypes
        } else if (species == "c_elegans" | species == "c_briggsae" | species == "c_tropicalis") {
            if (species == "c_elegans") {
                previous_isotypes = "${params.dataDir}/${species}/WI/concordance/20231213/isotype_groups.tsv"
            } else if (species == "c_briggsae") {
                previous_isotypes = "${params.dataDir}/${species}/WI/concordance/20240129/isotype_groups.tsv"
            } else {
                previous_isotypes = "${params.dataDir}/${species}/WI/concordance/20231201/isotype_groups.tsv"
            }
        } else {
            println """
            When using a species other than c_elegans, c_briggsae, or c_tropicalis,
            a previous isotypes file must be specified with the option --previous_isotypes
            """
            exit 1
        }
        if (params.cutoff != null) {
                cutoff = params.cutoff
        } else if (species == "c_elegans" | species == "c_briggsae" | species == "c_tropicalis") {
            if (species == "c_elegans") {
                cutoff = 0.9997
            } else if (species == "c_briggsae") {
                cutoff = 0.9997
            } else {
                cutoff = 0.999915
            }
        } else {
            println """
            When using a species other than c_elegans, c_briggsae, or c_tropicalis,
            a concordance cutoff must be specified with the option --cutoff
            """
            exit 1
        }
    }
} else if (params.help == true & params.debug == false) {
    species = params.species
    vcf_file = params.vcf_file
    bam_folder = params.bam_folder
    previous_isotypes = params.previous_isotypes
    cutoff = params.cutoff
}

def log_summary() {

    out =  '''
_____   ______     ___    _________  ____  ____  _______  ________       ____  _____  ________  
|_   _|.' ____ \\  .'   `. |  _   _  ||_  _||_  _||_   __ \\|_   __  |     |_   \\|_   _||_   __  | 
  | |  | (___ \\_|/  .-.  \\|_/ | | \\_|  \\ \\  / /    | |__) | | |_ \\_|______ |   \\ | |    | |_ \\_| 
  | |   _.____`. | |   | |    | |       \\ \\/ /     |  ___/  |  _| _|______|| |\\ \\| |    |  _|    
 _| |_ | \\____) |\\  `-'  /   _| |_      _|  |_    _| |_    _| |__/ |      _| |_\\   |_  _| |_     
|_____| \\______.' `.___.'   |_____|    |______|  |_____|  |________|     |_____|\\____||_____|    
                                                                                               
                                              
'''

out += """

To run the pipeline:

nextflow main.nf --help
nextflow main.nf --debug
nextflow main.nf --vcf_file=/path/to/vcf_file --species c_elegans -output-dir=/path/to/output
nextflow main.nf --vcf_file=/path/to/vcf_file --bam_location=/path/to/bams --previous_isotypes=/path/to/previous_isotype_file -output-dir=/path/to/output

    parameters                 description                           Set/Default
    ==========                 ===========                           ========================
    --debug                    Use --debug to indicate debug mode    ${params.debug}
    --vcf_file                 All strains VCF file                  ${vcf_file}

    --species                  Species to call isotypes from         ${species}
    and / or
    --cutoff                   Concordance cutoff for isotype calls  ${cutoff}
    --bam_location             Directory of BAM files                ${bam_folder}
    --previous_isotypes        File containing previous isotypes     ${previous_isotypes}

    username                                                         ${"whoami".execute().in.text}

---
"""
out
}

log.info(log_summary())

if (params.help == true) {
    exit 1
}

now = new Date()
timestamp = now.format("yyyyMMdd-HH-mm-ss")
log.info("Started running ${now}")

workflow {
    main:
    ch_versions = Channel.empty()

    // Open vcf file
    ch_vcf = Channel.fromPath(vcf_file, checkIfExists: true)
        .ifEmpty { exit 1, "vcf file not found" }
    ENCODE_VCF( ch_vcf,
                Channel.fromPath("${workflow.ProjectDir}/bin/encode_vcf.py") )

    // Perform gtcheck
    GTCHECK( ENCODE_VCF.out.encoded,
             Channel.fromPath("${workflow.ProjectDir}/bin/gtcheck.py") )

    // Find coverages
    ch_bams = Channel.fromPath(bam_folder, checkIfExists: true)
        .ifEmpty { exit 1, "bam folder not found" }
    FIND_COVERAGES( ch_bams )

    // Assign isotypes
    CALL_ISOTYPES( GTCHECK.out.gtcheck,
                   Channel.fromPath("${workflow.ProjectDir}/bin/call_isotypes.py"),
                   Channel.fromPath("${workflow.ProjectDir}/bin/compare_isotype_calls.py"),
                   Channel.fromPath(previous_isotypes),
                   FIND_COVERAGES.out.coverages,
                   cutoff )

    if (species == "c_elegans") {
        // If finding isotypes for c_elegans, need to add in unsequenced historical samples
        ADD_UNSEQ_STRAINS( CALL_ISOTYPES.out.groups,
                           Channel.fromPath( unsequenced_strains) )
        ch_all_groups = ADD_UNSEQ_STRAINS.out.groups
    } else {
        ch_all_groups = Channel.empty()
    }

    // Collate and save software versions
    ch_versions
        .collectFile(name: 'workflow_software_versions.txt', sort: true, newLine: true)
        .set { ch_collated_versions }

    publish:
    ENCODE_VCF.out.encoded        >> "."
    CALL_ISOTYPES.out.groups      >> "."
    CALL_ISOTYPES.out.comparison  >> "."
    CALL_ISOTYPES.out.summary     >> "."
    GTCHECK.out.gtcheck           >> "."
    ch_all_groups                 >> "."
    ch_collated_versions          >> "."
}

// Current bug that publish doesn't work without an output closure
output {
    "." {
        mode "copy"
    }
}