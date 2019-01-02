#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
               C a n c e r   S i g n a t u r e   W o r k f l o w
========================================================================================

----------------------------------------------------------------------------------------
*/
 

def usage_message() {
    log.info"""
    =========================================
              Cancer Signature v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run --bam_pair '*/*-{normal,tumor}.bam' --reference GRCh37.fa --project sens2010123
    Mandatory arguments:
      --bam_pair                    Path to input tumor-normal bam (must be surrounded with quotes)
      --reference                   Reference file used when the bams were aligned. Currently, the workflow assume that all bams were aligned using the same reference
      --project                     SLURM project code
    MSI
      --microsates                  Homopolymer and microsates file required by msisensor
      --save_msi                    Save the MSI files - not done by default
    Other options:
      --out_dir                      The output directory where the results will be saved
    """.stripIndent()
}

params.help = false
if (params.help) {
    usage_message()
    exit 0
}

// Configurable variables
params.name = false
//params.project = false
params.save_msi = false
params.out_dir = './results'
if (params.save_msi){
    msi_out_dir = file( "${params.out_dir}/msi/" )
    if( !msi_out_dir.exists() ) {
      msi_out_dir.mkdirs()
    }
}
else {
    msi_out_dir = null
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for normal-tumor bams
 */

Channel
    .fromFilePairs("$params.bam_pair")
    .ifEmpty { error "Cannot find any normal-tumor paris in ${params.bam_pair}" }  
    .set { ch_bam_pairs }


// log run detail

log.info "========================================="
log.info "          Cancer Signature v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Save MSI']       = params.save_msi ? 'Yes' : 'No'
summary['Working dir']    = workflow.workDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
summary['Bam pair ']      = params.bam_pair
summary['Output dir']     = params.out_dir
if(msi_out_dir) summary['MSI output dir'] = msi_out_dir


log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * STEP 2 - feature extraction
 */

/*
 * STEP 2.1.1 - scan normal-tumor bam for MSI loci
 */
process msisensor_msi {
    if (params.save_msi) {
        publishDir path: msi_out_dir,
                   mode: 'copy'
    }

    input:
    set pair_id, bam_pair from ch_bam_pairs

    output:
    set pair_id, file("${pair_id}"), file("${pair_id}_somatic"), file("${pair_id}_germline"), file("${pair_id}_dis") into msi_files

    script:
    """
    $baseDir/scripts/msi/msisensor_msi.sh \
     -d $params.microsates                \
     -n ${bam_pair[0]}                    \
     -t ${bam_pair[1]}                    \
     -o $pair_id
    """
}
/*
 * STEP 2.1.2 - extract MSI feature
 */
process extract_msi_feature {
    if (params.save_msi) {
        publishDir path: msi_out_dir,
                   mode: 'copy'
    }

    input:
    set pair_id, file(msi_score), file(msi_somatic), file(msi_germline), file(msi_dis) from msi_files

    output:
    set pair_id, file("${pair_id}.msi-features.txt") into msi_feature

    script:
    """
    $baseDir/scripts/msi/extract_msi_features.py \
     --raw_msisensor_report $msi_score           \
     --raw_msisensor_somatic $msi_somatic        \
     --sample_id $pair_id                        \
     --output_file ${pair_id}.msi-features.txt
    """
}
