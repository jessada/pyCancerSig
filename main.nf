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
    nextflow run --bam_pair '/path/to/folder/with/bam/pairs/*-{normal,tumor}.bam' --somatic_snv_vcf /path/to/folder/with/somatic/snv/vcfs/*.vcf.gz --somatic_sv_vcf /path/to/folder/with/somatic/sv/vcfs/*.vcf --reference GRCh37.fa --project sens2010123
    Mandatory arguments:
      --bam_pair                    Path to input tumor-normal bam (must be surrounded with quotes)
      --somatic_snv_vcf             Path to input somatic single nucleotide variant in vcf format
      --somatic_sv_vcf              Path to input somatic structural variant in vcf format
      --reference                   Reference file used when the bams were aligned. Currently, the workflow assume that all bams were aligned using the same reference
      --project                     SLURM project code
    MSI
      --microsates                  Homopolymer and microsates file required by msisensor
      --save_msi                    Save the MSI files - not done by default
    TiTv
      --save_titv                   Save the TiTv files - not done by default
    features
      --save_features               Save the join features files - not done by default
    Other options:
      --freq_cutoff                 Cut-off criteria for variants to be included in the profile
      --out_dir                     The output directory where the results will be saved
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
params.save_titv = false
params.save_features = false
params.freq_cutoff = 0.01
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
if (params.save_titv){
    titv_out_dir = file( "${params.out_dir}/titv/" )
    if( !titv_out_dir.exists() ) {
      titv_out_dir.mkdirs()
    }
}
else {
    titv_out_dir = null
}
if (params.save_features){
    features_out_dir = file( "${params.out_dir}/features/" )
    if( !features_out_dir.exists() ) {
      features_out_dir.mkdirs()
    }
}
else {
    features_out_dir = null
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
    .ifEmpty { error "Cannot find any normal-tumor pairs in ${params.bam_pair}" }  
    .set { ch_bam_pairs }
/*
 * Create a channel for somatic SNV vcf
 */
Channel
    .fromPath("$params.somatic_snv_vcf")
    .ifEmpty { error "Cannot find any somatic SNV vcf in ${params.somatic_snv_vcf}" }  
    .map { file -> tuple(file.baseName.replaceFirst(/.vcf/, ""), file) }
    .set { ch_somatic_snv_vcfs }
/*
 * Create a channel for somatic SV vcf
 */
Channel
    .fromPath("$params.somatic_sv_vcf")
    .ifEmpty { error "Cannot find any somatic SV vcf in ${params.somatic_sv_vcf}" }  
    .map { file -> tuple(file.baseName, file) }
    .set { ch_somatic_sv_vcfs }

// log run detail
log.info "========================================="
log.info "          Cancer Signature v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']                  = custom_runName ?: workflow.runName
summary['Save MSI']                  = params.save_msi ? 'Yes' : 'No'
summary['Save Ti/Tv']                = params.save_titv ? 'Yes' : 'No'
summary['Save Features']             = params.save_features ? 'Yes' : 'No'
summary['Variant frequency cut-off'] = params.freq_cutoff
summary['Working dir']               = workflow.workDir
summary['Script dir']                = workflow.projectDir
summary['Config Profile']            = workflow.profile
summary['Current home']              = "$HOME"
summary['Current user']              = "$USER"
summary['Current path']              = "$PWD"
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
summary['Bam pair ']                 = params.bam_pair
summary['Somatic SNV vcf ']          = params.somatic_snv_vcf
summary['Somatic SV vcf ']          = params.somatic_sv_vcf
summary['Genome reference ']         = params.reference
summary['Output dir']                = params.out_dir
if(msi_out_dir) summary['MSI output dir'] = msi_out_dir
if(titv_out_dir) summary['Ti/Tv output dir'] = titv_out_dir
if(features_out_dir) summary['Ti/Tv output dir'] = features_out_dir


log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "========================================="

/*
 * STEP 2 - feature extraction
 */

/*
 * STEP 2.1.1 - scan normal-tumor bam for MSI loci
 */
process msisensor_msi {
    tag "$sample_id"
    if (params.save_msi) {
        publishDir path: msi_out_dir,
                   mode: 'copy'
    }

    input:
    set sample_id, bam_pair from ch_bam_pairs

    output:
    set sample_id, file("${sample_id}"), file("${sample_id}_somatic"), file("${sample_id}_germline"), file("${sample_id}_dis") into msi_files

    script:
    """
    $baseDir/scripts/msi/msisensor_msi.sh \
     -d $params.microsates                \
     -n ${bam_pair[0]}                    \
     -t ${bam_pair[1]}                    \
     -o $sample_id
    """
}
/*
 * STEP 2.1.2 - extract MSI features
 */
process extract_msi_features {
    tag "$sample_id"

    input:
    set sample_id, file(msi_score), file(msi_somatic), file(msi_germline), file(msi_dis) from msi_files

    output:
    set sample_id, file("${sample_id}.msi-features.txt") into msi_features

    script:
    """
    $baseDir/scripts/msi/extract_msi_features.py \
     --raw_msisensor_report $msi_score           \
     --raw_msisensor_somatic $msi_somatic        \
     --sample_id $sample_id                      \
     --output_file ${sample_id}.msi-features.txt
    """
}
/*
 * STEP 2.2.1 - scan for transition and transversion in somotic SNV vcf
 */
process count_titv {
    tag "$sample_id"
    if (params.save_titv) {
        publishDir path: titv_out_dir,
                   mode: 'copy'
    }

    input:
    set sample_id, sample_vcf from ch_somatic_snv_vcfs

    output:
    set sample_id, file("${sample_id}.event"), file("${sample_id}.titv") into titv_files

    script:
    """
    $baseDir/scripts/titv/count_titv.sh \
     -i $sample_vcf                     \
     -g IGT                             \
     -r $params.reference               \
     -f $params.freq_cutoff             \
     -e ${sample_id}.event              \
     -o ${sample_id}.titv
    """
}
/*
 * STEP 2.3.1 - scan bam file for structural variants
 */
// Temporary skip due to FindSV dependency
//process run_FindSV {
//}
/*
 * STEP 2.3.2 - perform SV somatic call
 */
// Temporary skip due to FindSV dependency
//process somatic_sv_call {
//}
/*
 * STEP 2.3.3 - extract structural variation features
 */
process extract_sv_features {
    tag "$sample_id"

    input:
    set sample_id, somatic_sv_vcf from ch_somatic_sv_vcfs

    output:
    set sample_id, file("${sample_id}.sv-features.txt") into sv_features


    script:
    """
    $baseDir/scripts/sv/extract_sv_features.py \
     --vcf ${somatic_sv_vcf}                     \
     --frequency $params.freq_cutoff           \
     --out ${sample_id}.sv-features.txt
    """
}

all_features = msi_features.join(sv_features, by:0).join(titv_files, by: 0)

/*
 * STEP 3 - merge features
 */

/*
 * STEP 3.1 - join features
 */
process join_features {
    tag "$sample_id"

    input:
    set sample_id, file(msi_features), file(sv_features), file(titv_event), file(titv_titv) from all_features

    output:
    file "${sample_id}.join-features.txt" into join_features

    script:
    """
    $baseDir/scripts/features/parse_n_join_features.py \
     --sample_id $sample_id                            \
     --msi_features $msi_features                      \
     --titv_features $titv_titv                    \
     --sv_features $sv_features                        \
     --output_filename ${sample_id}.join-features.txt
    """
}

join_features.println()

