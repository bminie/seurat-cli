#!/usr/bin/env nextflow

/*
 * =============================================================================
 * Seurat CLI Nextflow Pipeline
 * =============================================================================
 *
 * Flexible pipeline for running Seurat CLI scripts in sequence.
 *
 * Usage:
 *   nextflow run main.nf --input /path/to/data --outdir results
 *   nextflow run main.nf --input /path/to/data --workflow sctransform
 *   nextflow run main.nf -profile demo
 *
 * =============================================================================
 */

nextflow.enable.dsl=2

// -----------------------------------------------------------------------------
// Parameters
// -----------------------------------------------------------------------------

params.help = false

// Workflow selection
params.workflow = "sctransform"  // basic, sctransform, integration, cell_cycle, full

// Input/Output
params.input = null
params.input_list = null  // For integration: file with sample paths
params.outdir = "results"
params.scripts_dir = "${projectDir}/../../scripts"

// Common parameters
params.species = "human"
params.threads = 4

// Demo mode
params.demo = false

// Basic analysis parameters
params.basic_min_features = 200
params.basic_max_features = 2500
params.basic_max_mt = 5
params.basic_resolution = 0.5
params.basic_dims = "1:10"

// SCTransform parameters
params.sct_n_features = 3000
params.sct_resolution = 0.5
params.sct_dims = "1:30"
params.sct_vars_to_regress = ""

// Integration parameters
params.int_method = "CCAIntegration"
params.int_normalization = "LogNormalize"
params.int_resolution = 0.5
params.int_dims = "1:30"
params.int_sample_names = ""
params.int_conditions = ""

// DE parameters
params.de_mode = "all_markers"
params.de_test = "wilcox"
params.de_logfc = 0.25
params.de_min_pct = 0.1

// Cell cycle parameters
params.cc_regress = false
params.cc_regress_difference = false

// Visualization parameters
params.viz_reduction = "umap"
params.viz_features = ""

// -----------------------------------------------------------------------------
// Help Message
// -----------------------------------------------------------------------------

def helpMessage() {
    log.info """
    =============================================================================
    Seurat CLI Nextflow Pipeline
    =============================================================================

    Usage:
      nextflow run main.nf --input /path/to/data --outdir results

    Required:
      --input           Input data path (10x directory, .h5, .rds, or .csv)
                        OR use --demo for built-in dataset

    Workflow Options:
      --workflow        Workflow type: basic, sctransform, integration, cell_cycle, full
                        Default: sctransform

    Common Options:
      --outdir          Output directory (default: results)
      --species         Species: human or mouse (default: human)
      --threads         Number of threads (default: 4)
      --demo            Use built-in demo dataset

    Workflow Descriptions:
      basic:        01_basic_analysis → 04_de → 06_visualization
      sctransform:  02_sctransform → 04_de → 06_visualization
      integration:  03_integration → 04_de → 06_visualization
      cell_cycle:   02_sctransform → 05_cell_cycle → 04_de → 06_visualization
      full:         02_sctransform → 05_cell_cycle → 04_de → 06_visualization

    Script-Specific Parameters:
      See nextflow.config for all available parameters.

    Examples:
      # Run with demo data
      nextflow run main.nf --demo --outdir demo_results

      # Basic analysis
      nextflow run main.nf --input /data/sample --workflow basic

      # SCTransform with cell cycle regression
      nextflow run main.nf --input /data/sample --workflow full --cc_regress_difference

      # Integration of multiple samples
      nextflow run main.nf --input_list samples.txt --workflow integration

    =============================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.demo && !params.input && !params.input_list) {
    helpMessage()
    log.error "Error: --input or --input_list is required (or use --demo)"
    exit 1
}

// -----------------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------------

process BASIC_ANALYSIS {
    tag "basic_analysis"
    publishDir "${params.outdir}/01_basic", mode: 'copy'

    input:
    path input_data

    output:
    path "seurat_analyzed.rds", emit: rds
    path "cluster_markers_all.csv", emit: markers
    path "*.csv"
    path "*.png"
    path "*.txt"

    script:
    def demo_flag = params.demo ? "--demo" : ""
    def input_arg = params.demo ? "" : "--input ${input_data}"
    """
    Rscript ${params.scripts_dir}/01_basic_analysis.R \
        ${input_arg} \
        --output . \
        --species ${params.species} \
        --min_features ${params.basic_min_features} \
        --max_features ${params.basic_max_features} \
        --max_mt ${params.basic_max_mt} \
        --resolution ${params.basic_resolution} \
        --dims_use "${params.basic_dims}" \
        --threads ${params.threads} \
        ${demo_flag}
    """
}

process SCTRANSFORM {
    tag "sctransform"
    publishDir "${params.outdir}/01_sctransform", mode: 'copy'

    input:
    path input_data

    output:
    path "seurat_sctransform.rds", emit: rds
    path "cluster_markers_all.csv", emit: markers
    path "*.csv"
    path "*.png"
    path "*.txt"

    script:
    def demo_flag = params.demo ? "--demo" : ""
    def input_arg = params.demo ? "" : "--input ${input_data}"
    def vars_regress = params.sct_vars_to_regress ? "--vars_to_regress \"${params.sct_vars_to_regress}\"" : ""
    """
    Rscript ${params.scripts_dir}/02_sctransform.R \
        ${input_arg} \
        --output . \
        --species ${params.species} \
        --n_variable_features ${params.sct_n_features} \
        --resolution ${params.sct_resolution} \
        --dims_use "${params.sct_dims}" \
        --threads ${params.threads} \
        ${vars_regress} \
        ${demo_flag}
    """
}

process INTEGRATION {
    tag "integration"
    publishDir "${params.outdir}/01_integration", mode: 'copy'

    input:
    path input_list

    output:
    path "seurat_integrated.rds", emit: rds
    path "cluster_markers_all.csv", emit: markers
    path "*.csv"
    path "*.png"
    path "*.txt"

    script:
    def demo_flag = params.demo ? "--demo" : ""
    def input_arg = params.demo ? "" : "--input_list ${input_list}"
    def sample_names = params.int_sample_names ? "--sample_names \"${params.int_sample_names}\"" : ""
    def conditions = params.int_conditions ? "--conditions \"${params.int_conditions}\"" : ""
    """
    Rscript ${params.scripts_dir}/03_integration.R \
        ${input_arg} \
        --output . \
        --species ${params.species} \
        --method ${params.int_method} \
        --normalization ${params.int_normalization} \
        --resolution ${params.int_resolution} \
        --dims_use "${params.int_dims}" \
        --threads ${params.threads} \
        ${sample_names} \
        ${conditions} \
        ${demo_flag}
    """
}

process CELL_CYCLE {
    tag "cell_cycle"
    publishDir "${params.outdir}/02_cell_cycle", mode: 'copy'

    input:
    path seurat_rds

    output:
    path "seurat_cell_cycle.rds", emit: rds
    path "cell_cycle_scores.csv", emit: scores
    path "*.csv"
    path "*.png"
    path "*.txt"

    script:
    def regress_flag = params.cc_regress ? "--regress_cc" : ""
    def regress_diff_flag = params.cc_regress_difference ? "--regress_cc_difference" : ""
    """
    Rscript ${params.scripts_dir}/05_cell_cycle.R \
        --input ${seurat_rds} \
        --output . \
        --species ${params.species} \
        ${regress_flag} \
        ${regress_diff_flag}
    """
}

process DIFFERENTIAL_EXPRESSION {
    tag "differential_expression"
    publishDir "${params.outdir}/03_differential_expression", mode: 'copy'

    input:
    path seurat_rds

    output:
    path "de_results_all.csv", emit: results
    path "de_results_significant.csv", emit: significant
    path "*.csv"
    path "*.png", optional: true
    path "*.txt", optional: true

    script:
    """
    Rscript ${params.scripts_dir}/04_differential_expression.R \
        --input ${seurat_rds} \
        --output . \
        --mode ${params.de_mode} \
        --test_use ${params.de_test} \
        --logfc_threshold ${params.de_logfc} \
        --min_pct ${params.de_min_pct} \
        --threads ${params.threads}
    """
}

process VISUALIZATION {
    tag "visualization"
    publishDir "${params.outdir}/visualization", mode: 'copy'

    input:
    path seurat_rds
    path marker_file

    output:
    path "visualization_summary.md", emit: summary
    path "*.png"
    path "*.csv"
    path "*.txt"

    script:
    def features_arg = params.viz_features ? "--features \"${params.viz_features}\"" : ""
    """
    Rscript ${params.scripts_dir}/06_visualization.R \
        --input ${seurat_rds} \
        --output . \
        --marker_file ${marker_file} \
        --reduction ${params.viz_reduction} \
        ${features_arg}
    """
}

// -----------------------------------------------------------------------------
// Workflows
// -----------------------------------------------------------------------------

workflow BASIC_WORKFLOW {
    take:
    input_ch

    main:
    BASIC_ANALYSIS(input_ch)
    DIFFERENTIAL_EXPRESSION(BASIC_ANALYSIS.out.rds)
    VISUALIZATION(BASIC_ANALYSIS.out.rds, DIFFERENTIAL_EXPRESSION.out.significant)

    emit:
    rds = BASIC_ANALYSIS.out.rds
    de = DIFFERENTIAL_EXPRESSION.out.results
    viz = VISUALIZATION.out.summary
}

workflow SCTRANSFORM_WORKFLOW {
    take:
    input_ch

    main:
    SCTRANSFORM(input_ch)
    DIFFERENTIAL_EXPRESSION(SCTRANSFORM.out.rds)
    VISUALIZATION(SCTRANSFORM.out.rds, DIFFERENTIAL_EXPRESSION.out.significant)

    emit:
    rds = SCTRANSFORM.out.rds
    de = DIFFERENTIAL_EXPRESSION.out.results
    viz = VISUALIZATION.out.summary
}

workflow INTEGRATION_WORKFLOW {
    take:
    input_ch

    main:
    INTEGRATION(input_ch)
    DIFFERENTIAL_EXPRESSION(INTEGRATION.out.rds)
    VISUALIZATION(INTEGRATION.out.rds, DIFFERENTIAL_EXPRESSION.out.significant)

    emit:
    rds = INTEGRATION.out.rds
    de = DIFFERENTIAL_EXPRESSION.out.results
    viz = VISUALIZATION.out.summary
}

workflow CELLCYCLE_WORKFLOW {
    take:
    input_ch

    main:
    SCTRANSFORM(input_ch)
    CELL_CYCLE(SCTRANSFORM.out.rds)
    DIFFERENTIAL_EXPRESSION(CELL_CYCLE.out.rds)
    VISUALIZATION(CELL_CYCLE.out.rds, DIFFERENTIAL_EXPRESSION.out.significant)

    emit:
    rds = CELL_CYCLE.out.rds
    de = DIFFERENTIAL_EXPRESSION.out.results
    viz = VISUALIZATION.out.summary
}

workflow FULL_WORKFLOW {
    take:
    input_ch

    main:
    SCTRANSFORM(input_ch)
    CELL_CYCLE(SCTRANSFORM.out.rds)
    DIFFERENTIAL_EXPRESSION(CELL_CYCLE.out.rds)
    VISUALIZATION(CELL_CYCLE.out.rds, DIFFERENTIAL_EXPRESSION.out.significant)

    emit:
    rds = CELL_CYCLE.out.rds
    de = DIFFERENTIAL_EXPRESSION.out.results
    viz = VISUALIZATION.out.summary
}

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------

workflow {
    // Create input channel
    if (params.demo) {
        input_ch = Channel.of(file("${projectDir}/dummy_input"))
    } else if (params.workflow == "integration") {
        input_ch = Channel.fromPath(params.input_list)
    } else {
        input_ch = Channel.fromPath(params.input)
    }

    // Run selected workflow
    if (params.workflow == "basic") {
        BASIC_WORKFLOW(input_ch)
    } else if (params.workflow == "sctransform") {
        SCTRANSFORM_WORKFLOW(input_ch)
    } else if (params.workflow == "integration") {
        INTEGRATION_WORKFLOW(input_ch)
    } else if (params.workflow == "cell_cycle") {
        CELLCYCLE_WORKFLOW(input_ch)
    } else if (params.workflow == "full") {
        FULL_WORKFLOW(input_ch)
    } else {
        log.error "Unknown workflow: ${params.workflow}"
        exit 1
    }
}
