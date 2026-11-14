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
params.workflow = "sctransform"  // basic, sctransform, integration, cell_cycle, full, sample_level, spatial

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

// Sample-level (scSLIDE) parameters; defaults match 03_integration.R metadata
params.sl_input_rds = null  // existing Seurat RDS instead of running integration
params.sl_sample_col = "sample"
params.sl_response = "condition"
params.sl_celltype_col = "seurat_clusters"
params.sl_reduction = "integrated"
params.sl_dims = "1:30"
params.sl_traj_dcs = "1"
params.sl_latent_vars = ""

// Spatial (Xenium/Atera) parameters
params.sp_platform = "auto"
params.sp_molecules = ""
params.sp_crop = ""
params.sp_normalization = "SCT"
params.sp_resolution = 0.3

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
      --input           Input data path (10x directory, .h5, .rds, or .csv;
                        a Xenium/Atera output folder for --workflow spatial)
                        OR use --demo for built-in dataset

    Workflow Options:
      --workflow        Workflow type: basic, sctransform, integration, cell_cycle,
                        full, sample_level, spatial
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
      sample_level: 03_integration → 07_sample_level (or 07 alone on --sl_input_rds)
      spatial:      08_spatial on --input (no demo dataset)

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

      # Sample-level analysis of integrated samples (needs >= 4 samples)
      nextflow run main.nf --input_list samples.txt --int_conditions "c,c,t,t" \
          --workflow sample_level

      # Xenium or Atera data
      nextflow run main.nf --input xenium_outs --workflow spatial --sp_molecules Gad1,Sst

    HPC (profiles: slurm or pbs, plus singularity or apptainer):
      --queue           Partition/queue (default: scheduler default)
      --cluster_options Extra submit options; write it with "=" so the value
                        isn't read as a flag: --cluster_options="--account=lab"
      --container       Image, or a local .sif file for offline compute nodes
                        (default: ghcr.io/bminie/seurat-cli:latest)
      --max_memory, --max_cpus, --max_time
                        Per-task caps (default: 16.GB, 8, 24.h)

      nextflow run main.nf -profile slurm,apptainer --input_list samples.txt \\
          --workflow integration --cluster_options="--account=lab"

    =============================================================================
    """.stripIndent()
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
    path "seurat_object.rds", emit: rds
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
    // The samples from --input_list, staged into the task directory so they
    // are readable inside containers and on any cluster node
    path samples, stageAs: 'sample*/*'

    output:
    path "seurat_integrated.rds", emit: rds
    path "cluster_markers_all.csv", emit: markers
    path "*.csv"
    path "*.png"
    path "*.txt"

    script:
    def demo_flag = params.demo ? "--demo" : ""
    def input_arg = params.demo ? "" : "--input_list samples.txt"
    def sample_paths = (samples instanceof List ? samples : [samples]).collect { p -> "'${p}'" }
    def write_list = params.demo ? "" : "printf '%s\\n' ${sample_paths.join(' ')} > samples.txt"
    def sample_names = params.int_sample_names ? "--sample_names \"${params.int_sample_names}\"" : ""
    def conditions = params.int_conditions ? "--conditions \"${params.int_conditions}\"" : ""
    """
    ${write_list}
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

process SAMPLE_LEVEL {
    tag "sample_level"
    publishDir "${params.outdir}/sample_level", mode: 'copy'

    input:
    path seurat_rds

    output:
    path "traj_de_summary.csv", emit: summary
    path "sample_embedding.csv", emit: embedding
    path "sample_object.rds"
    path "*.csv"
    path "*.png", optional: true
    path "*.txt"

    script:
    // In demo mode the script loads ifnb with donor IDs and sets its own columns
    def input_arg = params.demo ? "--demo" : "--input ${seurat_rds}"
    def celltype_arg = params.sl_celltype_col ? "--celltype_col \"${params.sl_celltype_col}\"" : ""
    def latent_arg = params.sl_latent_vars ? "--latent_vars \"${params.sl_latent_vars}\"" : ""
    """
    Rscript ${params.scripts_dir}/07_sample_level.R \
        ${input_arg} \
        --output . \
        --sample_col "${params.sl_sample_col}" \
        --response "${params.sl_response}" \
        ${celltype_arg} \
        --reduction ${params.sl_reduction} \
        --dims "${params.sl_dims}" \
        --traj_dcs "${params.sl_traj_dcs}" \
        ${latent_arg} \
        --threads ${params.threads}
    """
}

process SPATIAL {
    tag "spatial"
    publishDir "${params.outdir}/spatial", mode: 'copy'

    input:
    path data_dir

    output:
    path "cluster_assignments.csv", emit: clusters
    path "seurat_spatial.rds", emit: rds
    path "*.csv"
    path "*.png"
    path "*.txt"

    script:
    def molecules_arg = params.sp_molecules ? "--molecules \"${params.sp_molecules}\"" : ""
    def crop_arg = params.sp_crop ? "--crop \"${params.sp_crop}\"" : ""
    """
    Rscript ${params.scripts_dir}/08_spatial.R \
        --input ${data_dir} \
        --output . \
        --platform ${params.sp_platform} \
        --normalization ${params.sp_normalization} \
        --resolution ${params.sp_resolution} \
        ${molecules_arg} \
        ${crop_arg} \
        --threads ${params.threads}
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

workflow SAMPLE_LEVEL_WORKFLOW {
    take:
    input_ch

    main:
    // The demo and --sl_input_rds start from an existing object; otherwise
    // integrate the samples first
    def rds_ch = input_ch
    if (!params.demo && !params.sl_input_rds) {
        INTEGRATION(input_ch)
        rds_ch = INTEGRATION.out.rds
    }
    SAMPLE_LEVEL(rds_ch)

    emit:
    summary = SAMPLE_LEVEL.out.summary
}

workflow SPATIAL_WORKFLOW {
    take:
    input_ch

    main:
    SPATIAL(input_ch)

    emit:
    clusters = SPATIAL.out.clusters
}

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate inputs
    if (!params.demo && !params.input && !params.input_list && !params.sl_input_rds &&
        params.workflow != "spatial") {
        helpMessage()
        log.error "Error: --input or --input_list is required (or use --demo)"
        exit 1
    }
    // The spatial workflow has no demo dataset, so it ignores --demo (which
    // the demo and test profiles set) and always reads --input
    if (params.workflow == "spatial" && !params.input) {
        log.error "Error: the spatial workflow has no demo dataset; set --input to a Xenium or Atera output folder"
        exit 1
    }

    // "--cluster_options --account=x" (no "=") leaves the option without a
    // value, which Nextflow sets to true (or "true")
    if (params.cluster_options in [true, 'true']) {
        log.error 'Error: --cluster_options has no value. Write it with "=", e.g. --cluster_options="--account=lab"'
        exit 1
    }

    // Create input channel
    if (params.workflow == "spatial") {
        input_ch = channel.of(file(params.input, checkIfExists: true))
    } else if (params.demo) {
        input_ch = channel.of(file("${projectDir}/dummy_input"))
    } else if (params.workflow == "sample_level" && params.sl_input_rds) {
        input_ch = channel.fromPath(params.sl_input_rds, checkIfExists: true)
    } else if (params.workflow in ["integration", "sample_level"]) {
        if (!params.input_list) {
            log.error "Error: the ${params.workflow} workflow needs --input_list (or --demo)"
            exit 1
        }
        // One path per line (absolute, or relative to the launch directory);
        // blank lines and # comments are skipped. All samples go to one task.
        input_ch = channel.fromPath(params.input_list, checkIfExists: true)
            .splitText()
            .map { line -> line.trim() }
            .filter { line -> line && !line.startsWith('#') }
            .map { line -> file(line, checkIfExists: true) }
            .collect()
    } else {
        input_ch = channel.fromPath(params.input)
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
    } else if (params.workflow == "sample_level") {
        SAMPLE_LEVEL_WORKFLOW(input_ch)
    } else if (params.workflow == "spatial") {
        SPATIAL_WORKFLOW(input_ch)
    } else {
        log.error "Unknown workflow: ${params.workflow}"
        exit 1
    }
}
