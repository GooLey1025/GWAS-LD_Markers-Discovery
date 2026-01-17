nextflow.enable.dsl = 2

params.vcf = '705rice.graph.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz'
params.snp_vcf = '705rice.graph.0.5_0.05.full.snp.impute.biallelic.id.vcf.gz'  // SNP VCF for kinship calculation

// params.vcf = '705rice_graph.all.0.5_0.01.impute.vcf'
params.outdir = '705rice_graph.0.5_0.05.SNP_INDEL_SV'
params.phenotypes_dir = 'phenotypes/705rice'
params.assign_id = 'assign_id.sh'
params.split_vcf = false  // Set to true if VCF needs to be split (default: false, assuming VCF is already split (bcftools norm -m -both))

params.plink='plink'
params.gemma='gemma'

// general parameters
params.maf = 0.05
params.miss_threshold = 0.5

// GEMMA parameters
params.gemma_lmm_type = 1
params.hwe_threshold = 0

// PLINK parameters
// plink clump p1 is set by 1/independent tests via plink --indep-pairwise
params.clump_p2 = 0.05
params.clump_r2 = 0.1
params.clump_kb = 1000
params.plink_indep_r2 = 0.2

// params.binary_phenotypes = "LeafSheathColor_2012HZ,LeafSheathColor_2012YF"
workflow {
    // Create channels for VCF and phenotypes
    Channel.fromPath(params.vcf).set { vcf_ch }
    Channel.fromPath(params.snp_vcf).set { snp_vcf_ch }
    Channel.fromPath("${params.phenotypes_dir}/*.tsv").set { phenotypes_ch }

    plink_indep_ch = PLINK_INDEP_RUN(vcf_ch)
    threshold_ch = plink_indep_ch.threshold_file.map { it.text.trim()}

    // Optionally split VCF if needed
    if ( params.split_vcf ) {
        vcf_processed = SPLIT_VCF(vcf_ch)
    } else {
        vcf_processed = vcf_ch
    }

    // Step 1: Convert SNP VCF to PLINK format (without phenotype) for kinship calculation
    snp_plink = SNP_VCF_TO_PLINK(snp_vcf_ch)

    // Step 2: Calculate kinship matrix ONCE using SNP data (shared by all phenotypes)
    kinship = GEMMA_KINSHIP(snp_plink)

    // Step 3: Combine VCF with each phenotype for per-phenotype processing
    vcf_pheno_ch = vcf_processed.combine(phenotypes_ch)

    // Step 4: Per-phenotype PLINK conversion (with phenotype filtering)
    plink_files = PLINK_VCF_CONVERSION(vcf_pheno_ch)

    // Step 5: Combine each phenotype's PLINK files with the shared kinship
    plink_with_kinship = plink_files.combine(kinship)

    // Step 6: Per-phenotype GEMMA LMM association
    association = GEMMA_LMM_ASSOCIATION(plink_with_kinship)

    // Step 7: Generate Manhattan plots (with threshold)
    association_with_threshold = association.combine(threshold_ch)
    MANHATTAN_PLOT(association_with_threshold)

    // Step 8: Perform PLINK clumping
    clump_ch = association.combine(threshold_ch)
    PLINK_CLUMPING(clump_ch)
    
    // Step 9: Extract lead markers from all clumped files
    // Collect all clumped files (they are tuples of phenotype_name and clumped_file)
    clumped_collected = PLINK_CLUMPING.out
        .map { phenotype_name, clumped_file -> clumped_file }
        .collect()
    EXTRACT_LEAD_MARKERS(clumped_collected)
    GET_UNIQUE_MARKERS(EXTRACT_LEAD_MARKERS.out.markers_tsv)
    FILTER_VCF_BY_MARKERS(vcf_ch, GET_UNIQUE_MARKERS.out.markers_list)
}

process PLINK_INDEP_RUN {
    tag "plink_indep_run_${params.plink_indep_r2}_${params.maf}"
    publishDir "${params.outdir}/plink_indep", mode: 'copy', pattern: '*.{log,txt}'

    cpus 32
    memory '200 GB'

    input:
    path vcf

    output:
    path("threshold.txt"), emit: threshold_file
    path("${vcf.simpleName}.log"), emit: plink_indep_log
    script:
    """
    set -euo pipefail
    
    ${params.plink} --threads ${task.cpus} --vcf ${vcf} --double-id --make-bed --out effective_number \\
        --maf ${params.maf} --indep-pairwise 50 50 ${params.plink_indep_r2} \\
        --allow-extra-chr --allow-no-sex 2>&1 | tee ${vcf.simpleName}.log
    
    # Count independent markers and calculate threshold (1/n)
    if [ ! -f effective_number.prune.in ]; then
        echo "[ERROR] PLINK prune.in file not found" >&2
        exit 1
    fi
    
    n_indep=\$(wc -l < effective_number.prune.in)
    if [ \$n_indep -eq 0 ]; then
        echo "[ERROR] No independent markers found" >&2
        exit 1
    fi
    
    awk -v n=\$n_indep 'BEGIN{printf "%.12f\\n", 1/n}' > threshold.txt
    echo "[INFO] Independent markers: \$n_indep, Threshold: \$(cat threshold.txt)" >> ${vcf.simpleName}.log
    """
}


process SPLIT_VCF {
    publishDir "${params.outdir}/split_vcf/" 
    input:
    path vcf
    output:
    path("${vcf.baseName}.biallelic.vcf")
    script:
    """
    bcftools norm -m -both ${vcf} -o ${vcf.baseName}.biallelic.vcf
    """
}

process SNP_VCF_TO_PLINK {
    tag "SNP_kinship_prep"
    label 'process_medium'
    publishDir "${params.outdir}/kinship_plink", mode: 'copy'
    
    input:
    path snp_vcf
    
    output:
    tuple path("${snp_vcf.baseName}.bed"), 
          path("${snp_vcf.baseName}.bim"), 
          path("${snp_vcf.baseName}.fam"),
          val("${snp_vcf.baseName}")
    
    script:
    """
    ${params.plink} \\
        --vcf ${snp_vcf} \\
        --double-id \\
        --make-bed \\
        --out ${snp_vcf.baseName} \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --maf ${params.maf}
    """
}

process PLINK_VCF_CONVERSION {
    tag "${phenotype.baseName}"
    label 'process_medium'
    
    input:
    tuple path(vcf), path(phenotype)
    
    output:
    tuple val("${phenotype.baseName}"), 
          path("${phenotype.baseName}.bed"), 
          path("${phenotype.baseName}.bim"), 
          path("${phenotype.baseName}.fam"),
          path(phenotype)
    
    script:
    """
    ${params.plink} \\
        --vcf ${vcf} \\
        --double-id \\
        --pheno ${phenotype} \\
        --mpheno 1 \\
        --make-bed \\
        --out ${phenotype.baseName} \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --maf ${params.maf}
    """
}

process GEMMA_KINSHIP {
    tag "kinship_calculation"
    label 'process_medium'
    publishDir "${params.outdir}/kinship", 
               mode: 'copy', 
               pattern: "*.cXX.txt"
    
    input:
    tuple path(bed), path(bim), path(fam), val(bfile_prefix)
    
    output:
    path("kinship.cXX.txt")
    
    script:
    """
    # Create a dummy phenotype file for kinship calculation
    # GEMMA requires a phenotype file even for -gk mode
    # Format: FID IID phenotype (we use 1 as dummy value for all individuals)
    awk '{print \$1, \$2, 1}' ${fam} > dummy.pheno
    
    ${params.gemma} \\
        -bfile ${bfile_prefix} \\
        -p dummy.pheno \\
        -gk 1 \\
        -outdir . \\
        -o kinship \\
        -miss ${params.miss_threshold} \\
        -maf ${params.maf} \\
        -hwe ${params.hwe_threshold}
    """
}

process GEMMA_LMM_ASSOCIATION {
    tag "${phenotype_name}"
    label 'process_medium'
    publishDir "${params.outdir}/results/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.assoc.txt"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(phenotype), path(kinship)
    
    output:
    tuple val(phenotype_name), 
          path(bed), 
          path(bim), 
          path(fam),
          path("gemma_lmm.assoc.txt"),
          path(phenotype)
    
    script:
    """
    ${params.gemma} \\
        -bfile ${phenotype_name} \\
        -k ${kinship} \\
        -outdir . \\
        -o gemma_lmm \\
        -lmm ${params.gemma_lmm_type} \\
        -miss ${params.miss_threshold} \\
        -maf ${params.maf}
    """
}

process MANHATTAN_PLOT {
    tag "${phenotype_name}"
    label 'process_low'
    publishDir "${params.outdir}/plots/", 
               mode: 'copy', 
               pattern: "*.manhattan.png"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype), val(threshold)
    
    output:
    path("${phenotype_name}.manhattan.png")
    
    script:
    """
    python3 ${projectDir}/scripts/manhattan_plot.py \\
        ${assoc} \\
        ${phenotype_name}.manhattan.png \\
        ${phenotype_name} \\
        ${threshold}
    """
}

process PLINK_CLUMPING {
    tag "${phenotype_name}_${params.clump_p2}"
    label 'process_low'
    publishDir "${params.outdir}/lead_markers", 
               mode: 'copy', 
               pattern: "*.clumped"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype), val(threshold)
    
    output:
    tuple val(phenotype_name), path("${phenotype_name}.clumped")
    
    script:
    """
    set +e  # Disable exit on error for clumping
    ${params.plink} \\
        --bfile ${phenotype_name} \\
        --clump ${assoc} \\
        --clump-p1 ${threshold} \\
        --clump-p2 ${params.clump_p2} \\
        --clump-r2 ${params.clump_r2} \\
        --clump-kb ${params.clump_kb} \\
        --allow-no-sex \\
        --allow-extra-chr \\
        --out ${phenotype_name} \\
        --clump-snp-field rs \\
        --clump-field p_wald
    
    # If clumping failed or produced no output, create empty clumped file
    if [ ! -f ${phenotype_name}.clumped ]; then
        touch ${phenotype_name}.clumped
    fi
    set -e  # Re-enable exit on error
    """
}

process EXTRACT_LEAD_MARKERS {
    tag "extract_lead_markers"
    label 'process_low'
    publishDir "${params.outdir}/lead_markers", 
               mode: 'copy', 
               pattern: "all_lead_markers.tsv"
    
    input:
    path clumped_files
    
    output:
    path("all_lead_markers.tsv"), emit: markers_tsv
    
    script:
    """
    # Write header
    echo -e "Marker\\tTraits" > all_lead_markers.tsv
    
    # Process each clumped file
    # Files are copied to work directory, we can process them directly
    for clumped in *.clumped; do
        if [ -f "\$clumped" ] && [ -s "\$clumped" ]; then
            # Extract phenotype name from filename (remove .clumped extension)
            trait_name=\$(basename "\$clumped" .clumped)
            
            # Extract third column (SNP/marker ID), skip header
            # Use awk to handle whitespace-separated columns
            awk 'NR > 1 && \$3 != "" {print \$3 "\\t" trait}' trait="\$trait_name" "\$clumped" >> all_lead_markers.tsv
        fi
    done
    """
}

process GET_UNIQUE_MARKERS {
    tag "get_unique_markers"
    label 'process_low'
    publishDir "${params.outdir}/lead_markers", 
               mode: 'copy', 
               pattern: "unique_markers.list"
    
    input:
    path markers_tsv
    
    output:
    path("unique_markers.list"), emit: markers_list
    
    script:
    """
    # Extract first column (Marker), skip header, sort and uniq
    tail -n +2 ${markers_tsv} | cut -f1 | sort -V | uniq > unique_markers.list
    """
}

process FILTER_VCF_BY_MARKERS {
    tag "filter_vcf_by_markers"
    label 'process_medium'
    publishDir "${params.outdir}/lead_markers", 
               mode: 'copy', 
               pattern: "*.lead_markers.vcf.gz*"
    
    input:
    path vcf
    path markers_list
    
    output:
    path("${vcf.baseName}.lead_markers.vcf.gz"), emit: filtered_vcf
    path("${vcf.baseName}.lead_markers.vcf.gz.csi"), emit: filtered_vcf_index
    
    script:
    """
    # Filter VCF by marker IDs using bcftools
    # -i 'ID=@file' filters to keep only variants whose ID is in the file
    bcftools view -i 'ID=@${markers_list}' ${vcf} -Oz -o ${vcf.baseName}.lead_markers.vcf.gz
    
    # Index the filtered VCF
    bcftools index ${vcf.baseName}.lead_markers.vcf.gz
    """
}