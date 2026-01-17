# GWAS Lead Marker Discovery Pipeline

This repository provides a Nextflow (DSL2) workflow for GWAS-based lead marker discovery using PLINK and GEMMA. The pipeline performs per-phenotype linear mixed model (LMM) association analyses while sharing a single kinship matrix across all traits, followed by multiple-testing threshold estimation, Manhattan plotting, LD clumping, and final extraction of lead markers to generate a reduced VCF.

## Overview

The workflow implements a comprehensive GWAS pipeline with the following key features:

### Key Workflow Steps

1. **Independent Marker Estimation**: Calculate the number of independent markers using PLINK LD pruning to determine multiple-testing correction threshold
2. **Kinship Matrix Calculation**: Compute a single kinship matrix from SNP-only VCF data using GEMMA (shared across all phenotypes)
3. **Per-Phenotype Association Analysis**: Perform LMM association testing for each phenotype using GEMMA
4. **Manhattan Plot Generation**: Visualize association results with appropriate significance thresholds
5. **LD Clumping**: Identify lead markers using PLINK clumping based on association p-values
6. **Lead Marker Extraction**: Extract unique lead markers across all traits and generate filtered VCF

### Variant Type Support

- **SNP-only cohorts**: Analysis using SNP variants only
- **Multi-variant cohorts**: Joint analysis of SNP, INDEL, and SV variants
- The VCF must be biallelic and have properly populated ID columns indicating variant type (SNP-, INDEL-, SV-)

## Input Data Requirements

### VCF Inputs

#### Association Testing VCF (`--vcf`)

- **Purpose**: Used for per-phenotype GWAS analysis in GEMMA
- **Requirements**:
  - Must be biallelic
  - ID column must indicate variant type (SNP-, INDEL-, SV-)
  - Can be SNP-only or merged VCF containing SNPs, INDELs, and SVs together (example: `test.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz`)
- **Usage**: All variant types present in this VCF will be jointly tested in association analysis

#### Kinship Calculation VCF (`--snp_vcf`)

- **Purpose**: Used exclusively for kinship matrix calculation in GEMMA
- **Requirements**: SNP-only VCF
- **Rationale**: Ensures kinship matrix is estimated from high-confidence SNP markers only, while association testing can be performed on a broader variant set

### Phenotype Input Format

Phenotype files must be provided in the directory specified by `--phenotypes_dir`:

- **Format**: Tab-separated values (TSV) with header
- **File naming**: `*.tsv` extension
- **Required columns** (exactly 3):
  1. `FID` – Family ID
  2. `IID` – Individual ID
  3. Phenotype value (quantitative trait)

**Example phenotype file**:
```txt
FID     IID     Gel_consistency
ZJM81   ZJM81   82.5
ZJM84   ZJM84   85.0
ZJM85   ZJM85   44.0
ZJM87   ZJM87   80.0
ZJM88   ZJM88   48.0
```

### Trait Type Limitation

- **Currently supported**: Quantitative traits only
- **Not supported**: Categorical or binary traits
- **Note**: Using classification traits may lead to errors, particularly during Manhattan plot generation

## Requirements

The following software must be installed and available in the **system environment ($PATH)** before running the pipeline:

- `nextflow`
- `plink`
- `gemma` (GEMMA executable)
- `bcftools`
- `python3` (for Manhattan plot generation)

**Note**: The executable names are case-sensitive and must exactly match those listed above.

## Configuration

All parameters can be provided either via the command line or a Nextflow configuration file. 

### Input Parameters

| Parameter | Description | Required/Default |
|-----------|-------------|------------------|
| `--vcf` | Association testing VCF file path (biallelic, with ID column) | Required |
| `--snp_vcf` | SNP-only VCF file path for kinship calculation | Required |
| `--phenotypes_dir` | Directory containing phenotype TSV files | Required |
| `--outdir` | Output directory for all results | Required |
| `--split_vcf` | Set to true if VCF needs to be split (default: false, assuming VCF is already biallelic) | `false` |
| `--maf` | Minimum allele frequency threshold | `0.05` |
| `--miss_threshold` | Maximum missing data rate threshold | `0.5` |
| `--gemma_lmm_type` | GEMMA LMM type | `1` |
| `--hwe_threshold` | Hardy-Weinberg equilibrium threshold | `0` |
| `--plink_indep_r2` | R² threshold for PLINK independent marker estimation | `0.2` |
| `--clump_p2` | Secondary p-value threshold for clumping | `0.05` |
| `--clump_r2` | R² threshold for clumping | `0.1` |
| `--clump_kb` | Physical distance threshold for clumping (kb) | `1000` |

**Note**: The primary clumping p-value threshold (`clump_p1`) is automatically calculated as 1/number_of_independent_markers.

## Command Examples

### Basic Usage

Run the complete pipeline with default parameters:

```bash
nextflow run main.nf \
    --vcf 705rice.graph.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz \
    --snp_vcf 705rice.graph.0.5_0.05.full.snp.impute.biallelic.id.vcf.gz \
    --phenotypes_dir phenotypes/705rice \
    --outdir 705rice_graph.0.5_0.05.SNP_INDEL_SV
```

### Resume Execution

If the pipeline fails or is interrupted, resume from the last successful step:

```bash
nextflow run main.nf -resume [previous parameters...]
```

## Output Structure

All outputs are published into structured subfolders under the configured output directory:

- `plink_indep/` - PLINK independent marker estimation results and threshold calculation
- `kinship_plink/` - PLINK format files for kinship calculation
- `kinship/` - Computed kinship matrix (shared across all phenotypes)
- `results/<phenotype_name>/` - Per-phenotype GEMMA association results
- `plots/` - Manhattan plots for each phenotype
- `lead_markers/` - Clumped lead markers and filtered VCF files
  - `*.clumped` - PLINK clumping results per phenotype
  - `all_lead_markers.tsv` - Combined lead markers across all traits
  - `unique_markers.list` - Unique marker IDs
  - `*.lead_markers.vcf.gz` - Filtered VCF containing only lead markers

<hr style="border: 1px dashed #ccc;">

# Reproduce of Rice GraphPangenome GS

This section documents the complete workflow used for the Rice GraphPangenome GS paper submission. The pipeline generates comprehensive marker sets across multiple variant types (SNP, INDEL, and SV) using two complementary strategies:

1. **GWAS LD-based lead marker discovery**: Identifies trait-associated markers through GWAS analysis across multiple cohorts and phenotypes, followed by LD clumping to extract lead markers. This approach prioritizes markers with significant associations to phenotypic traits.

2. **Genome-wide LD-based lead marker discovery**: Performs genome-wide LD pruning to identify independent markers across the entire genome, regardless of trait associations. This approach ensures comprehensive genome coverage and marker representativeness.

3. **Public GWAS data Re-Calling and LD filter**: Re-genotypes previously published GWAS datasets (e.g., GWAS Atlas) using GATK3.7 UnifiedGenotyper in discovery mode (`-gt_mode DISCOVERY`) with predefined marker intervals. The re-called variants are then filtered using LD pruning to identify independent markers. This approach leverages existing GWAS findings and adapts them to the target population's LD structure.

4. **RiceNavi library QTN sites**: Incorporates curated quantitative trait nucleotide (QTN) sites from the RiceNavi database, which contains well-validated SNP and INDEL variants associated with important agronomic traits. These pre-processed sites provide high-confidence markers with established biological significance.

The workflow processes multiple rice cohorts with varying sample sizes and variant compositions:
- **SNP-only cohorts** (404rice, 1439rice, 176rice, 378rice, 3023rice, 532rice): Analyzed using SNP variants only
- **Multi-variant cohorts** (1171rice, 705rice): Analyzed using combined SNP, INDEL, and SV variants

All commands below reproduce the exact analysis pipeline used in the paper. 

## VCF prepared

### GWAS LD-based lead marker discovery across cohorts and traits
```sh
nextflow run main.nf --vcf input_vcfs/ws_A_Population.biallelic.id.vcf --snp_vcf input_vcfs/ws_A_Population.biallelic.id.vcf --outdir 404rice --phenotypes_dir phenotypes/404rice 
nextflow run main.nf --vcf input_vcfs/ws_Huang_NC2015_POP.biallelic.id.vcf --snp_vcf input_vcfs/ws_Huang_NC2015_POP.biallelic.id.vcf --outdir 1439rice --phenotypes_dir phenotypes/1439rice
nextflow run main.nf --vcf input_vcfs/ws_Masuta.biallelic.id.vcf --snp_vcf input_vcfs/ws_Masuta.biallelic.id.vcf --outdir 176rice --phenotypes_dir phenotypes/176rice
nextflow run main.nf --vcf input_vcfs/ws_Susan_NC2011_378rice.biallelic.id.vcf --snp_vcf input_vcfs/ws_Susan_NC2011_378rice.biallelic.id.vcf --outdir 378rice --phenotypes_dir phenotypes/378rice
nextflow run main.nf --vcf input_vcfs/ws_3k_Rice_3023rice.biallelic.id.vcf --snp_vcf input_vcfs/ws_3k_Rice_3023rice.biallelic.id.vcf --outdir 3023rice --phenotypes_dir phenotypes/3023rice

nextflow run main.nf --vcf input_vcfs/532rice.0.5_0.05.snp.impute.biallelic.id.vcf --snp_vcf input_vcfs/532rice.0.5_0.05.snp.impute.biallelic.id.vcf --phenotypes_dir phenotypes/532rice --outdir 532rice

nextflow run main.nf --vcf input_vcfs/1171rice.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz --snp_vcf input_vcfs/1171rice.0.5_0.05.full.snp.impute.biallelic.id.vcf.gz --outdir 1171rice --phenotypes_dir phenotypes/1171rice
nextflow run main.nf --vcf 705rice.graph.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz --snp_vcf 705rice.graph.0.5_0.05.full.snp.impute.biallelic.id.vcf.gz --outdir 705rice --phenotypes_dir phenotypes/705rice

cohorts="404rice 1439rice 176rice 378rice 3023rice 532rice 1171rice 705rice"
export markers_dir=markers_d

mkdir -p $markers_dir
mkdir -p tmp
parallel -j 12 '
awk "FNR>1 && \$3 != \"\" {print \$3}" {}/lead_markers/*.clumped \
  > '"$markers_dir"'/{}.GWAS_LD.markers.txt
' ::: $cohorts
```

#### SNP-only
```sh
mkdir -p tmp

snp_cohorts="404rice 1439rice 176rice 378rice 3023rice 532rice"
parallel -j 12 --halt now,fail=1 '
  cohort="{}"
  files=(${cohort}/lead_markers/*.vcf.gz)

  if (( ${#files[@]} == 0 )); then
    echo "[WARN] No input VCF found for ${cohort} under ${cohort}/lead_markers/*.vcf.gz" >&2
    exit 0
  fi

  if (( ${#files[@]} > 1 )); then
    echo "[WARN] More than one VCF found for ${cohort}; using the first one: ${files[0]}" >&2
  fi
  in_vcf="${files[0]}"

  # Decompress to a plain VCF (not gz)
  bcftools view -o "tmp/${cohort}.GWAS_LD.vcf.gz" "${in_vcf}"

  bcftools view \
    -i "ID=@${markers_dir}/${cohort}.GWAS_LD.markers.txt" \
    -o "${markers_dir}/${cohort}.snp.GWAS_LD.vcf" \
    "tmp/${cohort}.GWAS_LD.vcf.gz"
' ::: $snp_cohorts

```

#### SNP && INDEL && SV
```sh
all_cohorts="1171rice 705rice"

parallel -j 12 --halt now,fail=1 '
  cohort="{}"
  files=(${cohort}/lead_markers/*.vcf.gz)

  if (( ${#files[@]} == 0 )); then
    echo "[WARN] No input VCF found for ${cohort} under ${cohort}/lead_markers/*.vcf.gz" >&2
    exit 0
  fi

  if (( ${#files[@]} > 1 )); then
    echo "[WARN] More than one VCF found for ${cohort}; using the first one: ${files[0]}" >&2
  fi
  in_vcf="${files[0]}"

  bcftools view -Ov -o "tmp/${cohort}.GWAS_LD.vcf.gz" "${in_vcf}"

  # SNP
  bcftools view \
    -i "ID=@${markers_dir}/${cohort}.GWAS_LD.markers.txt && ID~\"^SNP-\"" \
    -Ov -o "${markers_dir}/${cohort}.snp.GWAS_LD.vcf" \
    "tmp/${cohort}.GWAS_LD.vcf.gz"

  # INDEL
  bcftools view \
    -i "ID=@${markers_dir}/${cohort}.GWAS_LD.markers.txt && ID~\"^INDEL-\"" \
    -Ov -o "${markers_dir}/${cohort}.indel.GWAS_LD.vcf" \
    "tmp/${cohort}.GWAS_LD.vcf.gz"

  # SV
  bcftools view \
    -i "ID=@${markers_dir}/${cohort}.GWAS_LD.markers.txt && ID~\"^SV-\"" \
    -Ov -o "${markers_dir}/${cohort}.sv.GWAS_LD.vcf" \
    "tmp/${cohort}.GWAS_LD.vcf.gz"
' ::: $all_cohorts
```

### Genome-wide LD-based lead markers discovery

#### SNP-only 
```sh
export markers_dir=markers_d
snp_cohorts="404rice 1439rice 176rice 378rice 3023rice 532rice"

mkdir -p plink_tmp_d
plink --threads 32 --vcf input_vcfs/ws_A_Population.biallelic.id.vcf --indep-pairwise 1000 50 0.02 --out plink_tmp_d/404rice.snp.WG_ld --double-id
plink --vcf input_vcfs/ws_Huang_NC2015_POP.biallelic.id.vcf --indep-pairwise 1000 50 0.02 --out plink_tmp_d/1439rice.snp.WG_ld --double-id
plink --vcf input_vcfs/ws_Matsuoka.biallelic.id.vcf --indep-pairwise 1000 50 0.02 --out plink_tmp_d/176rice.snp.WG_ld --double-id
plink --vcf input_vcfs/ws_Susan_NC2011_378rice.biallelic.id.vcf --indep-pairwise 1000 50 0.02 --out plink_tmp_d/378rice.snp.WG_ld --double-id
plink --vcf input_vcfs/ws_3k_Rice_3023rice.biallelic.id.vcf --indep-pairwise 1000 50 0.02 --out plink_tmp_d/3023rice.snp.WG_ld --double-id
plink --vcf input_vcfs/532rice.0.5_0.05.snp.impute.biallelic.id.vcf --indep-pairwise 1000 50 0.02 --out plink_tmp_d/532rice.snp.WG_ld --double-id

parallel -j 12 '
cat plink_tmp_d/{}.snp.WG_ld.prune.in > $markers_dir/{}.snp.WG_ld.markers.txt
' ::: $snp_cohorts


bcftools view -i "ID=@${markers_dir}/404rice.snp.WG_ld.markers.txt" input_vcfs/ws_A_Population.biallelic.id.vcf -o $markers_dir/404rice.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/1439rice.snp.WG_ld.markers.txt" input_vcfs/ws_Huang_NC2015_POP.biallelic.id.vcf -o $markers_dir/1439.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/176rice.snp.WG_ld.markers.txt" input_vcfs/ws_Matsuoka.biallelic.id.vcf -o $markers_dir/176rice.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/378rice.snp.WG_ld.markers.txt" input_vcfs/ws_Susan_NC2011_378rice.biallelic.id.vcf -o $markers_dir/378rice.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/3023rice.snp.WG_ld.markers.txt" input_vcfs/ws_3k_Rice_3023rice.biallelic.id.vcf -o $markers_dir/3023rice.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/532rice.snp.WG_ld.markers.txt" input_vcfs/532rice.0.5_0.05.snp.impute.biallelic.id.vcf -o $markers_dir/532rice.snp.WG.LeadMarkers.vcf
```

#### SNP && INDEL && SV
```sh
all_cohorts="1171rice 705rice"

parallel -j 12 '
bcftools view --threads 32 -i "ID~\"^SNP-\"" -Oz -o plink_tmp_d/{}.0.5_0.05.full.snp.impute.biallelic.id.format.vcf.gz input_vcfs/{}.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz
bcftools view --threads 32 -i "ID~\"^INDEL-\"" -Oz -o plink_tmp_d/{}.0.5_0.05.full.indel.impute.biallelic.id.format.vcf.gz input_vcfs/{}.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz
bcftools view --threads 32 -i "ID~\"^SV-\"" -Oz -o plink_tmp_d/{}.0.5_0.05.full.sv.impute.biallelic.id.format.vcf.gz input_vcfs/{}.0.5_0.05.full.all.impute.biallelic.id.format.vcf.gz

plink --vcf plink_tmp_d/{}.0.5_0.05.full.snp.impute.biallelic.id.format.vcf.gz --indep-pairwise 1000 50 0.02 --out plink_tmp_d/{}.snp.WG_ld --double-id
plink --vcf plink_tmp_d/{}.0.5_0.05.full.indel.impute.biallelic.id.format.vcf.gz --indep-pairwise 1000 50 0.02 --out plink_tmp_d/{}.indel.WG_ld --double-id
plink --vcf plink_tmp_d/{}.0.5_0.05.full.sv.impute.biallelic.id.format.vcf.gz --indep-pairwise 1000 50 0.02 --out plink_tmp_d/{}.sv.WG_ld --double-id

cat plink_tmp_d/{}.snp.WG_ld.prune.in > $markers_dir/{}.snp.WG_ld.markers.txt
cat plink_tmp_d/{}.indel.WG_ld.prune.in > $markers_dir/{}.indel.WG_ld.markers.txt
cat plink_tmp_d/{}.sv.WG_ld.prune.in > $markers_dir/{}.sv.WG_ld.markers.txt
' ::: $all_cohorts

bcftools view -i "ID=@${markers_dir}/1171rice.snp.WG_ld.markers.txt" plink_tmp_d/1171rice.0.5_0.05.full.snp.impute.biallelic.id.format.vcf.gz -o $markers_dir/1171rice.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/1171rice.indel.WG_ld.markers.txt" plink_tmp_d/1171rice.0.5_0.05.full.indel.impute.biallelic.id.format.vcf.gz -o $markers_dir/1171rice.indel.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/1171rice.sv.WG_ld.markers.txt" plink_tmp_d/1171rice.0.5_0.05.full.sv.impute.biallelic.id.format.vcf.gz -o $markers_dir/1171rice.sv.WG.LeadMarkers.vcf

bcftools view -i "ID=@${markers_dir}/705rice.snp.WG_ld.markers.txt" plink_tmp_d/705rice.0.5_0.05.full.snp.impute.biallelic.id.format.vcf.gz -o $markers_dir/705rice.snp.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/705rice.indel.WG_ld.markers.txt" plink_tmp_d/705rice.0.5_0.05.full.indel.impute.biallelic.id.format.vcf.gz -o $markers_dir/705rice.indel.WG.LeadMarkers.vcf
bcftools view -i "ID=@${markers_dir}/705rice.sv.WG_ld.markers.txt" plink_tmp_d/705rice.0.5_0.05.full.sv.impute.biallelic.id.format.vcf.gz -o $markers_dir/705rice.sv.WG.LeadMarkers.vcf
```

### Public data (GWAS Atlas) Re-Calling and LD-filter

For multi-variant cohorts (1171rice and 705rice), we re-genotyped previously published GWAS datasets using GATK3.7 UnifiedGenotyper. The re-called variants were then filtered through LD pruning to identify independent markers. These files are linked to standardized names for downstream analysis:

```sh
ln -s /data3/home/gulei/projects/GraphPan/UnifiedGenotyper_NF/705rice_vcf_out/705rice.public_data.snp.impute.biallelic.rename.id.ld_filter.vcf 705rice.snp.public_data_GWAS_LD.plink_prune.vcf

ln -s /data3/home/gulei/projects/GraphPan/UnifiedGenotyper_NF/1171rice_vcf_out/1171rice.public_data.snp.impute.biallelic.rename.id.ld_filter.vcf 1171rice.snp.public_data_GWAS_LD.plink_prune.vcf
```

### RiceNavi data sites VCF

The RiceNavi database provides curated quantitative trait nucleotide (QTN) sites with established associations to agronomic traits. Pre-processed site VCFs are used directly:

- `RiceNavi.snp.sites.vcf`: Curated SNP sites from RiceNavi
- `RiceNavi.indel.sites.vcf`: Curated INDEL sites from RiceNavi

These files have been pre-processed and are ready for use in the marker panel generation pipeline.

## Sites VCF and Intervals Prepared for Genotyping

Although we have generated various VCF files through different strategies, not all of them are suitable for generating a comprehensive population variant panel. The final marker selection integrates multiple complementary sources to ensure both trait relevance and genome-wide coverage.

### Marker Selection Strategy

For the 1171rice population (and similarly for 705rice), the final marker panel combines markers from four complementary sources:

1. **Multi-cohort, multi-trait GWAS lead variants (LD) markers**: Lead markers identified through GWAS analysis across multiple cohorts and traits, capturing trait-associated variants with strong statistical support.

2. **Population-specific genome-wide LD markers**: Independent markers identified through genome-wide LD pruning within the target population (1171rice), ensuring comprehensive genome coverage and population-specific LD structure representation.

3. **Re-genotyped public GWAS data LD markers**: Variants re-genotyped from previously published public datasets (e.g., GWAS Atlas) using the target population's samples, then filtered through LD pruning. This approach leverages existing GWAS findings while adapting them to the target population's genetic structure.

4. **RiceNavi curated QTN sites**: Well-validated quantitative trait nucleotide sites from the RiceNavi database, providing high-confidence markers with established biological significance for important agronomic traits.

### Marker Panel Generation

The `summary.py` script integrates all marker sources and generates the final marker panel. The following example demonstrates the command structure and input VCF files used for the 1171rice population:

```sh
P=1171rice # 705rice

python3 summary.py --snp-vcf $P.snp.public_data_GWAS_LD.plink_prune.vcf \
    RiceNavi.snp.sites.vcf \
    1171rice.snp.WG.LeadMarkers.vcf \
    705rice.snp.GWAS_LD.vcf \
    532rice.snp.GWAS_LD.vcf \
    404rice.snp.GWAS_LD.vcf \
    378rice.snp.GWAS_LD.vcf \
    3023rice.snp.GWAS_LD.vcf \
    176rice.snp.GWAS_LD.vcf \
    1439rice.snp.GWAS_LD.vcf \
    1171rice.snp.GWAS_LD.vcf \
    --indel-vcf $P.indel.WG.LeadMarkers.vcf \
    705rice.indel.GWAS_LD.vcf \
    1171rice.indel.GWAS_LD.vcf \
    RiceNavi.indel.sites.vcf \
    --sv-vcf 1171rice.sv.GWAS_LD.info.vcf \
    705rice.sv.GWAS_LD.info.vcf \
    $P.sv.WG.LeadMarkers.info.vcf \
    -o $P -p $P
```

### Output Files and Downstream Usage

The output files from `summary.py` are used as input for the [GATK-DELLY-Allele_based-Genotyping](https://github.com/GooLey1025/GATK-DELLY-Allele_based-Genotyping) pipeline, which performs allele-based genotyping across new samples.

**Important Notes:**

1. **Input SV VCF Requirements**: For structural variant (SV) VCFs, the `INFO` column must be retained from the original Delly-generated VCF. The INFO column contains critical structural variant annotations required for proper genotyping. If the INFO column is missing, variants will not be considered during genotyping, although the program will continue to run.

2. **INFO Column Restoration**: The INFO column is typically discarded during Beagle imputation. To restore it, use `bcftools annotate` to copy the INFO field from the original Delly VCF. Example:

```sh
bcftools annotate \
    -a /path/to/GATK-DELLY-VariantsCalling/1171_sv_run/1171rice/delly_no_germline/0.5_0.05.filtered.merged.genotype.id.samples_sort.multiallelic.vcf.gz \
    -c INFO \
    -o 1171rice.sv.GWAS_LD.info.vcf \
    1171rice.sv.GWAS_LD.vcf.gz
```

This command annotates the filtered SV VCF (`1171rice.sv.GWAS_LD.vcf.gz`) with INFO fields from the original Delly VCF, ensuring all necessary structural variant metadata is preserved for downstream genotyping.
