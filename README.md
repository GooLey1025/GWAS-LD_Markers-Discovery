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

# Reproduce of Rice GraphPangenome GS

This section documents the complete workflow used for the Rice GraphPangenome GS paper submission. The pipeline generates comprehensive marker sets across multiple variant types (SNP, INDEL, and SV) using two complementary strategies:

1. **GWAS LD-based lead marker discovery**: Identifies trait-associated markers through GWAS analysis across multiple cohorts and phenotypes, followed by LD clumping to extract lead markers. This approach prioritizes markers with significant associations to phenotypic traits.

2. **Genome-wide LD-based lead marker discovery**: Performs genome-wide LD pruning to identify independent markers across the entire genome, regardless of trait associations. This approach ensures comprehensive genome coverage and marker representativeness.

The workflow processes multiple rice cohorts with varying sample sizes and variant compositions:
- **SNP-only cohorts** (404rice, 1439rice, 176rice, 378rice, 3023rice, 532rice): Analyzed using SNP variants only
- **Multi-variant cohorts** (1171rice, 705rice): Analyzed using combined SNP, INDEL, and SV variants

All commands below reproduce the exact analysis pipeline used in the paper. 

#### GWAS LD-based lead marker discovery across cohorts and traits
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

##### SNP && INDEL && SV
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

#### Genome-wide LD-based lead markers discovery

##### SNP-only 
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

##### SNP && INDEL && SV
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
