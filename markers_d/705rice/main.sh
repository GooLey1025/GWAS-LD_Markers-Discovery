bcftools annotate -a /data3/home/gulei/projects/GraphPan/variants_calling-gatk_delly/1171_sv_run/1171rice/delly_no_germline/0.5_0.05.filtered.merged.genotype.id.samples_sort.multiallelic.vcf.gz -c INFO -o 1171rice.sv.GWAS_LD.info.vcf 1171rice.sv.GWAS_LD.vcf.gz 

python3 summary.py --snp-vcf 705rice.snp.public_data_GWAS_LD.plink_prune.vcf 705rice.snp.WG.LeadMarkers.vcf \
    705rice.snp.GWAS_LD.vcf 532rice.snp.GWAS_LD.vcf 404rice.snp.GWAS_LD.vcf 378rice.snp.GWAS_LD.vcf 3023rice.snp.GWAS_LD.vcf \
    176rice.snp.GWAS_LD.vcf 1439rice.snp.GWAS_LD.vcf 1171rice.snp.GWAS_LD.vcf RiceNavi.snp.sites.vcf \
    --indel-vcf 705rice.indel.WG.LeadMarkers.vcf 705rice.indel.GWAS_LD.vcf \
    1171rice.indel.GWAS_LD.vcf RiceNavi.indel.sites.vcf \
    --sv-vcf 1171rice.sv.GWAS_LD.info.vcf 705rice.sv.GWAS_LD.info.vcf 705rice.sv.WG.LeadMarkers.info.vcf \
    -o 705rice -p 705rice

