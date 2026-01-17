#!/usr/bin/env python3
"""
Summarize VCF files: extract chr and pos, and add metadata columns
based on filename patterns.

This script processes VCF files in a directory and generates:
- A summary TSV file with variant positions and metadata
- Interval files for SNP and INDEL markers
- Sites-only VCF files for SNP, INDEL, and SV variants
"""

import argparse
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Set, Optional

# Default blacklist of positions to exclude from intervals and VCF files
# Format: (chr_num, pos) tuples
DEFAULT_POSITION_BLACKLIST = {
    ('2', '35987830'),
    ('3', '36436858'),
    ('11', '29330602'),
    ('11', '30732291')
}

def parse_filename(filename):
    """
    Parse filename to extract source, type, and population.
    
    Returns: (source, type, population)
    """
    basename = os.path.basename(filename)
    
    # Determine source
    if 'RiceNavi' in basename:
        source = 'RiceNavi'
    elif 'public_data_GWAS_LD' in basename:
        source = 'Public_GWAS_LD'
    elif 'GWAS_LD' in basename:
        source = 'GWAS_LD'
    elif 'WG.LeadMarkers' in basename:
        source = 'WG_LD'
    else:
        source = 'Unknown'
    
    # Determine type
    if '.snp.' in basename:
        var_type = 'snp'
    elif '.indel.' in basename:
        var_type = 'indel'
    elif '.sv.' in basename:
        var_type = 'sv'
    else:
        var_type = 'unknown'
    
    # Determine population
    if source == 'RiceNavi' or source == 'Public_GWAS_LD':
        population = 'NA'
    else:
        # Extract population from filename (e.g., 705rice, 532rice, etc.)
        match = re.search(r'(\d+rice)', basename)
        if match:
            population = match.group(1)
        else:
            population = 'NA'
    
    return source, var_type, population

def extract_vcf_positions(vcf_file, is_rice_navi=False):
    """
    Extract chr, pos, and optionally type from VCF file.
    
    For RiceNavi files, also extract type from ID field.
    
    Returns: list of (chr, pos, type) tuples, where type is None for non-RiceNavi files
    """
    positions = []
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip header lines
                if line.startswith('##'):
                    continue
                # Skip empty lines
                if not line:
                    continue
                # Skip column header line (but it starts with #CHROM)
                if line.startswith('#CHROM'):
                    continue
                
                # Parse data line
                fields = line.split('\t')
                if len(fields) >= 2:
                    chr_col = fields[0]
                    pos_col = fields[1]
                    # Validate that pos is numeric
                    try:
                        pos = int(pos_col)
                        
                        # For RiceNavi, extract type from ID field
                        var_type = None
                        if is_rice_navi and len(fields) >= 3:
                            id_field = fields[2]
                            if id_field.startswith('SNP-'):
                                var_type = 'snp'
                            elif id_field.startswith('INDEL-'):
                                var_type = 'indel'
                            elif id_field.startswith('SV-'):
                                var_type = 'sv'
                            # If type not found in ID, try to infer from REF/ALT
                            if var_type is None and len(fields) >= 5:
                                ref = fields[3]
                                alt = fields[4]
                                # Simple heuristic: if lengths differ significantly, it's an indel
                                if abs(len(ref) - len(alt)) > 1:
                                    var_type = 'indel'
                                elif len(ref) == 1 and len(alt) == 1:
                                    var_type = 'snp'
                                else:
                                    var_type = 'unknown'
                        
                        positions.append((chr_col, pos, var_type))
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading {vcf_file}: {e}", file=sys.stderr)
    
    return positions

def extract_chr_number(chr_name):
    """
    Extract chromosome number from chromosome name, removing any prefix.
    Returns the chromosome number as a string (e.g., "1", "2", etc.)
    """
    # Remove common prefixes (case-insensitive)
    chr_name = re.sub(r'^[Cc]hr[oO]m[oO]s[oO]me[_-]?', '', chr_name, flags=re.IGNORECASE)
    chr_name = re.sub(r'^[Cc]hr[oO]m[oO]s[oO]me', '', chr_name, flags=re.IGNORECASE)
    chr_name = re.sub(r'^[Cc]hr[oO]m', '', chr_name, flags=re.IGNORECASE)
    chr_name = re.sub(r'^[Cc]hr', '', chr_name, flags=re.IGNORECASE)
    chr_name = re.sub(r'^chr', '', chr_name, flags=re.IGNORECASE)
    chr_name = re.sub(r'^Chr', '', chr_name)
    # Extract number if present
    match = re.search(r'(\d+)', chr_name)
    if match:
        return match.group(1)
    # If no number found, return as is (for non-standard chromosomes)
    return chr_name

def make_biallelic(ref, alt):
    """
    Convert multi-allelic REF/ALT to biallelic by keeping only the first ALT allele.
    
    Args:
        ref: Reference allele (should be single allele)
        alt: Alternate allele(s), may contain multiple alleles separated by commas
    
    Returns:
        tuple: (biallelic_ref, biallelic_alt)
    """
    # If REF contains multiple alleles (shouldn't happen normally, but handle it)
    if ',' in ref:
        ref = ref.split(',')[0]
    
    # If ALT contains multiple alleles, keep only the first one
    if ',' in alt:
        alt = alt.split(',')[0]
    
    return ref, alt

def load_blacklist(blacklist_file: Optional[str] = None) -> Set[Tuple[str, str]]:
    """
    Load position blacklist from file or return default blacklist.
    
    Args:
        blacklist_file: Path to blacklist file (one position per line: chr\tpos)
    
    Returns:
        Set of (chr_num, pos) tuples
    """
    if blacklist_file and os.path.exists(blacklist_file):
        blacklist = set()
        try:
            with open(blacklist_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    if len(fields) >= 2:
                        chr_col = fields[0]
                        pos_col = fields[1]
                        chr_num = extract_chr_number(chr_col)
                        blacklist.add((chr_num, pos_col))
            print(f"Loaded {len(blacklist)} positions from blacklist file: {blacklist_file}", file=sys.stderr)
            return blacklist
        except Exception as e:
            print(f"Warning: Error loading blacklist file {blacklist_file}: {e}. Using default blacklist.", file=sys.stderr)
            return DEFAULT_POSITION_BLACKLIST
    return DEFAULT_POSITION_BLACKLIST

def generate_intervals_file(tsv_file, output_file, variant_type='both', interval_prefix="FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.", chr_prefix="Chr", blacklist: Optional[Set[Tuple[str, str]]] = None):
    """
    Generate intervals file from TSV file.
    
    variant_type: 'snp', 'indel', or 'both' (default: 'both')
    Format: 1:1006357-1006357 (chromosome number without prefix)
    """
    if blacklist is None:
        blacklist = DEFAULT_POSITION_BLACKLIST
    
    intervals = []
    seen_intervals = set()  # Track (chr_num, pos) to avoid duplicates
    
    try:
        with open(tsv_file, 'r') as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chr_col = fields[0]
                    pos = fields[1]
                    var_type = fields[3]
                    
                    # Filter by variant type
                    if variant_type == 'both':
                        if var_type in ['snp', 'indel']:
                            include = True
                        else:
                            include = False
                    elif variant_type == 'snp':
                        include = (var_type == 'snp')
                    elif variant_type == 'indel':
                        include = (var_type == 'indel')
                    else:
                        include = False
                    
                    if include:
                        # Extract chromosome number (remove any prefix)
                        chr_num = extract_chr_number(chr_col)
                        
                        # Check blacklist
                        if (chr_num, pos) in blacklist:
                            continue
                        
                        # Check for duplicates using (chr_num, pos) tuple
                        interval_key = (chr_num, pos)
                        if interval_key not in seen_intervals:
                            seen_intervals.add(interval_key)
                            interval = f"{chr_num}:{pos}-{pos}"
                            intervals.append((chr_num, int(chr_num) if chr_num.isdigit() else 999999, int(pos), interval))
        
        # Sort by chr number (numeric) and pos
        intervals.sort(key=lambda x: (x[1], x[2]))
        
        # Write intervals file
        with open(output_file, 'w') as f:
            for _, _, _, interval in intervals:
                f.write(interval + '\n')
        
        print(f"Intervals file written to {output_file}", file=sys.stderr)
        print(f"Total intervals: {len(intervals)}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error generating intervals file: {e}", file=sys.stderr)

def extract_snp_from_vcf(snp_indel_vcf_files, output_file, chr_prefix="Chr", blacklist: Optional[Set[Tuple[str, str]]] = None):
    """
    Extract SNP variants from VCF files and generate snp.sites.vcf.
    
    snp_indel_vcf_files: list of (filename, real_path) tuples for SNP/INDEL VCF files
    """
    if blacklist is None:
        blacklist = DEFAULT_POSITION_BLACKLIST
    # VCF header for SNP sites
    file_date = datetime.now().strftime("%Y%m%d")
    vcf_header = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate={file_date}
##NOTE=This is a sites-only VCF file containing variant positions only. No sample genotype data is included. This file should not be used for genotype calling or analysis requiring sample data.
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
    
    snp_records = []
    seen_variants = set()  # Track (chr, pos) to avoid duplicates
    
    for filename, vcf_file in sorted(snp_indel_vcf_files):
        print(f"Extracting SNP from {filename}...", file=sys.stderr)
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip header lines
                    if line.startswith('##'):
                        continue
                    # Skip column header line
                    if line.startswith('#CHROM'):
                        continue
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Parse data line
                    fields = line.split('\t')
                    if len(fields) >= 8:
                        chr_col = fields[0]
                        pos_col = fields[1]
                        id_field = fields[2] if len(fields) > 2 else '.'
                        ref = fields[3] if len(fields) > 3 else 'N'
                        alt = fields[4] if len(fields) > 4 else '.'
                        qual = fields[5] if len(fields) > 5 else '.'
                        filter_field = fields[6] if len(fields) > 6 else '.'
                        info = fields[7] if len(fields) > 7 else '.'
                        
                        # Make biallelic (keep only first ALT allele)
                        ref, alt = make_biallelic(ref, alt)
                        
                        # Determine variant type from ID field (for RiceNavi) or INFO field
                        var_type = None
                        if id_field.startswith('SNP-'):
                            var_type = 'snp'
                        elif id_field.startswith('INDEL-'):
                            var_type = 'indel'
                        elif id_field.startswith('SV-'):
                            var_type = 'sv'
                        else:
                            # Try to infer from REF/ALT lengths
                            if abs(len(ref) - len(alt)) > 1:
                                var_type = 'indel'
                            elif len(ref) == 1 and len(alt) == 1:
                                var_type = 'snp'
                            # Check INFO field for SVTYPE
                            if 'SVTYPE=' in info:
                                var_type = 'sv'
                        
                        # Only include SNP variants, skip INDEL and SV
                        if var_type != 'snp':
                            continue
                        
                        # Extract chromosome number (remove any prefix)
                        chr_num = extract_chr_number(chr_col)
                        
                        # Check blacklist
                        if (chr_num, pos_col) in blacklist:
                            continue
                        
                        # Check for duplicates
                        pos = int(pos_col) if pos_col.isdigit() else 0
                        variant_key = (chr_num, pos)
                        if variant_key not in seen_variants:
                            seen_variants.add(variant_key)
                            snp_records.append({
                                'chr': chr_num,
                                'chr_num': int(chr_num) if chr_num.isdigit() else 999999,  # For sorting
                                'pos': pos_col,
                                'id': id_field,
                                'ref': ref,
                                'alt': alt,
                                'qual': qual,
                                'filter': filter_field,
                                'info': info
                            })
        except Exception as e:
            print(f"Error reading SNP VCF {vcf_file}: {e}", file=sys.stderr)
    
    # Sort by chr number (numeric) and pos
    snp_records.sort(key=lambda x: (x['chr_num'], int(x['pos']) if x['pos'].isdigit() else 0))
    
    # Write output VCF file
    try:
        with open(output_file, 'w') as f:
            f.write(vcf_header)
            for record in snp_records:
                f.write(f"{record['chr']}\t{record['pos']}\t{record['id']}\t{record['ref']}\t{record['alt']}\t{record['qual']}\t{record['filter']}\t{record['info']}\n")
        
        print(f"SNP sites VCF written to {output_file}", file=sys.stderr)
        print(f"Total SNP records: {len(snp_records)}", file=sys.stderr)
    except Exception as e:
        print(f"Error writing SNP VCF: {e}", file=sys.stderr)

def extract_indel_from_vcf(snp_indel_vcf_files, output_file, chr_prefix="Chr", blacklist: Optional[Set[Tuple[str, str]]] = None):
    """
    Extract INDEL variants from VCF files and generate indel.sites.vcf.
    
    snp_indel_vcf_files: list of (filename, real_path) tuples for SNP/INDEL VCF files
    """
    if blacklist is None:
        blacklist = DEFAULT_POSITION_BLACKLIST
    # VCF header for INDEL sites
    file_date = datetime.now().strftime("%Y%m%d")
    vcf_header = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate={file_date}
##NOTE=This is a sites-only VCF file containing variant positions only. No sample genotype data is included. This file should not be used for genotype calling or analysis requiring sample data.
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
    
    indel_records = []
    seen_variants = set()  # Track (chr, pos) to avoid duplicates
    
    for filename, vcf_file in sorted(snp_indel_vcf_files):
        print(f"Extracting INDEL from {filename}...", file=sys.stderr)
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip header lines
                    if line.startswith('##'):
                        continue
                    # Skip column header line
                    if line.startswith('#CHROM'):
                        continue
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Parse data line
                    fields = line.split('\t')
                    if len(fields) >= 8:
                        chr_col = fields[0]
                        pos_col = fields[1]
                        id_field = fields[2] if len(fields) > 2 else '.'
                        ref = fields[3] if len(fields) > 3 else 'N'
                        alt = fields[4] if len(fields) > 4 else '.'
                        qual = fields[5] if len(fields) > 5 else '.'
                        filter_field = fields[6] if len(fields) > 6 else '.'
                        info = fields[7] if len(fields) > 7 else '.'
                        
                        # Make biallelic (keep only first ALT allele)
                        ref, alt = make_biallelic(ref, alt)
                        
                        # Determine variant type from ID field (for RiceNavi) or INFO field
                        var_type = None
                        if id_field.startswith('SNP-'):
                            var_type = 'snp'
                        elif id_field.startswith('INDEL-'):
                            var_type = 'indel'
                        elif id_field.startswith('SV-'):
                            var_type = 'sv'
                        else:
                            # Try to infer from REF/ALT lengths
                            if abs(len(ref) - len(alt)) > 1:
                                var_type = 'indel'
                            elif len(ref) == 1 and len(alt) == 1:
                                var_type = 'snp'
                            # Check INFO field for SVTYPE
                            if 'SVTYPE=' in info:
                                var_type = 'sv'
                        
                        # Only include INDEL variants, skip SNP and SV
                        if var_type != 'indel':
                            continue
                        
                        # Extract chromosome number (remove any prefix)
                        chr_num = extract_chr_number(chr_col)
                        
                        # Check blacklist
                        if (chr_num, pos_col) in blacklist:
                            continue
                        
                        # Check for duplicates
                        pos = int(pos_col) if pos_col.isdigit() else 0
                        variant_key = (chr_num, pos)
                        if variant_key not in seen_variants:
                            seen_variants.add(variant_key)
                            indel_records.append({
                                'chr': chr_num,
                                'chr_num': int(chr_num) if chr_num.isdigit() else 999999,  # For sorting
                                'pos': pos_col,
                                'id': id_field,
                                'ref': ref,
                                'alt': alt,
                                'qual': qual,
                                'filter': filter_field,
                                'info': info
                            })
        except Exception as e:
            print(f"Error reading INDEL VCF {vcf_file}: {e}", file=sys.stderr)
    
    # Sort by chr number (numeric) and pos
    indel_records.sort(key=lambda x: (x['chr_num'], int(x['pos']) if x['pos'].isdigit() else 0))
    
    # Write output VCF file
    try:
        with open(output_file, 'w') as f:
            f.write(vcf_header)
            for record in indel_records:
                f.write(f"{record['chr']}\t{record['pos']}\t{record['id']}\t{record['ref']}\t{record['alt']}\t{record['qual']}\t{record['filter']}\t{record['info']}\n")
        
        print(f"INDEL sites VCF written to {output_file}", file=sys.stderr)
        print(f"Total INDEL records: {len(indel_records)}", file=sys.stderr)
    except Exception as e:
        print(f"Error writing INDEL VCF: {e}", file=sys.stderr)

def extract_sv_from_vcf(sv_vcf_files, output_file, chr_prefix="Chr", blacklist: Optional[Set[Tuple[str, str]]] = None):
    """
    Extract SV variants from VCF files and generate delly.sv.sites.vcf.
    
    sv_vcf_files: list of (filename, real_path) tuples for SV VCF files
    """
    if blacklist is None:
        blacklist = DEFAULT_POSITION_BLACKLIST
    # VCF header for delly SV sites
    file_date = datetime.now().strftime("%Y%m%d")
    delly_header = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate={file_date}
##NOTE=This is a sites-only VCF file containing variant positions only. No sample genotype data is included. This file should not be used for genotype calling or analysis requiring sample data.
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of PEs and SRs.">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2 coordinate in case of an inter-chromosomal translocation">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2 in case of an inter-chromosomal translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description="Median mapping quality of split-reads">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CONSBP,Number=1,Type=Integer,Description="Consensus SV breakpoint position">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
    
    sv_records = []
    seen_svs = set()  # Track (chr, pos) to avoid duplicates
    
    for filename, vcf_file in sorted(sv_vcf_files):
        print(f"Extracting SV from {filename}...", file=sys.stderr)
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip header lines
                    if line.startswith('##'):
                        continue
                    # Skip column header line
                    if line.startswith('#CHROM'):
                        continue
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Parse data line
                    fields = line.split('\t')
                    if len(fields) >= 8:
                        chr_col = fields[0]
                        pos_col = fields[1]
                        id_field = fields[2] if len(fields) > 2 else '.'
                        ref = fields[3] if len(fields) > 3 else 'N'
                        alt = fields[4] if len(fields) > 4 else '.'
                        qual = fields[5] if len(fields) > 5 else '.'
                        filter_field = fields[6] if len(fields) > 6 else '.'
                        info = fields[7] if len(fields) > 7 else '.'
                        
                        # Extract chromosome number (remove any prefix)
                        chr_num = extract_chr_number(chr_col)
                        
                        # Check blacklist
                        if (chr_num, pos_col) in blacklist:
                            continue
                        
                        # Check for duplicates
                        pos = int(pos_col) if pos_col.isdigit() else 0
                        sv_key = (chr_num, pos)
                        if sv_key not in seen_svs:
                            seen_svs.add(sv_key)
                            sv_records.append({
                                'chr': chr_num,
                                'chr_num': int(chr_num) if chr_num.isdigit() else 999999,  # For sorting
                                'pos': pos_col,
                                'id': id_field,
                                'ref': ref,
                                'alt': alt,
                                'qual': qual,
                                'filter': filter_field,
                                'info': info
                            })
        except Exception as e:
            print(f"Error reading SV VCF {vcf_file}: {e}", file=sys.stderr)
    
    # Sort by chr number (numeric) and pos
    sv_records.sort(key=lambda x: (x['chr_num'], int(x['pos']) if x['pos'].isdigit() else 0))
    
    # Write temporary VCF file, then convert to BCF
    
    # Determine output format based on file extension
    output_is_bcf = output_file.endswith('.bcf')
    temp_vcf_file = output_file.replace('.bcf', '.vcf') if output_is_bcf else output_file
    
    try:
        # Write VCF file
        with open(temp_vcf_file, 'w') as f:
            f.write(delly_header)
            for record in sv_records:
                f.write(f"{record['chr']}\t{record['pos']}\t{record['id']}\t{record['ref']}\t{record['alt']}\t{record['qual']}\t{record['filter']}\t{record['info']}\n")
        
        # Convert to BCF if needed
        if output_is_bcf:
            print(f"Converting VCF to BCF format...", file=sys.stderr)
            try:
                # Use bcftools to convert VCF to BCF
                result = subprocess.run(
                    ['bcftools', 'view', '-O', 'b', '-o', output_file, temp_vcf_file],
                    check=True,
                    capture_output=True,
                    text=True
                )
                # Remove temporary VCF file
                os.remove(temp_vcf_file)
                print(f"SV sites BCF written to {output_file}", file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error converting VCF to BCF: {e}", file=sys.stderr)
                print(f"bcftools stderr: {e.stderr}", file=sys.stderr)
                raise
            except FileNotFoundError:
                print("Error: bcftools not found. Please install bcftools to generate BCF files.", file=sys.stderr)
                raise
        else:
            print(f"SV sites VCF written to {temp_vcf_file}", file=sys.stderr)
        
        print(f"Total SV records: {len(sv_records)}", file=sys.stderr)
    except Exception as e:
        print(f"Error writing SV file: {e}", file=sys.stderr)
        # Clean up temporary file if it exists
        if output_is_bcf and os.path.exists(temp_vcf_file):
            try:
                os.remove(temp_vcf_file)
            except:
                pass
        raise

def validate_and_resolve_vcf_files(file_list: List[str]) -> List[Tuple[str, str]]:
    """
    Validate and resolve VCF file paths (including symlinks).
    
    Args:
        file_list: List of VCF file paths
    
    Returns:
        List of (filename, real_path) tuples
    """
    resolved_files = []
    seen_files = set()
    
    for file_path in file_list:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"VCF file not found: {file_path}")
        
        if not file_path.endswith('.vcf'):
            print(f"Warning: File does not have .vcf extension: {file_path}", file=sys.stderr)
        
        # Resolve symlink to get real path
        real_path = os.path.realpath(file_path)
        filename = os.path.basename(file_path)
        
        # Avoid duplicates
        if real_path not in seen_files:
            resolved_files.append((filename, real_path))
            seen_files.add(real_path)
        else:
            print(f"Warning: Skipping duplicate file: {file_path} (resolves to same path)", file=sys.stderr)
    
    return resolved_files

def process_vcf_files(vcf_files: List[Tuple[str, str]], input_dir: str) -> List[dict]:
    """
    Process VCF files and extract variant positions with metadata.
    
    Args:
        vcf_files: List of (filename, real_path) tuples
        input_dir: Input directory (for relative paths)
    
    Returns:
        List of dictionaries with variant information
    """
    all_data = []
    
    for filename, vcf_file in sorted(vcf_files):
        print(f"Processing {filename}...", file=sys.stderr)
        
        # Parse filename to get metadata (use original filename, not resolved path)
        source, var_type, population = parse_filename(filename)
        
        # Check if this is a RiceNavi file (needs per-variant type detection)
        is_rice_navi = (source == 'RiceNavi')
        
        # Extract positions (use real path for reading)
        positions = extract_vcf_positions(vcf_file, is_rice_navi=is_rice_navi)
        
        # Add to all_data
        for pos_tuple in positions:
            chr_col = pos_tuple[0]
            pos = pos_tuple[1]
            # For RiceNavi, use per-variant type; otherwise use filename-based type
            variant_type = pos_tuple[2] if is_rice_navi else var_type
            
            all_data.append({
                'chr': chr_col,
                'pos': pos,
                'source': source,
                'type': variant_type,
                'population': population
            })
    
    # Sort by chr and pos
    all_data.sort(key=lambda x: (x['chr'], x['pos']))
    return all_data

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Summarize VCF files: extract variant positions and generate summary files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process VCF files by specifying each type
  %(prog)s --snp-vcf file1.snp.vcf file2.snp.vcf --indel-vcf file1.indel.vcf --sv-vcf file1.sv.vcf -o output -p 705rice
  
  # Process with custom blacklist
  %(prog)s --snp-vcf *.snp.vcf --indel-vcf *.indel.vcf --sv-vcf *.sv.vcf -o output -p myprefix --blacklist blacklist.txt
  
  # Skip interval file generation
  %(prog)s --snp-vcf file.snp.vcf --indel-vcf file.indel.vcf -o output -p 705rice --no-intervals
        """
    )
    
    parser.add_argument(
        '--snp-vcf',
        type=str,
        nargs='+',
        default=[],
        metavar='FILE',
        help='SNP VCF file(s). Can specify multiple files.'
    )
    
    parser.add_argument(
        '--indel-vcf',
        type=str,
        nargs='+',
        default=[],
        metavar='FILE',
        help='INDEL VCF file(s). Can specify multiple files.'
    )
    
    parser.add_argument(
        '--sv-vcf',
        type=str,
        nargs='+',
        default=[],
        metavar='FILE',
        help='SV VCF file(s). Can specify multiple files.'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        default='.',
        help='Output directory for generated files (default: current directory)'
    )
    
    parser.add_argument(
        '-p', '--prefix',
        type=str,
        default='705rice',
        help='Prefix for output files (default: 705rice)'
    )
    
    parser.add_argument(
        '--blacklist',
        type=str,
        default=None,
        help='Path to blacklist file (one position per line: chr\\tpos). If not provided, uses default blacklist.'
    )
    
    parser.add_argument(
        '--interval-prefix',
        type=str,
        default='FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.',
        help='Prefix for interval names (default: FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.)'
    )
    
    parser.add_argument(
        '--chr-prefix',
        type=str,
        default='Chr',
        help='Chromosome prefix to remove (default: Chr)'
    )
    
    parser.add_argument(
        '--no-intervals',
        action='store_true',
        help='Skip generating interval files'
    )
    
    parser.add_argument(
        '--no-snp-vcf',
        action='store_true',
        help='Skip generating SNP sites VCF file'
    )
    
    parser.add_argument(
        '--no-indel-vcf',
        action='store_true',
        help='Skip generating INDEL sites VCF file'
    )
    
    parser.add_argument(
        '--no-sv-vcf',
        action='store_true',
        help='Skip generating SV sites VCF file'
    )
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_args()
    
    # Validate that at least one VCF file type is provided
    if not args.snp_vcf and not args.indel_vcf and not args.sv_vcf:
        print("Error: At least one VCF file type must be specified (--snp-vcf, --indel-vcf, or --sv-vcf)", file=sys.stderr)
        sys.exit(1)
    
    # Set up output directory
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Load blacklist
    blacklist = load_blacklist(args.blacklist)
    
    # Validate and resolve VCF files
    print("Validating VCF files...", file=sys.stderr)
    snp_vcf_files = []
    indel_vcf_files = []
    sv_vcf_files = []
    
    if args.snp_vcf:
        snp_vcf_files = validate_and_resolve_vcf_files(args.snp_vcf)
        print(f"  - SNP VCF files: {len(snp_vcf_files)}", file=sys.stderr)
    
    if args.indel_vcf:
        indel_vcf_files = validate_and_resolve_vcf_files(args.indel_vcf)
        print(f"  - INDEL VCF files: {len(indel_vcf_files)}", file=sys.stderr)
    
    if args.sv_vcf:
        sv_vcf_files = validate_and_resolve_vcf_files(args.sv_vcf)
        print(f"  - SV VCF files: {len(sv_vcf_files)}", file=sys.stderr)
    
    # Combine all VCF files for summary generation
    all_vcf_files = snp_vcf_files + indel_vcf_files + sv_vcf_files
    # For SNP/INDEL extraction, combine SNP and INDEL files
    snp_indel_vcf_files = snp_vcf_files + indel_vcf_files
    
    # Process VCF files and collect data
    print("\nProcessing VCF files...", file=sys.stderr)
    all_data = process_vcf_files(all_vcf_files, output_dir)
    
    if not all_data:
        print("Warning: No variant positions extracted from VCF files", file=sys.stderr)
        sys.exit(0)
    
    # Write TSV summary file
    tsv_output_file = os.path.join(output_dir, f'{args.prefix}.markers.tsv')
    with open(tsv_output_file, 'w') as f:
        # Write header
        f.write('chr\tpos\tsource\ttype\tpopulation\n')
        
        # Write data
        for row in all_data:
            f.write(f"{row['chr']}\t{row['pos']}\t{row['source']}\t{row['type']}\t{row['population']}\n")
    
    print(f"\nSummary written to: {tsv_output_file}", file=sys.stderr)
    print(f"Total records: {len(all_data)}", file=sys.stderr)
    
    # Generate interval files
    if not args.no_intervals:
        print("\nGenerating interval files...", file=sys.stderr)
        snp_intervals_file = os.path.join(output_dir, f'{args.prefix}.snp.markers.intervals')
        generate_intervals_file(
            tsv_output_file, snp_intervals_file,
            variant_type='snp',
            interval_prefix=args.interval_prefix,
            chr_prefix=args.chr_prefix,
            blacklist=blacklist
        )
        
        indel_intervals_file = os.path.join(output_dir, f'{args.prefix}.indel.markers.intervals')
        generate_intervals_file(
            tsv_output_file, indel_intervals_file,
            variant_type='indel',
            interval_prefix=args.interval_prefix,
            chr_prefix=args.chr_prefix,
            blacklist=blacklist
        )
    
    # Generate sites VCF files
    print("\nGenerating sites VCF files...", file=sys.stderr)
    
    if not args.no_snp_vcf and snp_vcf_files:
        snp_output_file = os.path.join(output_dir, f'{args.prefix}.snp.sites.vcf')
        extract_snp_from_vcf(snp_vcf_files, snp_output_file, chr_prefix=args.chr_prefix, blacklist=blacklist)
    elif args.no_snp_vcf:
        print("Skipping SNP sites VCF generation (--no-snp-vcf)", file=sys.stderr)
    elif not snp_vcf_files:
        print("No SNP VCF files provided, skipping SNP sites VCF generation", file=sys.stderr)
    
    if not args.no_indel_vcf and indel_vcf_files:
        indel_output_file = os.path.join(output_dir, f'{args.prefix}.indel.sites.vcf')
        extract_indel_from_vcf(indel_vcf_files, indel_output_file, chr_prefix=args.chr_prefix, blacklist=blacklist)
    elif args.no_indel_vcf:
        print("Skipping INDEL sites VCF generation (--no-indel-vcf)", file=sys.stderr)
    elif not indel_vcf_files:
        print("No INDEL VCF files provided, skipping INDEL sites VCF generation", file=sys.stderr)
    
    if not args.no_sv_vcf and sv_vcf_files:
        sv_output_file = os.path.join(output_dir, f'{args.prefix}.delly.sv.sites.bcf')
        extract_sv_from_vcf(sv_vcf_files, sv_output_file, chr_prefix=args.chr_prefix, blacklist=blacklist)
    elif args.no_sv_vcf:
        print("Skipping SV sites VCF generation (--no-sv-vcf)", file=sys.stderr)
    elif not sv_vcf_files:
        print("No SV VCF files provided, skipping SV sites VCF generation", file=sys.stderr)
    
    print("\nDone!", file=sys.stderr)

if __name__ == '__main__':
    main()

