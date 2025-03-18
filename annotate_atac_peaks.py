#!/usr/bin/env python3

import os
import pandas as pd
import gzip
import re
import glob
from collections import defaultdict
import argparse

def parse_gtf(gtf_file):
    """Parse GTF file and extract gene information."""
    print(f"Parsing GTF file: {gtf_file}")
    genes = defaultdict(dict)
    promoters = defaultdict(list)
    exons = defaultdict(list)
    
    # Define promoter region (e.g., 2kb upstream of TSS)
    promoter_size = 2000
    
    with gzip.open(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            # Convert to 0-based coordinates for easier comparison
            start = int(start) - 1
            end = int(end)
            
            # Extract gene_id and gene_name
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
            
            if not gene_id_match:
                continue
                
            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1) if gene_name_match else gene_id
            
            # Store gene information
            if feature == 'gene':
                genes[gene_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'name': gene_name
                }
                
                # Define promoter region based on strand
                if strand == '+':
                    promoter_start = max(0, start - promoter_size)
                    promoter_end = start + 100  # Include a bit downstream of TSS
                    promoters[chrom].append((promoter_start, promoter_end, gene_id, gene_name))
                else:  # strand == '-'
                    promoter_start = end - 100  # Include a bit upstream of TSS
                    promoter_end = end + promoter_size
                    promoters[chrom].append((promoter_start, promoter_end, gene_id, gene_name))
            
            # Store exon information
            elif feature == 'exon':
                exons[chrom].append((start, end, gene_id, gene_name))
    
    print(f"Parsed {len(genes)} genes")
    return genes, promoters, exons

def annotate_peak(peak, promoters, exons, genes):
    """Annotate a single peak with genomic features."""
    chrom, start, end = peak['seqnames'], int(peak['start']), int(peak['end'])
    
    # Initialize annotation
    annotation = {
        'nearest_gene': '',
        'nearest_gene_distance': float('inf'),
        'feature': 'intergenic',
        'overlapping_genes': []
    }
    
    # Check if peak overlaps with promoters
    if chrom in promoters:
        for prom_start, prom_end, gene_id, gene_name in promoters[chrom]:
            if start <= prom_end and end >= prom_start:
                annotation['feature'] = 'promoter'
                annotation['overlapping_genes'].append(gene_name)
                
                # Calculate distance to TSS
                gene = genes[gene_id]
                tss = gene['start'] if gene['strand'] == '+' else gene['end']
                peak_center = (start + end) // 2
                distance = abs(peak_center - tss)
                
                if distance < annotation['nearest_gene_distance']:
                    annotation['nearest_gene'] = gene_name
                    annotation['nearest_gene_distance'] = distance
    
    # If not in promoter, check if peak overlaps with exons
    if annotation['feature'] == 'intergenic' and chrom in exons:
        for exon_start, exon_end, gene_id, gene_name in exons[chrom]:
            if start <= exon_end and end >= exon_start:
                annotation['feature'] = 'exon'
                annotation['overlapping_genes'].append(gene_name)
                
                # Calculate distance to gene
                gene = genes[gene_id]
                gene_center = (gene['start'] + gene['end']) // 2
                peak_center = (start + end) // 2
                distance = abs(peak_center - gene_center)
                
                if distance < annotation['nearest_gene_distance']:
                    annotation['nearest_gene'] = gene_name
                    annotation['nearest_gene_distance'] = distance
    
    # If not in promoter or exon, find nearest gene
    if annotation['feature'] == 'intergenic':
        for gene_id, gene in genes.items():
            if gene['chrom'] != chrom:
                continue
                
            # Check if peak is within gene body (intron)
            if start <= gene['end'] and end >= gene['start']:
                annotation['feature'] = 'intron'
                annotation['overlapping_genes'].append(gene['name'])
                
                # Calculate distance to gene center
                gene_center = (gene['start'] + gene['end']) // 2
                peak_center = (start + end) // 2
                distance = abs(peak_center - gene_center)
                
                if distance < annotation['nearest_gene_distance']:
                    annotation['nearest_gene'] = gene['name']
                    annotation['nearest_gene_distance'] = distance
            else:
                # Calculate distance to gene
                peak_center = (start + end) // 2
                if peak_center < gene['start']:
                    distance = gene['start'] - peak_center
                else:
                    distance = peak_center - gene['end']
                
                if distance < annotation['nearest_gene_distance']:
                    annotation['nearest_gene'] = gene['name']
                    annotation['nearest_gene_distance'] = distance
    
    # Remove duplicates from overlapping genes
    annotation['overlapping_genes'] = list(set(annotation['overlapping_genes']))
    
    return annotation

def process_peak_file(peak_file, promoters, exons, genes, output_dir):
    """Process a single peak file and annotate all peaks."""
    print(f"Processing peak file: {peak_file}")
    
    # Read peak file
    peaks = pd.read_csv(peak_file, sep='\t')
    
    # Annotate each peak
    annotations = []
    for _, peak in peaks.iterrows():
        annotation = annotate_peak(peak, promoters, exons, genes)
        annotations.append(annotation)
    
    # Add annotations to peaks dataframe
    for key in annotations[0].keys():
        if key == 'overlapping_genes':
            peaks[key] = [','.join(a[key]) if a[key] else '' for a in annotations]
        else:
            peaks[key] = [a[key] for a in annotations]
    
    # Create output filename
    basename = os.path.basename(peak_file)
    output_file = os.path.join(output_dir, f"annotated_{basename}")
    
    # Save annotated peaks
    peaks.to_csv(output_file, sep='\t', index=False)
    print(f"Saved annotated peaks to: {output_file}")
    
    return peaks

def main():
    parser = argparse.ArgumentParser(description='Annotate ATAC-seq peaks with genomic features')
    parser.add_argument('--gtf', required=True, help='Path to GTF annotation file (gzipped)')
    parser.add_argument('--peaks_dir', required=True, help='Directory containing narrowPeak files')
    parser.add_argument('--output_dir', default='results/annotated_peaks', help='Output directory for annotated peaks')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse GTF file
    genes, promoters, exons = parse_gtf(args.gtf)
    
    # Get all peak files
    peak_files = glob.glob(os.path.join(args.peaks_dir, '*.narrowPeak'))
    
    # Process each peak file
    all_annotated_peaks = []
    for peak_file in peak_files:
        annotated_peaks = process_peak_file(peak_file, promoters, exons, genes, args.output_dir)
        all_annotated_peaks.append(annotated_peaks)
    
    # Optionally, create a merged file with all annotated peaks
    if all_annotated_peaks:
        merged_peaks = pd.concat(all_annotated_peaks)
        merged_output = os.path.join(args.output_dir, "all_annotated_peaks.tsv")
        merged_peaks.to_csv(merged_output, sep='\t', index=False)
        print(f"Saved merged annotated peaks to: {merged_output}")
    
    print("Peak annotation completed successfully!")

if __name__ == "__main__":
    main() 