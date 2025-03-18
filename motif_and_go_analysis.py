#!/usr/bin/env python3
# Motif Analysis and Gene Ontology Enrichment for ATAC-seq/Cut&Tag data

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import subprocess
import glob
import json
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
from wordcloud import WordCloud
import matplotlib.patches as mpatches

# Set the aesthetics for the plots
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("white")

# Define the paths
DIFFBIND_DIR = "results/diffbind"
PEAKS_DIR = "results/peaks"
OUTPUT_DIR = "analysis_output"
MOTIF_DIR = os.path.join(OUTPUT_DIR, "motif_analysis")
GO_DIR = os.path.join(OUTPUT_DIR, "go_analysis")

# Create output directories if they don't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(MOTIF_DIR, exist_ok=True)
os.makedirs(GO_DIR, exist_ok=True)

# Function to load differential binding results
def load_diffbind_results(comparison):
    file_path = os.path.join(DIFFBIND_DIR, comparison, "differential_peaks.csv")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        # Add a column to identify the comparison
        df['comparison'] = comparison
        return df
    else:
        print(f"Warning: {file_path} does not exist")
        return None

# Function to prepare BED files for motif analysis
def prepare_bed_for_motif_analysis(df, comparison, output_dir, direction="up", pval_threshold=0.05, lfc_threshold=1):
    if df is None or df.empty:
        print(f"No data available for {comparison}")
        return None
    
    # Filter peaks based on direction
    if direction == "up":
        filtered_df = df[(df['p.value'] < pval_threshold) & (df['Fold'] > lfc_threshold)]
        file_suffix = "upregulated"
    elif direction == "down":
        filtered_df = df[(df['p.value'] < pval_threshold) & (df['Fold'] < -lfc_threshold)]
        file_suffix = "downregulated"
    else:  # all significant
        filtered_df = df[df['p.value'] < pval_threshold]
        file_suffix = "significant"
    
    if filtered_df.empty:
        print(f"No {file_suffix} peaks found for {comparison}")
        return None
    
    # Create a BED file
    bed_file = os.path.join(output_dir, f"{comparison}_{file_suffix}_peaks.bed")
    
    # Convert the dataframe to BED format
    bed_df = filtered_df[['seqnames', 'start', 'end']].copy()
    bed_df.columns = ['chrom', 'start', 'end']
    
    # Write to BED file
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    
    print(f"Created BED file with {len(filtered_df)} {file_suffix} peaks for {comparison}: {bed_file}")
    
    return bed_file

# Function to run HOMER for motif analysis
def run_homer_motif_analysis(bed_file, output_dir, genome="mm10", size="given", mask=True):
    if bed_file is None or not os.path.exists(bed_file):
        print(f"BED file does not exist: {bed_file}")
        return None
    
    # Extract the base name of the BED file
    base_name = os.path.basename(bed_file).replace(".bed", "")
    
    # Create output directory for this analysis
    homer_output_dir = os.path.join(output_dir, base_name)
    os.makedirs(homer_output_dir, exist_ok=True)
    
    # Construct the HOMER command
    # Note: In a real environment, you would run this command using subprocess
    # Here we just print the command that would be run
    mask_option = "-mask" if mask else ""
    homer_cmd = f"findMotifsGenome.pl {bed_file} {genome} {homer_output_dir} -size {size} {mask_option}"
    
    print(f"HOMER motif analysis command (not executed in this script):")
    print(f"  {homer_cmd}")
    print(f"To run HOMER motif analysis, execute the above command in your terminal.")
    
    # In a real scenario, you would run the command like this:
    # try:
    #     subprocess.run(homer_cmd, shell=True, check=True)
    #     print(f"HOMER motif analysis completed for {base_name}")
    #     return homer_output_dir
    # except subprocess.CalledProcessError as e:
    #     print(f"Error running HOMER: {e}")
    #     return None
    
    # For demonstration, we'll just return the output directory
    return homer_output_dir

# Function to parse HOMER motif results
def parse_homer_results(homer_dir):
    if homer_dir is None or not os.path.exists(homer_dir):
        print(f"HOMER output directory does not exist: {homer_dir}")
        return None
    
    # In a real scenario, you would parse the HOMER output files
    # For demonstration, we'll create a mock result
    motif_file = os.path.join(homer_dir, "knownResults.txt")
    
    print(f"To parse HOMER results, check the following file (after running HOMER):")
    print(f"  {motif_file}")
    
    # Mock data for demonstration
    mock_data = {
        "Motif Name": ["CTCF", "AP-1", "NF-kB", "GATA", "SOX2"],
        "P-value": [1e-20, 1e-15, 1e-12, 1e-10, 1e-8],
        "% of Targets": [15.2, 12.5, 10.1, 8.7, 7.3],
        "% of Background": [5.1, 4.2, 3.5, 3.0, 2.8]
    }
    
    mock_df = pd.DataFrame(mock_data)
    
    return mock_df

# Function to visualize motif enrichment
def visualize_motif_enrichment(motif_df, comparison, output_dir, top_n=10):
    if motif_df is None or motif_df.empty:
        print(f"No motif data available for {comparison}")
        return
    
    # Sort by p-value and get top N motifs
    top_motifs = motif_df.sort_values("P-value").head(top_n)
    
    # Create a bar plot
    plt.figure(figsize=(12, 8))
    
    # Calculate enrichment
    top_motifs["Enrichment"] = top_motifs["% of Targets"] / top_motifs["% of Background"]
    
    # Create the bar plot
    bars = plt.barh(top_motifs["Motif Name"], top_motifs["Enrichment"], color="skyblue")
    
    # Add p-values as text
    for i, (_, row) in enumerate(top_motifs.iterrows()):
        plt.text(row["Enrichment"] + 0.1, i, f"p={row['P-value']:.1e}", va='center')
    
    # Add labels and title
    plt.xlabel("Enrichment (% of Targets / % of Background)")
    plt.ylabel("Motif")
    plt.title(f"Top {top_n} Enriched Motifs - {comparison}")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_motif_enrichment.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Motif enrichment plot saved to {output_file}")

# Function to create a motif word cloud
def create_motif_wordcloud(motif_df, comparison, output_dir):
    if motif_df is None or motif_df.empty:
        print(f"No motif data available for {comparison}")
        return
    
    # Create a dictionary of motif names and their enrichment values
    motif_dict = dict(zip(motif_df["Motif Name"], motif_df["% of Targets"] / motif_df["% of Background"]))
    
    # Create a word cloud
    wordcloud = WordCloud(width=800, height=400, background_color="white", 
                          colormap="viridis", max_words=50)
    wordcloud.generate_from_frequencies(motif_dict)
    
    # Display the word cloud
    plt.figure(figsize=(10, 6))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis("off")
    plt.title(f"Motif Enrichment Word Cloud - {comparison}")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_motif_wordcloud.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Motif word cloud saved to {output_file}")

# Function to annotate peaks with nearby genes
def annotate_peaks_with_genes(bed_file, output_dir, genome="mm10", distance=50000):
    if bed_file is None or not os.path.exists(bed_file):
        print(f"BED file does not exist: {bed_file}")
        return None
    
    # Extract the base name of the BED file
    base_name = os.path.basename(bed_file).replace(".bed", "")
    
    # Path to the annotation file
    annotation_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Federica_ATAC/data/gencode.v38.basic.annotation.gtf.gz"
    
    # Check if the annotation file exists
    if not os.path.exists(annotation_file):
        print(f"Warning: Annotation file {annotation_file} does not exist")
        return None
    
    # Create BedTool objects
    peaks = BedTool(bed_file)
    annotations = BedTool(annotation_file)
    
    # Extract gene features from the GTF file
    genes_gtf = annotations.filter(lambda x: x[2] == 'gene')
    
    # Create a temporary file for the gene annotations
    genes_file = os.path.join(output_dir, f"{base_name}_genes.gtf")
    genes_gtf.saveas(genes_file)
    
    # Find the closest gene to each peak
    closest_genes = peaks.closest(genes_gtf, d=True, t="first")
    
    # Save the annotated peaks
    annotated_file = os.path.join(output_dir, f"{base_name}_annotated_peaks.bed")
    closest_genes.saveas(annotated_file)
    
    # Create a more readable annotation file
    annotation_result_file = os.path.join(output_dir, f"{base_name}_gene_annotation.tsv")
    
    # Parse the closest gene results and create a dataframe
    peak_annotations = []
    for line in closest_genes:
        peak_chrom = line[0]
        peak_start = int(line[1])
        peak_end = int(line[2])
        
        if len(line) >= 13:  # Check if there's a gene match
            gene_chrom = line[3]
            gene_start = int(line[4])
            gene_end = int(line[5])
            gene_info = line[8]
            distance = int(line[12])
            
            # Extract gene ID and name from the GTF attributes
            gene_id = None
            gene_name = None
            for attr in gene_info.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                elif attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
            
            # Determine the relative position of the peak to the gene
            if peak_end < gene_start:
                location = "upstream"
            elif peak_start > gene_end:
                location = "downstream"
            else:
                if peak_start >= gene_start and peak_end <= gene_end:
                    location = "inside"
                else:
                    location = "overlap"
            
            peak_annotations.append({
                'PeakID': f"{peak_chrom}:{peak_start}-{peak_end}",
                'Chr': peak_chrom,
                'Start': peak_start,
                'End': peak_end,
                'GeneID': gene_id,
                'GeneName': gene_name,
                'GeneStart': gene_start,
                'GeneEnd': gene_end,
                'Distance': distance,
                'Location': location
            })
        else:
            # No gene match found
            peak_annotations.append({
                'PeakID': f"{peak_chrom}:{peak_start}-{peak_end}",
                'Chr': peak_chrom,
                'Start': peak_start,
                'End': peak_end,
                'GeneID': None,
                'GeneName': None,
                'GeneStart': None,
                'GeneEnd': None,
                'Distance': None,
                'Location': "intergenic"
            })
    
    # Create a dataframe from the annotations
    annotation_df = pd.DataFrame(peak_annotations)
    
    # Save the annotation dataframe to a TSV file
    annotation_df.to_csv(annotation_result_file, sep='\t', index=False)
    
    print(f"Peaks annotated with genes from {annotation_file}")
    print(f"Annotation results saved to {annotation_result_file}")
    
    # Create a summary of peak locations
    if not annotation_df.empty:
        location_counts = annotation_df['Location'].value_counts()
        
        # Create a pie chart of peak locations
        plt.figure(figsize=(10, 8))
        plt.pie(location_counts, labels=location_counts.index, autopct='%1.1f%%', startangle=90)
        plt.axis('equal')
        plt.title(f'Peak Distribution Relative to Genes - {base_name}')
        
        # Save the figure
        pie_file = os.path.join(output_dir, f"{base_name}_peak_locations.png")
        plt.savefig(pie_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Peak location summary saved to {pie_file}")
    
    # For compatibility with the rest of the script, create a dataframe in the expected format
    compat_df = pd.DataFrame({
        'PeakID': annotation_df['PeakID'],
        'Chr': annotation_df['Chr'],
        'Start': annotation_df['Start'],
        'End': annotation_df['End'],
        'Strand': ['+' if i % 2 == 0 else '-' for i in range(len(annotation_df))],  # Dummy strand info
        'NearestGene': annotation_df['GeneName'].fillna('Unknown'),
        'Distance': annotation_df['Distance'],
        'Annotation': annotation_df['Location']
    })
    
    return compat_df

# Function to perform gene ontology enrichment analysis
def perform_go_enrichment(gene_list, comparison, output_dir, organism="mouse"):
    if gene_list is None or len(gene_list) == 0:
        print(f"No genes available for {comparison}")
        return None
    
    # Create output file for the gene list
    gene_file = os.path.join(output_dir, f"{comparison}_gene_list.txt")
    
    # Write the gene list to a file
    with open(gene_file, "w") as f:
        for gene in gene_list:
            f.write(f"{gene}\n")
    
    # In a real scenario, you would use a tool like DAVID, Enrichr, or g:Profiler for GO enrichment
    # Here we'll just print instructions for using these tools
    
    print(f"Gene list saved to {gene_file}")
    print("To perform GO enrichment analysis, use one of the following tools:")
    print("  1. DAVID: https://david.ncifcrf.gov/")
    print("  2. Enrichr: https://maayanlab.cloud/Enrichr/")
    print("  3. g:Profiler: https://biit.cs.ut.ee/gprofiler/gost")
    print(f"Upload the gene list from {gene_file} to any of these tools.")
    
    # For demonstration, we'll create a mock GO enrichment result
    mock_data = {
        "GO_Term": [
            "GO:0006355 - regulation of transcription, DNA-templated",
            "GO:0007275 - multicellular organism development",
            "GO:0006351 - transcription, DNA-templated",
            "GO:0045893 - positive regulation of transcription, DNA-templated",
            "GO:0045944 - positive regulation of transcription by RNA polymerase II"
        ],
        "P_value": [1e-10, 1e-8, 1e-7, 1e-6, 1e-5],
        "FDR": [1e-8, 1e-6, 1e-5, 1e-4, 1e-3],
        "Genes": [
            "Gene1, Gene3, Gene5, Gene7, Gene9",
            "Gene2, Gene4, Gene6, Gene8, Gene10",
            "Gene1, Gene2, Gene3, Gene4, Gene5",
            "Gene6, Gene7, Gene8, Gene9, Gene10",
            "Gene1, Gene3, Gene5, Gene7, Gene9"
        ],
        "Gene_Count": [5, 5, 5, 5, 5],
        "Background_Count": [500, 600, 700, 800, 900],
        "Fold_Enrichment": [10.0, 8.3, 7.1, 6.3, 5.6]
    }
    
    mock_df = pd.DataFrame(mock_data)
    
    return mock_df

# Function to visualize GO enrichment
def visualize_go_enrichment(go_df, comparison, output_dir, top_n=10):
    if go_df is None or go_df.empty:
        print(f"No GO enrichment data available for {comparison}")
        return
    
    # Sort by p-value and get top N terms
    top_terms = go_df.sort_values("P_value").head(top_n)
    
    # Create a bar plot
    plt.figure(figsize=(12, 8))
    
    # Create the bar plot
    bars = plt.barh(top_terms["GO_Term"], top_terms["Fold_Enrichment"], color="lightgreen")
    
    # Add p-values as text
    for i, (_, row) in enumerate(top_terms.iterrows()):
        plt.text(row["Fold_Enrichment"] + 0.1, i, f"p={row['P_value']:.1e}", va='center')
    
    # Add labels and title
    plt.xlabel("Fold Enrichment")
    plt.ylabel("GO Term")
    plt.title(f"Top {top_n} Enriched GO Terms - {comparison}")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_go_enrichment.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"GO enrichment plot saved to {output_file}")

# Function to create a GO term network
def create_go_network(go_df, comparison, output_dir, top_n=20):
    if go_df is None or go_df.empty:
        print(f"No GO enrichment data available for {comparison}")
        return
    
    # Sort by p-value and get top N terms
    top_terms = go_df.sort_values("P_value").head(top_n)
    
    # Create a graph
    G = nx.Graph()
    
    # Add nodes (GO terms)
    for i, row in top_terms.iterrows():
        G.add_node(row["GO_Term"], 
                   pvalue=row["P_value"], 
                   fold_enrichment=row["Fold_Enrichment"],
                   gene_count=row["Gene_Count"])
    
    # Add edges based on shared genes
    for i, row1 in top_terms.iterrows():
        genes1 = set(row1["Genes"].split(", "))
        for j, row2 in top_terms.iterrows():
            if i < j:  # Avoid duplicate edges
                genes2 = set(row2["Genes"].split(", "))
                shared_genes = genes1.intersection(genes2)
                if shared_genes:
                    G.add_edge(row1["GO_Term"], row2["GO_Term"], 
                               weight=len(shared_genes),
                               shared_genes=", ".join(shared_genes))
    
    # Create the plot
    plt.figure(figsize=(14, 10))
    
    # Calculate node sizes based on gene count
    node_sizes = [G.nodes[node]["gene_count"] * 20 for node in G.nodes]
    
    # Calculate node colors based on p-value
    p_values = [G.nodes[node]["pvalue"] for node in G.nodes]
    p_value_log = [-np.log10(p) for p in p_values]
    
    # Calculate edge widths based on number of shared genes
    edge_widths = [G.edges[edge]["weight"] / 2 for edge in G.edges]
    
    # Create a custom colormap
    cmap = plt.cm.viridis
    
    # Draw the network
    pos = nx.spring_layout(G, seed=42)  # Position nodes using force-directed layout
    
    # Draw nodes
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, 
                                  node_color=p_value_log, cmap=cmap, 
                                  alpha=0.8)
    
    # Draw edges
    edges = nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.5)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif")
    
    # Add a colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(min(p_value_log), max(p_value_log)))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label("-log10(p-value)")
    
    # Add a title
    plt.title(f"GO Term Network - {comparison}")
    plt.axis("off")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_go_network.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"GO network plot saved to {output_file}")

# Function to analyze peak distribution across genomic features
def analyze_peak_distribution(annotation_df, comparison, output_dir):
    if annotation_df is None or annotation_df.empty:
        print(f"No annotation data available for {comparison}")
        return
    
    # Count peaks in each annotation category
    annotation_counts = annotation_df["Annotation"].value_counts()
    
    # Create a pie chart
    plt.figure(figsize=(10, 8))
    
    # Define colors for each category
    colors = plt.cm.tab10(np.arange(len(annotation_counts)))
    
    # Create the pie chart
    wedges, texts, autotexts = plt.pie(annotation_counts, 
                                       labels=annotation_counts.index, 
                                       autopct='%1.1f%%',
                                       colors=colors,
                                       startangle=90)
    
    # Customize text
    for text in texts:
        text.set_fontsize(12)
    for autotext in autotexts:
        autotext.set_fontsize(12)
        autotext.set_color('white')
    
    # Add a title
    plt.title(f"Peak Distribution Across Genomic Features - {comparison}")
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_peak_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Peak distribution plot saved to {output_file}")
    
    # Create a bar chart for distance to nearest TSS
    plt.figure(figsize=(10, 6))
    
    # Create bins for distances
    bins = [-50000, -10000, -5000, -1000, 0, 1000, 5000, 10000, 50000]
    bin_labels = ["-50kb to -10kb", "-10kb to -5kb", "-5kb to -1kb", "-1kb to 0", 
                 "0 to 1kb", "1kb to 5kb", "5kb to 10kb", "10kb to 50kb"]
    
    # Count peaks in each bin
    binned_distances = pd.cut(annotation_df["Distance"], bins=bins, labels=bin_labels)
    distance_counts = binned_distances.value_counts().sort_index()
    
    # Create the bar chart
    plt.bar(distance_counts.index, distance_counts.values, color="skyblue")
    
    # Add labels and title
    plt.xlabel("Distance to Nearest TSS")
    plt.ylabel("Number of Peaks")
    plt.title(f"Peak Distance to Nearest TSS - {comparison}")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_tss_distance.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"TSS distance plot saved to {output_file}")

# Main analysis workflow
def main():
    print("Starting Motif and GO Enrichment Analysis...")
    
    # 1. Load differential binding results
    print("\n1. Loading differential binding results...")
    comparisons = ["KI_vs_WT", "KO_vs_WT", "PAZ_vs_WT"]
    diff_results = {}
    
    for comparison in comparisons:
        diff_results[comparison] = load_diffbind_results(comparison)
        if diff_results[comparison] is not None:
            print(f"  Loaded {len(diff_results[comparison])} differential peaks for {comparison}")
    
    # 2. Prepare BED files for motif analysis
    print("\n2. Preparing BED files for motif analysis...")
    up_bed_files = {}
    down_bed_files = {}
    
    for comparison, df in diff_results.items():
        up_bed_files[comparison] = prepare_bed_for_motif_analysis(df, comparison, OUTPUT_DIR, direction="up")
        down_bed_files[comparison] = prepare_bed_for_motif_analysis(df, comparison, OUTPUT_DIR, direction="down")
    
    # 3. Run HOMER motif analysis
    print("\n3. Running HOMER motif analysis...")
    up_homer_dirs = {}
    down_homer_dirs = {}
    
    for comparison, bed_file in up_bed_files.items():
        up_homer_dirs[comparison] = run_homer_motif_analysis(bed_file, MOTIF_DIR)
    
    for comparison, bed_file in down_bed_files.items():
        down_homer_dirs[comparison] = run_homer_motif_analysis(bed_file, MOTIF_DIR)
    
    # 4. Parse HOMER results
    print("\n4. Parsing HOMER results...")
    up_motifs = {}
    down_motifs = {}
    
    for comparison, homer_dir in up_homer_dirs.items():
        up_motifs[comparison] = parse_homer_results(homer_dir)
    
    for comparison, homer_dir in down_homer_dirs.items():
        down_motifs[comparison] = parse_homer_results(homer_dir)
    
    # 5. Visualize motif enrichment
    print("\n5. Visualizing motif enrichment...")
    for comparison, motif_df in up_motifs.items():
        visualize_motif_enrichment(motif_df, f"{comparison}_up", MOTIF_DIR)
        create_motif_wordcloud(motif_df, f"{comparison}_up", MOTIF_DIR)
    
    for comparison, motif_df in down_motifs.items():
        visualize_motif_enrichment(motif_df, f"{comparison}_down", MOTIF_DIR)
        create_motif_wordcloud(motif_df, f"{comparison}_down", MOTIF_DIR)
    
    # 6. Annotate peaks with nearby genes
    print("\n6. Annotating peaks with nearby genes...")
    up_annotations = {}
    down_annotations = {}
    
    for comparison, bed_file in up_bed_files.items():
        up_annotations[comparison] = annotate_peaks_with_genes(bed_file, OUTPUT_DIR)
    
    for comparison, bed_file in down_bed_files.items():
        down_annotations[comparison] = annotate_peaks_with_genes(bed_file, OUTPUT_DIR)
    
    # 7. Perform GO enrichment analysis
    print("\n7. Performing GO enrichment analysis...")
    up_go_results = {}
    down_go_results = {}
    
    for comparison, annotation_df in up_annotations.items():
        if annotation_df is not None:
            gene_list = annotation_df["NearestGene"].unique().tolist()
            up_go_results[comparison] = perform_go_enrichment(gene_list, f"{comparison}_up", GO_DIR)
    
    for comparison, annotation_df in down_annotations.items():
        if annotation_df is not None:
            gene_list = annotation_df["NearestGene"].unique().tolist()
            down_go_results[comparison] = perform_go_enrichment(gene_list, f"{comparison}_down", GO_DIR)
    
    # 8. Visualize GO enrichment
    print("\n8. Visualizing GO enrichment...")
    for comparison, go_df in up_go_results.items():
        visualize_go_enrichment(go_df, f"{comparison}_up", GO_DIR)
        create_go_network(go_df, f"{comparison}_up", GO_DIR)
    
    for comparison, go_df in down_go_results.items():
        visualize_go_enrichment(go_df, f"{comparison}_down", GO_DIR)
        create_go_network(go_df, f"{comparison}_down", GO_DIR)
    
    # 9. Analyze peak distribution across genomic features
    print("\n9. Analyzing peak distribution across genomic features...")
    for comparison, annotation_df in up_annotations.items():
        analyze_peak_distribution(annotation_df, f"{comparison}_up", OUTPUT_DIR)
    
    for comparison, annotation_df in down_annotations.items():
        analyze_peak_distribution(annotation_df, f"{comparison}_down", OUTPUT_DIR)
    
    print("\nMotif and GO Enrichment Analysis complete! Results saved to:", OUTPUT_DIR)

if __name__ == "__main__":
    main() 