#!/usr/bin/env python3
# ATAC-seq/Cut&Tag Analysis Script

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pybedtools
from pybedtools import BedTool
import glob
from matplotlib_venn import venn2, venn3
import matplotlib.patches as mpatches
from adjustText import adjust_text
from upsetplot import from_memberships, plot

# Set the aesthetics for the plots
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("white")

# Define the paths
DIFFBIND_DIR = "results/diffbind"
PEAKS_DIR = "results/peaks"
OUTPUT_DIR = "analysis_output"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

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

# Function to load peak files
def load_peak_files(sample_name):
    file_path = os.path.join(PEAKS_DIR, f"{sample_name}_peaks.narrowPeak")
    if os.path.exists(file_path):
        # narrowPeak format: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak
        cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                'signalValue', 'pValue', 'qValue', 'peak']
        df = pd.read_csv(file_path, sep='\t', header=None, names=cols)
        return df
    else:
        print(f"Warning: {file_path} does not exist")
        return None

# Function to create volcano plot
def create_volcano_plot(df, comparison, output_dir, pval_threshold=0.05, lfc_threshold=1):
    if df is None or df.empty:
        print(f"No data available for {comparison}")
        return
    
    # Create a new figure
    plt.figure(figsize=(10, 8))
    
    # Calculate -log10 of p-values
    df['-log10(p-value)'] = -np.log10(df['p.value'])
    
    # Determine colors based on significance and fold change
    df['color'] = 'grey'  # Default color
    df.loc[(df['p.value'] < pval_threshold) & (df['Fold'] > lfc_threshold), 'color'] = 'red'  # Upregulated
    df.loc[(df['p.value'] < pval_threshold) & (df['Fold'] < -lfc_threshold), 'color'] = 'blue'  # Downregulated
    
    # Create the scatter plot
    plt.scatter(df['Fold'], df['-log10(p-value)'], c=df['color'], alpha=0.6, s=30)
    
    # Add lines for thresholds
    plt.axhline(-np.log10(pval_threshold), color='gray', linestyle='--')
    plt.axvline(lfc_threshold, color='gray', linestyle='--')
    plt.axvline(-lfc_threshold, color='gray', linestyle='--')
    
    # Add labels for some of the most significant points
    texts = []
    for idx, row in df.sort_values('p.value').head(10).iterrows():
        # Create a label with chromosome and position
        label = f"{row['seqnames']}:{row['start']}-{row['end']}"
        texts.append(plt.text(row['Fold'], -np.log10(row['p.value']), label, fontsize=8))
    
    # Adjust text to avoid overlapping
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))
    
    # Add legend
    red_patch = mpatches.Patch(color='red', label='Upregulated')
    blue_patch = mpatches.Patch(color='blue', label='Downregulated')
    grey_patch = mpatches.Patch(color='grey', label='Not significant')
    plt.legend(handles=[red_patch, blue_patch, grey_patch], loc='upper right')
    
    # Add labels and title
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.title(f'Volcano Plot - {comparison}')
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_volcano_plot.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Volcano plot saved to {output_file}")
    
    # Return summary statistics
    significant = df[df['p.value'] < pval_threshold]
    upregulated = significant[significant['Fold'] > lfc_threshold]
    downregulated = significant[significant['Fold'] < -lfc_threshold]
    
    return {
        'total': len(df),
        'significant': len(significant),
        'upregulated': len(upregulated),
        'downregulated': len(downregulated)
    }

# Function to create MA plot
def create_ma_plot(df, comparison, output_dir):
    if df is None or df.empty:
        print(f"No data available for {comparison}")
        return
    
    plt.figure(figsize=(10, 8))
    
    # Calculate average concentration
    df['AveExpr'] = (df['Conc_WT'] + df[f"Conc_{comparison.split('_')[0]}"]) / 2
    
    # Determine colors based on significance
    df['color'] = 'grey'  # Default color
    df.loc[df['FDR'] < 0.05, 'color'] = 'red'  # Significant
    
    # Create the scatter plot
    plt.scatter(df['AveExpr'], df['Fold'], c=df['color'], alpha=0.6, s=30)
    
    # Add a horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add labels and title
    plt.xlabel('Average Expression')
    plt.ylabel('Log2 Fold Change')
    plt.title(f'MA Plot - {comparison}')
    
    # Add legend
    red_patch = mpatches.Patch(color='red', label='FDR < 0.05')
    grey_patch = mpatches.Patch(color='grey', label='Not significant')
    plt.legend(handles=[red_patch, grey_patch], loc='upper right')
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_ma_plot.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"MA plot saved to {output_file}")

# Function to create heatmap of top differential peaks
def create_heatmap(df, comparison, output_dir, top_n=50):
    if df is None or df.empty:
        print(f"No data available for {comparison}")
        return
    
    # Get the top N differential peaks by p-value
    top_peaks = df.sort_values('p.value').head(top_n)
    
    # Create a heatmap of the fold changes
    plt.figure(figsize=(10, 12))
    
    # Create a new dataframe for the heatmap
    heatmap_data = pd.DataFrame({
        'Peak': [f"{row['seqnames']}:{row['start']}-{row['end']}" for _, row in top_peaks.iterrows()],
        'Log2FC': top_peaks['Fold'],
        '-log10(p-value)': -np.log10(top_peaks['p.value'])
    })
    
    # Set the peak as the index
    heatmap_data = heatmap_data.set_index('Peak')
    
    # Create the heatmap
    sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, annot=False, 
                linewidths=0.5, yticklabels=True)
    
    plt.title(f'Top {top_n} Differential Peaks - {comparison}')
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{comparison}_top_peaks_heatmap.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Heatmap saved to {output_file}")

    
# Function to compare peaks between conditions
def compare_peaks(peak_files, output_dir):
    # Load all peak files
    peak_dfs = {}
    for sample in peak_files:
        peak_dfs[sample] = load_peak_files(sample)
    
    # Create BED files for each sample
    bed_files = {}
    for sample, df in peak_dfs.items():
        if df is not None and not df.empty:
            bed_file = os.path.join(output_dir, f"{sample}_peaks.bed")
            df[['chrom', 'start', 'end']].to_csv(bed_file, sep='\t', header=False, index=False)
            bed_files[sample] = bed_file
    
    # Create BedTool objects
    bed_tools = {sample: BedTool(file) for sample, file in bed_files.items()}
    
    # Compare peaks between samples (pairwise)
    for i, (sample1, bed1) in enumerate(bed_tools.items()):
        for sample2, bed2 in list(bed_tools.items())[i+1:]:
            # Find overlapping peaks
            overlap = bed1.intersect(bed2, wa=True)
            
            # Count overlapping peaks
            overlap_count = len(overlap)
            
            # Calculate Jaccard index
            jaccard = bed1.jaccard(bed2)
            
            print(f"Comparison between {sample1} and {sample2}:")
            print(f"  Overlapping peaks: {overlap_count}")
            print(f"  Jaccard index: {jaccard['jaccard']}")
    
    # Create Venn diagrams for groups
    # Group samples by condition
    conditions = {}
    for sample in peak_files:
        condition = sample.split('-')[0]  # Assuming format like "WT1-ATAC"
        if condition not in conditions:
            conditions[condition] = []
        conditions[condition].append(sample)
    
    # Create a Venn diagram for each pair of conditions
    for i, (cond1, samples1) in enumerate(conditions.items()):
        for cond2, samples2 in list(conditions.items())[i+1:]:
            # Merge peaks within each condition
            merged_bed1 = None
            for sample in samples1:
                if sample in bed_tools:
                    if merged_bed1 is None:
                        merged_bed1 = bed_tools[sample]
                    else:
                        merged_bed1 = merged_bed1.cat(bed_tools[sample], postmerge=True)
            
            merged_bed2 = None
            for sample in samples2:
                if sample in bed_tools:
                    if merged_bed2 is None:
                        merged_bed2 = bed_tools[sample]
                    else:
                        merged_bed2 = merged_bed2.cat(bed_tools[sample], postmerge=True)
            
            if merged_bed1 is not None and merged_bed2 is not None:
                # Count peaks in each set
                count1 = len(merged_bed1)
                count2 = len(merged_bed2)
                
                # Count overlapping peaks
                overlap = merged_bed1.intersect(merged_bed2, wa=True)
                overlap_count = len(overlap)
                
                # Create Venn diagram
                plt.figure(figsize=(8, 8))
                venn2(subsets=(count1 - overlap_count, count2 - overlap_count, overlap_count),
                      set_labels=(cond1, cond2))
                plt.title(f'Peak Overlap between {cond1} and {cond2}')
                
                # Save the figure
                output_file = os.path.join(output_dir, f"{cond1}_vs_{cond2}_venn.png")
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"Venn diagram saved to {output_file}")

# Function to create an UpSet plot for peak overlaps
def create_upset_plot(peak_files, output_dir):
    # Load all peak files
    peak_dfs = {}
    for sample in peak_files:
        peak_dfs[sample] = load_peak_files(sample)
    
    # Create BED files for each sample
    bed_files = {}
    for sample, df in peak_dfs.items():
        if df is not None and not df.empty:
            bed_file = os.path.join(output_dir, f"{sample}_peaks.bed")
            df[['chrom', 'start', 'end']].to_csv(bed_file, sep='\t', header=False, index=False)
            bed_files[sample] = bed_file
    
    # Create BedTool objects
    bed_tools = {sample: BedTool(file) for sample, file in bed_files.items()}
    
    # Group samples by condition
    conditions = {}
    for sample in peak_files:
        condition = sample.split('-')[0]  # Assuming format like "WT1-ATAC"
        if condition not in conditions:
            conditions[condition] = []
        conditions[condition].append(sample)
    
    # Merge peaks within each condition
    merged_beds = {}
    for cond, samples in conditions.items():
        merged_bed = None
        for sample in samples:
            if sample in bed_tools:
                if merged_bed is None:
                    merged_bed = bed_tools[sample]
                else:
                    merged_bed = merged_bed.cat(bed_tools[sample], postmerge=True)
        if merged_bed is not None:
            merged_beds[cond] = merged_bed
    
    # Create a universe of all peaks
    all_peaks = None
    for bed in merged_beds.values():
        if all_peaks is None:
            all_peaks = bed
        else:
            all_peaks = all_peaks.cat(bed, postmerge=True)
    
    if all_peaks is None:
        print("No peaks found in any sample")
        return
    
    # Create a dictionary to store peak memberships
    memberships = {cond: set() for cond in merged_beds.keys()}
    
    # For each peak in the universe, check which conditions it belongs to
    for i, peak in enumerate(all_peaks):
        peak_str = f"{peak.chrom}:{peak.start}-{peak.end}"
        for cond, bed in merged_beds.items():
            if len(bed.intersect(BedTool([peak]), u=True)) > 0:
                memberships[cond].add(peak_str)
    
    # Create the UpSet plot
    plt.figure(figsize=(12, 8))
    upset_data = from_memberships(memberships)
    plot(upset_data, sort_by='cardinality')
    plt.title('Peak Overlap between Conditions')
    
    # Save the figure
    output_file = os.path.join(output_dir, "peak_overlap_upset.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"UpSet plot saved to {output_file}")

# Function to create a PCA plot of samples
def create_pca_plot(peak_files, output_dir):
    """
    Create a PCA plot of samples based on peak intensity data.
    
    Args:
        peak_files: List of sample names for which peak files exist
        output_dir: Directory to save the output plot
    """
    print("Creating PCA plot of samples based on peak intensity...")
    
    # Create a count matrix from the peak files
    # First, create a union of all peaks
    all_peaks = None
    sample_peak_dfs = {}
    
    # Load all peak files
    for sample in peak_files:
        df = load_peak_files(sample)
        if df is not None and not df.empty:
            sample_peak_dfs[sample] = df
            
            # Convert to BedTool
            bed_file = os.path.join(output_dir, f"{sample}_peaks.bed")
            df[['chrom', 'start', 'end']].to_csv(bed_file, sep='\t', header=False, index=False)
            
            if all_peaks is None:
                all_peaks = pybedtools.BedTool(bed_file)
            else:
                all_peaks = all_peaks.cat(pybedtools.BedTool(bed_file), postmerge=True)
    
    if all_peaks is None:
        print("No peaks found in any sample")
        return
    
    # Create a merged peak file
    merged_peaks_file = os.path.join(output_dir, "merged_peaks.bed")
    all_peaks = all_peaks.sort().merge()
    all_peaks.saveas(merged_peaks_file)
    
    # Create a count matrix
    count_matrix = pd.DataFrame(index=range(len(all_peaks)))
    
    # Add coordinates to the count matrix
    coords = []
    for i, peak in enumerate(all_peaks):
        coords.append(f"{peak.chrom}:{peak.start}-{peak.end}")
    count_matrix['coordinates'] = coords
    count_matrix = count_matrix.set_index('coordinates')
    
    # For each sample, count reads in each peak
    for sample in peak_files:
        if sample in sample_peak_dfs:
            # Get the filtered BAM file path
            bam_file = f"results/bam/{sample}.filtered.bam"
            
            if os.path.exists(bam_file):
                # Count reads in each peak
                counts = []
                bam = pybedtools.BedTool(bam_file)
                
                for peak in all_peaks:
                    peak_bed = pybedtools.BedTool([peak])
                    overlap = bam.intersect(peak_bed, c=True)
                    count = int(str(overlap[0]).strip().split('\t')[-1])
                    counts.append(count)
                
                count_matrix[sample] = counts
            else:
                print(f"Warning: BAM file {bam_file} not found")
    
    # Save the count matrix
    count_matrix.to_csv(os.path.join(output_dir, "peak_count_matrix.csv"))
    
    # Perform PCA
    if len(count_matrix.columns) < 2:
        print("Not enough samples for PCA")
        return
    
    # Normalize counts (log2 transform with pseudocount)
    normalized_counts = np.log2(count_matrix.values + 1)
    
    # Run PCA
    pca = stats.PCA(n_components=2)
    pca_result = pca.fit_transform(normalized_counts.T)
    
    # Create a dataframe with PCA results
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df['sample'] = count_matrix.columns
    
    # Add condition information
    pca_df['condition'] = [sample.split('-')[0] for sample in pca_df['sample']]
    
    # Create the PCA plot
    plt.figure(figsize=(10, 8))
    
    # Create a color palette for conditions
    conditions = pca_df['condition'].unique()
    palette = sns.color_palette("Set1", n_colors=len(conditions))
    color_dict = dict(zip(conditions, palette))
    
    # Plot each point
    for condition in conditions:
        subset = pca_df[pca_df['condition'] == condition]
        plt.scatter(subset['PC1'], subset['PC2'], 
                   label=condition, 
                   color=color_dict[condition],
                   s=100, alpha=0.7)
    
    # Add sample labels
    for i, row in pca_df.iterrows():
        plt.annotate(row['sample'], 
                    (row['PC1'], row['PC2']),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=10)
    
    # Add axis labels with variance explained
    variance_explained = pca.explained_variance_ratio_ * 100
    plt.xlabel(f'PC1 ({variance_explained[0]:.2f}%)')
    plt.ylabel(f'PC2 ({variance_explained[1]:.2f}%)')
    
    # Add legend and title
    plt.legend(title='Condition')
    plt.title('PCA of ATAC-seq Samples')
    plt.grid(True, alpha=0.3)
    
    # Add a tight layout
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "pca_plot.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"PCA plot saved to {output_file}")
    
    # Also save the PCA data
    pca_df.to_csv(os.path.join(output_dir, "pca_data.csv"), index=False)
    
    return pca_df

# Main analysis workflow
def main():
    print("Starting ATAC-seq/Cut&Tag analysis...")
    
    # 1. Load differential binding results
    print("\n1. Loading differential binding results...")
    comparisons = ["KI_vs_WT", "KO_vs_WT", "PAZ_vs_WT"]
    diff_results = {}
    
    for comparison in comparisons:
        diff_results[comparison] = load_diffbind_results(comparison)
        if diff_results[comparison] is not None:
            print(f"  Loaded {len(diff_results[comparison])} differential peaks for {comparison}")
    
    # 2. Create volcano plots
    print("\n2. Creating volcano plots...")
    volcano_stats = {}
    for comparison, df in diff_results.items():
        volcano_stats[comparison] = create_volcano_plot(df, comparison, OUTPUT_DIR)
        if volcano_stats[comparison]:
            print(f"  {comparison}: {volcano_stats[comparison]['significant']} significant peaks " +
                  f"({volcano_stats[comparison]['upregulated']} up, {volcano_stats[comparison]['downregulated']} down)")
    
    # 3. Create MA plots
    print("\n3. Creating MA plots...")
    for comparison, df in diff_results.items():
        create_ma_plot(df, comparison, OUTPUT_DIR)
    
    # 4. Create heatmaps of top differential peaks
    print("\n4. Creating heatmaps of top differential peaks...")
    for comparison, df in diff_results.items():
        create_heatmap(df, comparison, OUTPUT_DIR)
    
    # 5. Annotate peaks with genomic features
    print("\n5. Annotating peaks with genomic features...")
    for comparison, df in diff_results.items():
        annotate_peaks(df, comparison, OUTPUT_DIR)
    
    # 6. Compare peaks between conditions
    print("\n6. Comparing peaks between conditions...")
    # Get all peak files
    peak_files = [os.path.basename(f).replace("_peaks.narrowPeak", "") 
                 for f in glob.glob(os.path.join(PEAKS_DIR, "*_peaks.narrowPeak"))]
    
    compare_peaks(peak_files, OUTPUT_DIR)
    
    # 7. Create UpSet plot for peak overlaps
    print("\n7. Creating UpSet plot for peak overlaps...")
    create_upset_plot(peak_files, OUTPUT_DIR)
    
    # 8. Create PCA plot
    print("\n8. Creating PCA plot...")
    create_pca_plot(peak_files, OUTPUT_DIR)
    
    print("\nAnalysis complete! Results saved to:", OUTPUT_DIR)

if __name__ == "__main__":
    main() 