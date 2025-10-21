#!/usr/bin/env python3

"""
Script to filter and match variants from VA Million Veteran Program fine-mapping data.

This script:
1. Filters variants for EUR population
2. Extracts highest PIP per variant (identified by MVP ID)
3. Selects variants with max PIP >= 0.95
4. Matches variants with max PIP < 0.05 to the high PIP set based on:
   - Minor allele frequency distribution
   - P-value distribution
   - VEP annotation distribution
5. Outputs both sets to a single file with specified columns

Usage:
    python filter_and_match_variants.py -i Data_S1.xlsx -o matched_variants_output.csv
    python filter_and_match_variants.py --input Data_S1.xlsx --output results.csv
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import argparse
import sys
import warnings
warnings.filterwarnings('ignore')

def calculate_maf(eaf):
    """Convert effect allele frequency to minor allele frequency."""
    return np.where(eaf > 0.5, 1 - eaf, eaf)

def match_distributions(high_pip_df, low_pip_df, n_neighbors=5):
    """
    Match low PIP variants to high PIP variants based on MAF, p-value, and VEP annotation.
    
    Parameters:
    -----------
    high_pip_df : DataFrame
        Variants with max PIP >= 0.95
    low_pip_df : DataFrame
        Candidate variants with max PIP < 0.05
    n_neighbors : int
        Number of neighbors to consider for matching
    
    Returns:
    --------
    DataFrame : Matched low PIP variants
    """
    # Create MAF column
    high_pip_df['MAF'] = calculate_maf(high_pip_df['EAF Population'])
    low_pip_df['MAF'] = calculate_maf(low_pip_df['EAF Population'])
    
    # Handle missing VEP annotations - fill with 'Unknown'
    high_pip_df['VEP Annotation'] = high_pip_df['VEP Annotation'].fillna('Unknown')
    low_pip_df['VEP Annotation'] = low_pip_df['VEP Annotation'].fillna('Unknown')
    
    # Convert VEP annotations to numerical categories
    vep_categories = sorted(set(high_pip_df['VEP Annotation'].unique()) | 
                          set(low_pip_df['VEP Annotation'].unique()))
    vep_to_num = {vep: i for i, vep in enumerate(vep_categories)}
    
    high_pip_df['VEP_num'] = high_pip_df['VEP Annotation'].map(vep_to_num)
    low_pip_df['VEP_num'] = low_pip_df['VEP Annotation'].map(vep_to_num)
    
    # Log-transform p-values for better matching (handle p=0)
    high_pip_df['log_pval'] = -np.log10(high_pip_df['P-Value Population'].clip(lower=1e-300))
    low_pip_df['log_pval'] = -np.log10(low_pip_df['P-Value Population'].clip(lower=1e-300))
    
    # Create feature matrix for matching
    # Normalize features to similar scales
    from sklearn.preprocessing import StandardScaler
    
    high_features = high_pip_df[['MAF', 'log_pval', 'VEP_num']].values
    low_features = low_pip_df[['MAF', 'log_pval', 'VEP_num']].values
    
    scaler = StandardScaler()
    high_features_scaled = scaler.fit_transform(high_features)
    low_features_scaled = scaler.transform(low_features)
    
    # Use k-nearest neighbors to find matches
    # We want to find one low PIP variant for each high PIP variant
    n_high = len(high_pip_df)
    n_low = len(low_pip_df)
    
    if n_low < n_high:
        print(f"Warning: Only {n_low} low PIP variants available to match {n_high} high PIP variants")
        print("Returning all available low PIP variants")
        return low_pip_df
    
    # Build KNN model
    knn = NearestNeighbors(n_neighbors=min(n_neighbors, n_low), metric='euclidean')
    knn.fit(low_features_scaled)
    
    # Find nearest neighbors for each high PIP variant
    distances, indices = knn.kneighbors(high_features_scaled)
    
    # Select matched variants (without replacement)
    selected_indices = []
    available_indices = set(range(n_low))
    
    # Sort by distance to prioritize best matches
    match_pairs = []
    for i in range(n_high):
        for j, idx in enumerate(indices[i]):
            if idx in available_indices:
                match_pairs.append((distances[i][j], idx))
                available_indices.remove(idx)  # FIX: Remove immediately after selection
                break
    
    # Sort by distance and select
    match_pairs.sort()
    selected_indices = [idx for _, idx in match_pairs[:n_high]]
    
    matched_low_pip = low_pip_df.iloc[selected_indices].copy()
    
    return matched_low_pip

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter and match variants from VA Million Veteran Program fine-mapping data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i Data_S1.xlsx -o matched_variants.csv
  %(prog)s --input /path/to/Data_S1.xlsx --output /path/to/results.csv
"""
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Path to the input Excel file (Data_S1.xlsx)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='matched_variants_output.csv',
        help='Path to the output CSV file (default: matched_variants_output.csv)'
    )
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    input_file = args.input
    output_file = args.output
    
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print()
    
    # Check if input file exists
    import os
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found!")
        sys.exit(1)
    
    print("Loading data...")
    try:
        df = pd.read_excel(input_file, header=1)
    except Exception as e:
        print(f"Error loading input file: {e}")
        sys.exit(1)
    
    print(f"Total rows loaded: {len(df)}")
    print(f"Columns: {df.columns.tolist()}")
    
    # Filter for EUR population
    print("\nFiltering for EUR population...")
    df_eur = df[df['Population'] == 'EUR'].copy()
    print(f"EUR population rows: {len(df_eur)}")
    
    # Get highest PIP per variant (identified by MVP ID)
    print("\nExtracting highest PIP per variant...")
    idx_max_pip = df_eur.groupby('MVP ID')['Overall PIP'].idxmax()
    df_max_pip = df_eur.loc[idx_max_pip].copy()
    print(f"Unique variants: {len(df_max_pip)}")
    
    # Separate high PIP (>= 0.95) and low PIP (< 0.05) variants
    print("\nSeparating high and low PIP variants...")
    high_pip = df_max_pip[df_max_pip['Overall PIP'] >= 0.95].copy()
    low_pip_candidates = df_max_pip[df_max_pip['Overall PIP'] < 0.05].copy()
    
    print(f"High PIP variants (>= 0.95): {len(high_pip)}")
    print(f"Low PIP candidates (< 0.05): {len(low_pip_candidates)}")
    
    if len(high_pip) == 0:
        print("Error: No high PIP variants found!")
        sys.exit(1)
    
    if len(low_pip_candidates) == 0:
        print("Error: No low PIP variants found!")
        sys.exit(1)
    
    # Match low PIP variants to high PIP variants
    print("\nMatching low PIP variants to high PIP distribution...")
    matched_low_pip = match_distributions(high_pip, low_pip_candidates)
    print(f"Matched low PIP variants: {len(matched_low_pip)}")
    
    # Combine the two sets
    print("\nCombining high and matched low PIP variants...")
    high_pip['Set'] = 'High_PIP'
    matched_low_pip['Set'] = 'Low_PIP'
    
    combined = pd.concat([high_pip, matched_low_pip], ignore_index=True)
    
    # Select only the required columns
    output_columns = [
        'MVP ID', 'RSID', 'BP', 'BP38', 'VEP Annotation', 'CHR',
        'EAF Population', 'Beta Population', 'P-Value Population',
        'Overall PIP', 'CS-Level Pip', 'mu', 'Set'
    ]
    
    output_df = combined[output_columns].copy()
    
    # Sort by Set (High PIP first) then by Overall PIP (descending)
    output_df = output_df.sort_values(['Set', 'Overall PIP'], 
                                     ascending=[False, False])
    
    # Save to CSV
    try:
        output_df.to_csv(output_file, index=False)
        print(f"\nOutput saved to: {output_file}")
    except Exception as e:
        print(f"Error saving output file: {e}")
        sys.exit(1)
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    for set_name in ['High_PIP', 'Low_PIP']:
        set_df = output_df[output_df['Set'] == set_name]
        maf = calculate_maf(set_df['EAF Population'])
        
        print(f"\n{set_name}:")
        print(f"  N variants: {len(set_df)}")
        print(f"  PIP range: [{set_df['Overall PIP'].min():.4f}, {set_df['Overall PIP'].max():.4f}]")
        print(f"  MAF mean ± std: {maf.mean():.4f} ± {maf.std():.4f}")
        print(f"  P-value median: {set_df['P-Value Population'].median():.2e}")
        print(f"  VEP annotations: {set_df['VEP Annotation'].nunique()} unique")
        print(f"  Top VEP annotations:")
        vep_counts = set_df['VEP Annotation'].value_counts().head(5)
        for vep, count in vep_counts.items():
            print(f"    {vep}: {count}")
    
    print("\n" + "="*60)
    print("Distribution matching comparison:")
    print("="*60)
    
    high_set = output_df[output_df['Set'] == 'High_PIP']
    low_set = output_df[output_df['Set'] == 'Low_PIP']
    
    # MAF comparison
    high_maf = calculate_maf(high_set['EAF Population'])
    low_maf = calculate_maf(low_set['EAF Population'])
    
    # KS test for MAF
    ks_maf = stats.ks_2samp(high_maf, low_maf)
    print(f"MAF KS test p-value: {ks_maf.pvalue:.4f}")
    
    # KS test for p-values (log scale)
    log_pval_high = -np.log10(high_set['P-Value Population'].clip(lower=1e-300))
    log_pval_low = -np.log10(low_set['P-Value Population'].clip(lower=1e-300))
    ks_pval = stats.ks_2samp(log_pval_high, log_pval_low)
    print(f"P-value (log) KS test p-value: {ks_pval.pvalue:.4f}")
    
    # VEP annotation comparison
    high_vep = set(high_set['VEP Annotation'].unique())
    low_vep = set(low_set['VEP Annotation'].unique())
    vep_overlap = len(high_vep & low_vep) / len(high_vep | low_vep)
    print(f"VEP annotation overlap (Jaccard): {vep_overlap:.4f}")
    
    print("\nDone!")

if __name__ == '__main__':
    main()
