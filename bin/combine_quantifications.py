#!/usr/bin/env python
import os
import argparse
import pathlib
from pathlib import Path
import pandas as pd

if os.name == 'nt':
    pathlib.PosixPath = pathlib.WindowsPath

def get_args():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Combine spot2cell and mcquant outputs based on CellID column.')
    parser.add_argument('-s', '--spot2cell', type=str, required=True, help='Path to the spot2cell output file with spot counts.')
    parser.add_argument('-m', '--mcquant', type=str, required=True, help='Path to the quantification file with mean intensities.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output csv file.')
    args = parser.parse_args()

    # Standardize paths
    # Convert WindowsPath to PosixPath
    args.spot2cell = Path(args.spot2cell).resolve()
    args.mcquant = Path(args.mcquant).resolve()
    args.output = Path(args.output).resolve()

    return args

def combine_quantifications(spot2cell, mcquant, output_path):
    try:
        quantification = pd.read_csv(mcquant)
        spots = pd.read_csv(spot2cell)

        # Debug: Print DataFrame heads to check contents
        print("Quantification DataFrame head:\n", quantification.head())
        print("Spots DataFrame head:\n", spots.head())

        combined = pd.merge(quantification, spots, on='CellID')
        combined.rename(columns={'spot_count': 'cGAS_spot_count'}, inplace=True)

        # Find the index of 'X_centroid' and insert 'cGAS_spot_count' before it
        x_centroid_index = combined.columns.get_loc('X_centroid')
        cols = combined.columns.tolist()
        cols.insert(x_centroid_index, cols.pop(cols.index('cGAS_spot_count')))
        combined = combined[cols]

        # Debug: Print combined DataFrame head to check contents
        print("Combined DataFrame head with cGAS_spot_count inserted:\n", combined.head())

        combined.to_csv(output_path, index=False)
        print(f"DataFrame saved to {output_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    # Get the arguments
    args = get_args()

    combine_quantifications(args.spot2cell, args.mcquant, args.output)

if __name__ == '__main__':
    main()
