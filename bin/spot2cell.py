#!/usr/bin/env python
import os
import argparse
import pathlib
from pathlib import Path
import numpy as np
from tifffile import imread
if os.name == 'nt':
    pathlib.PosixPath = pathlib.WindowsPath

def get_args():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Assign spots to cells.')
    parser.add_argument('-s', '--spots', type=str, help='Path to the spot table csv file.')
    parser.add_argument('-m', '--mask', type=str, help='Path to the cell mask tif file.')
    parser.add_argument('-o', '--output', type=str, help='Path to the output csv file.')
    args = parser.parse_args()

    # Standardize paths
    # Convert WindowsPath to PosixPath
    args.spots = Path(args.spots).resolve()
    args.mask = Path(args.mask).resolve()
    args.output = Path(args.output).resolve()

    return args


def Spot2Cell(spots_path, mask_path, output_path):
    """
    Assign spots to cells.
    :param spots_path: path to the spot table csv file
    :param mask_path: Path to the cell mask tif file
    :param output_path: Path to the output csv file
    :return: None
    """
    # Read spot table using numpy
    spot_table = np.genfromtxt(spots_path, delimiter=',', skip_header=1, dtype=np.uint32)
    print(f"Total number of spots             = {spot_table.shape[0]}")

    # Read cell mask using tifffile
    cell_mask = imread(mask_path)
    print(f"Total number of cells             = {cell_mask.max()}")

    # Index mask file at the spot table coordinate to get the cell id of each spot (0 for background)
    spots = cell_mask[spot_table[:, 0], spot_table[:, 1]]

    # Get the unique cell ids and their counts
    cell_ids, counts = np.unique(spots, return_counts=True)

    # Get the number of background spots
    background_spots = counts[0]
    print(f"Total number of cells with spots  = {cell_ids.shape[0]-1}")
    print(f"Total number of background spots  = {background_spots}")

    # Get the non-zero cell ids and their counts
    cell_ids = cell_ids[1:]
    counts = counts[1:]

    # Create a zero numpy array the size of the number of cell ids in the mask
    cell_spots = np.zeros((cell_mask.max(), 2), dtype=np.uint32)

    # Fill the first column of cell_counts with consecutive cell ids
    cell_spots[:, 0] = np.arange(1, cell_mask.max() + 1)

    # Assign the counts to the corresponding cell ids
    cell_spots[cell_ids - 1, 1] = counts

    # Save the cell spots to a csv file with a header row
    np.savetxt(output_path, cell_spots, delimiter=',', header='CellID,spot_count', comments='', fmt='%d')


def main():
    # Get the arguments
    args = get_args()

    # Create an instance of the Spot2Cell class.
    Spot2Cell(args.spots, args.mask, args.output)


if __name__ == '__main__':
    main()
