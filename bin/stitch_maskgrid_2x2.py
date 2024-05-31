#!/usr/bin/env python

import tifffile
import numpy as np
import skimage
import matplotlib.pyplot as plt
from os.path import abspath
from argparse import ArgumentParser as AP
import argparse
import time
import copy
import skimage.measure

def longest_contiguous_subsequence(indices):
    if len(indices) == 0:
        return np.array([])

    longest_start = 0
    longest_length = 1
    current_start = 0
    current_length = 1

    for i in range(1, len(indices)):
        if indices[i] == indices[i-1] + 1:
            current_length += 1
        else:
            if current_length > longest_length:
                longest_length = current_length
                longest_start = current_start
            current_start = i
            current_length = 1

    if current_length > longest_length:
        longest_length = current_length
        longest_start = current_start

    return indices[longest_start:longest_start + longest_length]

def stitch_tiles_vertical(tile1, tile2, overlap):
    overlap1 = tile1[-overlap:, :]
    overlap2 = tile2[:overlap, :]

    keep_overlap = np.logical_and(overlap1 == 0, overlap2 == 0)

    image1_no_overlap = tile1[:-overlap, :]
    image2_no_overlap = tile2[overlap:, :]

    label1_no_overlap = skimage.measure.label(image1_no_overlap, background=0)
    maxlabel = label1_no_overlap.max()

    label2_no_overlap = skimage.measure.label(image2_no_overlap, background=0)
    label2_no_overlap[label2_no_overlap != 0] += maxlabel
    maxlabel = label2_no_overlap.max()

    label1_overlap = skimage.measure.label(overlap1, background=0)
    label1_overlap[label1_overlap != 0] += maxlabel
    maxlabel = label1_overlap.max()

    label2_overlap = skimage.measure.label(overlap2, background=0)
    label2_overlap[label2_overlap != 0] += maxlabel

    final_overlap = np.zeros((overlap, tile1.shape[1]), dtype=np.uint32)

    regions1_overlap = skimage.measure.regionprops(label1_overlap)
    centroids1_overlap = np.array([[int(np.floor(r.centroid)[0]), int(np.floor(r.centroid)[1])] for r in regions1_overlap])
    if centroids1_overlap.size == 0:
        centroids1_overlap = np.empty((0, 2), dtype=int)

    centroids1_overlap_above = centroids1_overlap[centroids1_overlap[:, 0] < overlap // 2]
    for centroid in centroids1_overlap_above:
        if label1_overlap[centroid[0], centroid[1]] != 0:
            final_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = label1_overlap[centroid[0], centroid[1]]
            label1_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = 0
            label2_overlap[label2_overlap == label2_overlap[centroid[0], centroid[1]]] = 0

    label2_overlap = skimage.measure.label(label2_overlap, background=0)
    label2_overlap[label2_overlap != 0] += maxlabel

    regions2_overlap = skimage.measure.regionprops(label2_overlap)
    centroids2_overlap = np.array([[int(np.floor(r.centroid)[0]), int(np.floor(r.centroid)[1])] for r in regions2_overlap])
    if centroids2_overlap.size == 0:
        centroids2_overlap = np.empty((0, 2), dtype=int)

    for centroid in centroids2_overlap:
        if label2_overlap[centroid[0], centroid[1]] != 0:
            final_overlap[label2_overlap == label2_overlap[centroid[0], centroid[1]]] = label2_overlap[centroid[0], centroid[1]]
            label1_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = 0
            label2_overlap[label2_overlap == label2_overlap[centroid[0], centroid[1]]] = 0

    regions1_leftover = skimage.measure.regionprops(label1_overlap)
    centroids1_leftover = np.array([[int(np.floor(r.centroid)[0]), int(np.floor(r.centroid)[1])] for r in regions1_leftover])
    if centroids1_leftover.size == 0:
        centroids1_leftover = np.empty((0, 2), dtype=int)

    for centroid in centroids1_leftover:
        if label1_overlap[centroid[0], centroid[1]] != 0:
            final_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = label1_overlap[centroid[0], centroid[1]]

    final_overlap = skimage.measure.label(final_overlap, background=0)

    upper_edge_unique = np.unique(final_overlap[0, :])
    upper_edge_unique = upper_edge_unique[upper_edge_unique != 0]

    if len(upper_edge_unique) > 0:
        for unique_edge in upper_edge_unique:
            matching_indices = np.where(final_overlap[0, :] == unique_edge)[0]
            matching_indices = longest_contiguous_subsequence(matching_indices)
            if len(matching_indices) > 0:
                position_to_check = int(np.floor(np.mean(matching_indices)))
                if image1_no_overlap[-1, position_to_check] != 0:
                    final_overlap[~keep_overlap] = np.where(final_overlap[~keep_overlap] == final_overlap[0, position_to_check], image1_no_overlap[-1, position_to_check], final_overlap[~keep_overlap])

    image_combined = np.concatenate((image1_no_overlap, final_overlap), axis=0)
    lower_edge_unique = np.unique(image_combined[-1, :])
    lower_edge_unique = lower_edge_unique[lower_edge_unique != 0]

    if len(lower_edge_unique) > 0:
        for unique_edge in lower_edge_unique:
            matching_indices = np.where(image_combined[-1, :] == unique_edge)[0]
            matching_indices = longest_contiguous_subsequence(matching_indices)
            if len(matching_indices) > 0:
                position_to_check = int(np.floor(np.mean(matching_indices)))
                if image2_no_overlap[0, position_to_check] != 0:
                    image_combined[image_combined == image_combined[-1, position_to_check]] = image2_no_overlap[0, position_to_check]

    image_combined = np.concatenate((image_combined, image2_no_overlap), axis=0)
    image_combined = skimage.measure.label(image_combined, background=0)

    return image_combined

def stitch_tiles_horizontal(tile1, tile2, overlap):
    overlap1 = tile1[:, -overlap:]
    overlap2 = tile2[:, :overlap]

    keep_overlap = np.logical_and(overlap1 == 0, overlap2 == 0)

    image1_no_overlap = tile1[:, :-overlap]
    image2_no_overlap = tile2[:, overlap:]

    label1_no_overlap = skimage.measure.label(image1_no_overlap, background=0)
    maxlabel = label1_no_overlap.max()

    label2_no_overlap = skimage.measure.label(image2_no_overlap, background=0)
    label2_no_overlap[label2_no_overlap != 0] += maxlabel
    maxlabel = label2_no_overlap.max()

    label1_overlap = skimage.measure.label(overlap1, background=0)
    label1_overlap[label1_overlap != 0] += maxlabel
    maxlabel = label1_overlap.max()

    label2_overlap = skimage.measure.label(overlap2, background=0)
    label2_overlap[label2_overlap != 0] += maxlabel

    final_overlap = np.zeros((tile1.shape[0], overlap), dtype=np.uint32)

    regions1_overlap = skimage.measure.regionprops(label1_overlap)
    centroids1_overlap = np.array([[int(np.floor(r.centroid)[0]), int(np.floor(r.centroid)[1])] for r in regions1_overlap])
    if centroids1_overlap.size == 0:
        centroids1_overlap = np.empty((0, 2), dtype=int)

    centroids1_overlap_left = centroids1_overlap[centroids1_overlap[:, 1] < overlap // 2]
    for centroid in centroids1_overlap_left:
        if label1_overlap[centroid[0], centroid[1]] != 0:
            final_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = label1_overlap[centroid[0], centroid[1]]
            label1_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = 0
            label2_overlap[label2_overlap == label2_overlap[centroid[0], centroid[1]]] = 0

    label2_overlap = skimage.measure.label(label2_overlap, background=0)
    label2_overlap[label2_overlap != 0] += maxlabel

    regions2_overlap = skimage.measure.regionprops(label2_overlap)
    centroids2_overlap = np.array([[int(np.floor(r.centroid)[0]), int(np.floor(r.centroid)[1])] for r in regions2_overlap])
    if centroids2_overlap.size == 0:
        centroids2_overlap = np.empty((0, 2), dtype=int)

    for centroid in centroids2_overlap:
        if label2_overlap[centroid[0], centroid[1]] != 0:
            final_overlap[label2_overlap == label2_overlap[centroid[0], centroid[1]]] = label2_overlap[centroid[0], centroid[1]]
            label1_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = 0
            label2_overlap[label2_overlap == label2_overlap[centroid[0], centroid[1]]] = 0

    regions1_leftover = skimage.measure.regionprops(label1_overlap)
    centroids1_leftover = np.array([[int(np.floor(r.centroid)[0]), int(np.floor(r.centroid)[1])] for r in regions1_leftover])
    if centroids1_leftover.size == 0:
        centroids1_leftover = np.empty((0, 2), dtype=int)

    for centroid in centroids1_leftover:
        if label1_overlap[centroid[0], centroid[1]] != 0:
            final_overlap[label1_overlap == label1_overlap[centroid[0], centroid[1]]] = label1_overlap[centroid[0], centroid[1]]

    final_overlap = skimage.measure.label(final_overlap, background=0)

    left_edge_unique = np.unique(final_overlap[:, 0])
    left_edge_unique = left_edge_unique[left_edge_unique != 0]

    if len(left_edge_unique) > 0:
        for unique_edge in left_edge_unique:
            matching_indices = np.where(final_overlap[:, 0] == unique_edge)[0]
            matching_indices = longest_contiguous_subsequence(matching_indices)
            if len(matching_indices) > 0:
                position_to_check = int(np.floor(np.mean(matching_indices)))
                if image1_no_overlap[position_to_check, -1] != 0:
                    final_overlap[~keep_overlap] = np.where(final_overlap[~keep_overlap] == final_overlap[position_to_check, 0], image1_no_overlap[position_to_check, -1], final_overlap[~keep_overlap])

    image_combined = np.concatenate((image1_no_overlap, final_overlap), axis=1)
    right_edge_unique = np.unique(image_combined[:, -1])
    right_edge_unique = right_edge_unique[right_edge_unique != 0]

    if len(right_edge_unique) > 0:
        for unique_edge in right_edge_unique:
            matching_indices = np.where(image_combined[:, -1] == unique_edge)[0]
            matching_indices = longest_contiguous_subsequence(matching_indices)
            if len(matching_indices) > 0:
                position_to_check = int(np.floor(np.mean(matching_indices)))
                if image2_no_overlap[position_to_check, 0] != 0:
                    image_combined[image_combined == image_combined[position_to_check, -1]] = image2_no_overlap[position_to_check, 0]

    image_combined = np.concatenate((image_combined, image2_no_overlap), axis=1)
    image_combined = skimage.measure.label(image_combined, background=0)

    return image_combined



def get_args():
    # Script description
    description = """Stitch segmentation masks from 2x2 grid"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Sections
    inputs = parser.add_argument_group(
        title="Required Input", description="Path to required input files"
    )
    inputs.add_argument(
        "--tile_0",
        dest="tile_0",
        action="store",
        required=True,
        help="File path to input image.",
    )
    inputs.add_argument(
        "--tile_1",
        dest="tile_1",
        action="store",
        required=True,
        help="File path to input image.",
    )
    inputs.add_argument(
        "--tile_2",
        dest="tile_2",
        action="store",
        required=True,
        help="File path to input image.",
    )
    inputs.add_argument(
        "--tile_3",
        dest="tile_3",
        action="store",
        required=True,
        help="File path to input image.",
    )
    inputs.add_argument(
        "--overlap",
        dest="overlap",
        action="store",
        required=True,
        help="Overlap between tiles.",
    )
    inputs.add_argument(
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Output file path.",
    )
    inputs.add_argument("--version", action="version", version="0.0.1")
    arg = parser.parse_args()

    arg.tile_0 = abspath(arg.tile_0)
    arg.tile_1 = abspath(arg.tile_1)
    arg.tile_2 = abspath(arg.tile_2)
    arg.tile_3 = abspath(arg.tile_3)
    arg.output = abspath(arg.output)

    return arg



def main(args):
    tile_0 = tifffile.imread(args.tile_0)
    tile_1 = tifffile.imread(args.tile_1)
    tile_2 = tifffile.imread(args.tile_2)
    tile_3 = tifffile.imread(args.tile_3)
    overlap = int(args.overlap)

    left = stitch_tiles_vertical(tile_0, tile_2, overlap)
    tile_0 = None
    tile_2 = None
    right = stitch_tiles_vertical(tile_1, tile_3, overlap)
    tile_1 = None
    tile_3 = None
    stitched_image = stitch_tiles_horizontal(left, right, overlap)

    image_combined = skimage.morphology.remove_small_objects(stitched_image, min_size=10, connectivity=0)

    stitched_image = skimage.measure.label(stitched_image, background=0)

    tifffile.imwrite(args.output, stitched_image, compression='zlib')

if __name__ == "__main__":
    # Read in arguments
    args = get_args()

    # Run script
    st = time.time()
    main(args)
    rt = time.time() - st

    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
