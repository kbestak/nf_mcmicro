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

def stitch_tiles_vertical(tile1, tile2, overlap):
    image1 = copy.copy(tile1)
    image2 = copy.copy(tile2)
    overlap1 = image1[-overlap:, :]
    overlap2 = image2[:overlap, :]

    image1_no_overlap = image1[:-overlap, :]
    image2_no_overlap = image2[overlap:, :]

    label1_no_overlap = skimage.measure.label(image1_no_overlap, background=0)
    maxlabel = label1_no_overlap.max()

    label2_no_overlap = skimage.measure.label(image2_no_overlap, background=0)
    label2_no_overlap = np.where(label2_no_overlap == 0, 0, label2_no_overlap+maxlabel)
    maxlabel = label2_no_overlap.max()

    label1_overlap = skimage.measure.label(overlap1, background=0)
    label1_overlap = np.where(label1_overlap == 0, 0, label1_overlap+maxlabel)
    maxlabel = label1_overlap.max()

    label2_overlap = skimage.measure.label(overlap2, background=0)
    label2_overlap = np.where(label2_overlap == 0, 0, label2_overlap+maxlabel)
    maxlabel = label2_overlap.max()

    regions1_overlap = skimage.measure.regionprops(label1_overlap)
    centroids1_overlap = np.array([np.floor(r.centroid).astype(np.uint32) for r in regions1_overlap])
    final_overlap = np.zeros((overlap, image1.shape[1]), dtype=np.uint32)

    centroids1_overlap_above = centroids1_overlap[centroids1_overlap[:,0] < overlap//2,:]
    if len(centroids1_overlap_above) > 0:
        for centroid in centroids1_overlap_above:
            if label1_overlap[centroid[0], centroid[1]] != 0:
                final_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], label1_overlap, final_overlap)
                label1_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], 0, label1_overlap)
                label2_overlap = np.where(label2_overlap == label2_overlap[centroid[0], centroid[1]], 0, label2_overlap)
    label2_overlap = skimage.measure.label(label2_overlap, background=0)
    label2_overlap = np.where(label2_overlap == 0, 0, label2_overlap+maxlabel)

    regions2_overlap = skimage.measure.regionprops(label2_overlap)
    centroids2_overlap = np.array([np.round(r.centroid).astype(np.uint32) for r in regions2_overlap])

    if len(centroids2_overlap) > 0:
        for centroid in centroids2_overlap:
            if label2_overlap[centroid[0], centroid[1]] != 0:
                final_overlap = np.where(label2_overlap == label2_overlap[centroid[0], centroid[1]], label2_overlap, final_overlap)
                label1_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], 0, label1_overlap)
                label2_overlap = np.where(label2_overlap == label2_overlap[centroid[0], centroid[1]], 0, label2_overlap)

    regions1_leftover = skimage.measure.regionprops(label1_overlap)
    centroids1_leftover = np.array([np.round(r.centroid).astype(np.uint32) for r in regions1_leftover])

    if len(centroids1_leftover) > 0:
        for centroid in centroids1_leftover:
            if label1_overlap[centroid[0], centroid[1]] != 0:
                final_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], label1_overlap, final_overlap)

    final_overlap = skimage.measure.label(final_overlap, background=0)

    upper_edge_unique = np.unique(final_overlap[0,:])
    upper_edge_unique = upper_edge_unique[upper_edge_unique != 0]
    lower_edge_unique = np.unique(final_overlap[-1,:])
    lower_edge_unique = lower_edge_unique[lower_edge_unique != 0]

    if len(upper_edge_unique) > 0:
        for unique_edge in upper_edge_unique:
            position_to_check = np.rint(np.median(np.where(final_overlap[0,:] == unique_edge))).astype(np.uint32)
            if image1_no_overlap[-1,position_to_check] != 0:
                print('position_to_check:', position_to_check)
                print('unique_edge:', unique_edge)
                print('value', image1_no_overlap[-1,position_to_check])
                final_overlap = np.where(final_overlap == final_overlap[0, position_to_check], image1_no_overlap[-1,position_to_check], final_overlap)

    if len(lower_edge_unique) > 0:
        for unique_edge in lower_edge_unique:
            position_to_check = np.rint(np.median(np.where(final_overlap[-1,:] == unique_edge))).astype(np.uint32)
            if image2_no_overlap[0,position_to_check] != 0:
                print('position_to_check:', position_to_check)
                print('unique_edge:', unique_edge)
                print('value', image2_no_overlap[0,position_to_check])
                final_overlap = np.where(final_overlap == final_overlap[-1, position_to_check], image2_no_overlap[0,position_to_check], final_overlap)

    image_combined = np.concatenate((image1_no_overlap, final_overlap, image2_no_overlap), axis=0)
    image_combined = skimage.measure.label(image_combined, background=0)

    return image_combined


def stitch_tiles_horizontal(tile1, tile2, overlap):
    image1 = copy.copy(tile1)
    image2 = copy.copy(tile2)
    overlap1 = image1[:,-overlap:]
    overlap2 = image2[:,:overlap]

    image1_no_overlap = image1[:,:-overlap]
    image2_no_overlap = image2[:,overlap:]

    label1_no_overlap = skimage.measure.label(image1_no_overlap, background=0)
    maxlabel = label1_no_overlap.max()

    label2_no_overlap = skimage.measure.label(image2_no_overlap, background=0)
    label2_no_overlap = np.where(label2_no_overlap == 0, 0, label2_no_overlap+maxlabel)
    maxlabel = label2_no_overlap.max()

    label1_overlap = skimage.measure.label(overlap1, background=0)
    label1_overlap = np.where(label1_overlap == 0, 0, label1_overlap+maxlabel)
    maxlabel = label1_overlap.max()

    label2_overlap = skimage.measure.label(overlap2, background=0)
    label2_overlap = np.where(label2_overlap == 0, 0, label2_overlap+maxlabel)
    maxlabel = label2_overlap.max()

    regions1_overlap = skimage.measure.regionprops(label1_overlap)

    centroids1_overlap = np.array([np.floor(r.centroid).astype(np.uint32) for r in regions1_overlap])

    final_overlap = np.zeros((image1.shape[0],overlap), dtype=np.uint32)


    centroids1_overlap_above = centroids1_overlap[:,centroids1_overlap[0,:] < overlap//2]
    if len(centroids1_overlap_above) > 0:
        for centroid in centroids1_overlap_above:
            if label1_overlap[centroid[0], centroid[1]] != 0:
                final_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], label1_overlap, final_overlap)
                label1_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], 0, label1_overlap)
                label2_overlap = np.where(label2_overlap == label2_overlap[centroid[0], centroid[1]], 0, label2_overlap)
    label2_overlap = skimage.measure.label(label2_overlap, background=0)
    label2_overlap = np.where(label2_overlap == 0, 0, label2_overlap+maxlabel)

    regions2_overlap = skimage.measure.regionprops(label2_overlap)
    centroids2_overlap = np.array([np.round(r.centroid).astype(np.uint32) for r in regions2_overlap])

    if len(centroids2_overlap) > 0:
        for centroid in centroids2_overlap:
            if label2_overlap[centroid[0], centroid[1]] != 0:
                final_overlap = np.where(label2_overlap == label2_overlap[centroid[0], centroid[1]], label2_overlap, final_overlap)
                label1_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], 0, label1_overlap)
                label2_overlap = np.where(label2_overlap == label2_overlap[centroid[0], centroid[1]], 0, label2_overlap)

    regions1_leftover = skimage.measure.regionprops(label1_overlap)
    centroids1_leftover = np.array([np.round(r.centroid).astype(np.uint32) for r in regions1_leftover])

    if len(centroids1_leftover) > 0:
        for centroid in centroids1_leftover:
            if label1_overlap[centroid[0], centroid[1]] != 0:
                final_overlap = np.where(label1_overlap == label1_overlap[centroid[0], centroid[1]], label1_overlap, final_overlap)

    final_overlap = skimage.measure.label(final_overlap, background=0)

    upper_edge_unique = np.unique(final_overlap[:,0])
    upper_edge_unique = upper_edge_unique[upper_edge_unique != 0]
    lower_edge_unique = np.unique(final_overlap[:,-1])
    lower_edge_unique = lower_edge_unique[lower_edge_unique != 0]

    if len(upper_edge_unique) > 0:
        for unique_edge in upper_edge_unique:
            position_to_check = np.rint(np.median(np.where(final_overlap[:,0] == unique_edge))).astype(np.uint32)
            if image1_no_overlap[position_to_check,-1] != 0:
                ##print('position_to_check:', position_to_check)
                ##print('unique_edge:', unique_edge)
                ##print('value', image1_no_overlap[position_to_check,-1])
                final_overlap = np.where(final_overlap == final_overlap[position_to_check,0], image1_no_overlap[position_to_check,-1], final_overlap)

    if len(lower_edge_unique) > 0:
        for unique_edge in lower_edge_unique:
            position_to_check = np.rint(np.median(np.where(final_overlap[:,-1] == unique_edge))).astype(np.uint32)
            if image2_no_overlap[position_to_check,0] != 0:
                ##print('position_to_check:', position_to_check)
                ##print('unique_edge:', unique_edge)
                ##print('value', image2_no_overlap[position_to_check,0])
                final_overlap = np.where(final_overlap == final_overlap[position_to_check,-1], image2_no_overlap[position_to_check,0], final_overlap)

    final_overlap = skimage.morphology.remove_small_objects(final_overlap, min_size=10, connectivity=0)

    image_combined = np.concatenate((image1_no_overlap, final_overlap, image2_no_overlap), axis=1)
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
    right = stitch_tiles_vertical(tile_1, tile_3, overlap)
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
