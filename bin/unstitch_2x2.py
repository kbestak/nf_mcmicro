#!/usr/bin/env python

import argparse
import palom.pyramid
import palom.reader
import pandas as pd
import numpy as np
from loguru import logger
from argparse import ArgumentParser as AP
from os.path import abspath
import time
import os


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
        "--input",
        dest="input",
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
        "--pixel_size",
        dest="pixel_size",
        action="store",
        required=False,
        type=float,
        help="Pixel size in microns.",
    )
    inputs.add_argument(
        "--markers",
        dest="markers",
        action="store",
        required=True,
        help="Path to marker file with 'marker_name' column.",
    )
    inputs.add_argument(
        "--channels",
        dest="channels",
        action="store",
        required=False,
        nargs="+",
        help="List of indices of the channels to be extracted.",
    )
    inputs.add_argument(
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Output directory.",
    )
    inputs.add_argument("--version", action="version", version="0.0.1")
    arg = parser.parse_args()

    arg.input = abspath(arg.input)
    arg.markers = abspath(arg.markers)
    arg.output = abspath(arg.output)

    return arg


def main(args):
    input_name = os.path.basename(args.input).split('.')[0]
    reader = palom.reader.OmePyramidReader(args.input)
    pyramid = reader.pyramid[0]
    output_dir = args.output
    overlap = int(args.overlap)
    if args.pixel_size is None:
        pixel_size = 0.4977523
    else:
        pixel_size = args.pixel_size

    mid_h = pyramid.shape[1] // 2
    mid_w = pyramid.shape[2] // 2

    tile_0 = pyramid[:,0:(mid_h+overlap//2),0:(mid_w+overlap//2)]
    tile_1 = pyramid[:,0:(mid_h+overlap//2),(mid_w-overlap//2):]
    tile_2 = pyramid[:,(mid_h-overlap//2):,0:(mid_w+overlap//2)]
    tile_3 = pyramid[:,(mid_h-overlap//2):,(mid_w-overlap//2):]
    tile_list = [tile_0, tile_1, tile_2, tile_3]
    outpath_list = [os.path.join(output_dir, f'{input_name}_tile_0.ome.tif'), os.path.join(output_dir, f'{input_name}_tile_1.ome.tif'), os.path.join(output_dir, f'{input_name}_tile_2.ome.tif'), os.path.join(output_dir, f'{input_name}_tile_3.ome.tif')]
    print(outpath_list)
    if not os.path.exists(abspath(output_dir)):
        os.makedirs(abspath(output_dir))
    markers = pd.read_csv(args.markers)
    if args.channels is None:
        channels = list(range(len(markers)))
    else:
        channels = args.channels
        for i in range(len(channels)):
            channels[i] = int(channels[i])
            if channels[i] < 0 or channels[i] >= len(markers):
                raise ValueError(f"Invalid channel index: {channels[i]}")


    channel_names = markers.marker_name.values[channels]
    print(channel_names)

    for i in range(4):
        palom.pyramid.write_pyramid(
            [tile_list[i][channels]], outpath_list[i], channel_names=channel_names, downscale_factor=2, pixel_size=pixel_size, tile_size=512, compression='zlib'
        )

if __name__ == "__main__":
    # Read in arguments
    args = get_args()

    # Run script
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
