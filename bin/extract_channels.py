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
import tifffile


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
    # add a true-false argument called "rescale"
    inputs.add_argument(
        "--rescale",
        dest="rescale",
        action="store_true",
        required=False,
        help="Rescale the channel by ",
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
    inputs.add_argument(
        "--output_tag",
        dest="output_tag",
        action="store",
        required=True,
        help="Tag used for naming output files - 'Tag_ChannelName.ome.tif'",
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
    if args.pixel_size is None:
        pixel_size = 0.4977523
    else:
        pixel_size = args.pixel_size

    markers = pd.read_csv(args.markers)

    if args.channels is None:
        channels = list(range(len(markers)))
    else:
        channels = args.channels
        for i in range(len(channels)):
            channels[i] = int(channels[i])
            if channels[i] < 0 or channels[i] >= len(markers):
                raise ValueError(f"Invalid channel index: {channels[i]}")

    print(pyramid[channels[0]].max().compute())

    channel_names = markers.marker_name.values[channels]
    outpath_list = [os.path.join(output_dir, f'{args.output_tag}_{channel_names[i]}.tif') for i in range(len(channels))]
    print(f'Pixel size is: {pixel_size}')
    print(outpath_list)
    for idx, channel in enumerate(channels):
        img_out = pyramid[int(channel)]
        out_path = outpath_list[idx]
        tifffile.imwrite(out_path, img_out.compute(), imagej=False, resolution=(10000 / pixel_size, 10000 / pixel_size, "centimeter"))

if __name__ == "__main__":
    # Read in arguments
    args = get_args()

    # Run script
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
