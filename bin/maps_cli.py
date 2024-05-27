#!/usr/bin/env python

from argparse import ArgumentParser as AP
import argparse
from maps.cell_phenotyping import Predictor
from os.path import abspath
import time
import pandas as pd

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
        "--input_csv",
        dest="input_csv",
        action="store",
        required=True,
        help="File path to input quantification table.",
    )
    inputs.add_argument(
        "--pretrained_model",
        dest="pretrained_model",
        action="store",
        required=True,
        help="File path to best checkpoint (pretrained model)",
    )
    inputs.add_argument(
        "--labelID",
        dest="labelID",
        action="store",
        required=True,
        help="File path to csv containing mapping between indices (from pretrained model) and labels (phenotypes)",
    )
    inputs.add_argument(
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Output file (csv) path.",
    )
    inputs.add_argument("--version", action="version", version="0.0.1")
    arg = parser.parse_args()

    arg.input_csv = abspath(arg.input_csv)
    arg.pretrained_model = abspath(arg.pretrained_model)
    arg.labelID = abspath(arg.labelID)
    arg.output = abspath(arg.output)

    return arg

def main(args):
    pretrained_model_checkpoint_path = args.pretrained_model

    labelID = pd.read_csv(args.labelID)
    num_classes = len(labelID)
    prep_table = pd.read_csv(args.input_csv)
    prep_table.rename(columns={'Area': 'cellSize'}, inplace=True)
    prep_table.drop(columns=['CellID', 'DAPI', 'AF', 'X_centroid', 'Y_centroid', 'MajorAxisLength', 'MinorAxisLength', 'Extent', 'Orientation'], inplace=True)
    prep_table.to_csv(f'{args.output}_prep_table.csv', index=False)

    data_path = f'{args.output}_prep_table.csv'
    model = Predictor(model_checkpoint_path=pretrained_model_checkpoint_path,
        num_features=8,
        num_classes=6,
        batch_size=128
        )

    pred_labels, pred_probs = model.predict(data_path)

    prep_table = pd.read_csv(args.input_csv)
    prep_table["labelID"] = pred_labels
    #prep_table["prediction_prob"] = pred_probs

    labelID_mapping = pd.Series(labelID.cell_label.values, index=labelID.labelID).to_dict()
    prep_table['phenotype_MAPS'] = prep_table['labelID'].map(labelID_mapping)

    prep_table.to_csv(args.output, index=False)


if __name__ == "__main__":
    # Read in arguments
    args = get_args()

    # Run script
    st = time.time()
    main(args)
    rt = time.time() - st

    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
