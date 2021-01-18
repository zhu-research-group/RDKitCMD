from rdkit import Chem
from typing import List
import pandas as pd
import numpy as np
import argparse, config, os

AVAILABLE_FPS = [('sag', 'Saagar')]

def calc_saagar(
                mols: List[Chem.Mol],
                saagar_filename=config.SAAGAR_SMARTS,
                counts=True,
                to_file=None
                ):

    """
    returns a pandas dataframe of Saagar finger prints
    mols: a list of Chem.Mol
    saagar_filename: path of Saagar fingerprints.
    counts: if True return number of matches for give mol, else return a 1 if any hits at all
    to_file: if None returns a data frame, else expects a path and file name to write to
    """

    if not os.path.exists(saagar_filename):
        raise Exception('File {} not found'.format(saagar_filename))
    else:
        saagar_df = pd.read_csv(saagar_filename, sep='\t')

    fps = []
    codes = saagar_df.Saagar_Label.values.tolist()

    mol_list = [mol if mol else None for mol in mols]
    num_fps = saagar_df.shape[0]
    for i, mol in enumerate(mol_list):
        print("Processing {} molecule # {} of ".format(i+1, len(mol_list)))
        fp = []
        if mol:
            for row, saagar_data in saagar_df.iterrows():
                label = saagar_data.Saagar_Label
                smarts = saagar_data.SMARTS
                smarts_mol = Chem.MolFromSmarts(smarts)
                matches = len(mol.GetSubstructMatches(smarts_mol))
                if not counts:
                    matches = int(matches > 0)
                    fp.append(matches)
                else:
                    fp.append(matches)

        else:
            fp.append([np.nan]*num_fps)
        fps.append(fp)

    fp_frame = pd.DataFrame(fps, columns=codes)

    if to_file == None:
        return fp_frame
    else:
        fp_frame.to_csv(to_file)






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process RDKit descriptors from smiles.txt file.')
    parser.add_argument('--in_file', metavar='inf', type=str,
                        help='the filename of an sdf file')

    parser.add_argument('--out_file', metavar='outf', type=str,
                        help='the filename of an sdf file')

    parser.add_argument('--fp', metavar='fp', type=str,
                        help='the filename of an sdf file')


    args = parser.parse_args()

    if args.fp not in [fp[0] for fp in AVAILABLE_FPS]:
        e_s = "{} not a valid fingerpint choice\nAvailable fingerprints are:\n".format(args.fp)
        for arg, name in AVAILABLE_FPS:
            e_s = e_s + "{}: {}\n".format(arg, name)
        raise Exception(e_s)

    mols = Chem.SDMolSupplier(args.in_file)

    if args.fp == 'sag':
        calc_saagar(mols, to_file=args.out_file)






