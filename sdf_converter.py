"""

Converts a tab delimited text file to an SDFile. The text file must have a column named smiles.  It also must have
headers for each column. This column will be used for structural information. All other columns will be stored as fields in the sdf file.
The SDFile will keep the name as text file.

example usage:

python sdf_converter.py database.txt


"""

from rdkit import Chem
import pandas as pd
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create an SDFile from a text file.')
    parser.add_argument('filename', metavar='F', type=str,
                        help='the filename of the text file')

    args = parser.parse_args()

    filename = args.filename
    basename = os.path.basename(filename).split('.')[0]

    db = pd.read_csv(filename, sep='\t')

    db.columns = list(map(lambda s: s.lower(), db.columns.astype(str)))

    if 'smiles' not in db.columns:
        raise Exception("There is no smiles column in the text file.")

    sdf_file = Chem.SDWriter('{}.sdf'.format(basename))
    mol_counter = 0

    for row, data in db.iterrows():
        smiles = data['smiles']
        mol = Chem.MolFromSmiles(smiles)

        if mol:
            for key, val in data[list(map(lambda s: s != 'smiles', data.index))].iteritems():
                mol.SetProp(key, str(val))
            sdf_file.write(mol)
            mol_counter += 1
        else:
            print("Could not create RDKit molecule for {}".format(smiles))

    sdf_file.close()

    print("Finished. Wrote {} molecules to an SDFile".format(mol_counter))

