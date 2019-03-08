from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
import pandas as pd
import argparse




def calc_descriptors_from_file(filename):
    """ calculates rdkit descriptors from a smiles.txt file """
    df = pd.read_table(filename, header=None)
    df['rdkit'] = [Chem.MolFromSmiles(smi) for smi in df[0]]
    df[pd.isnull(df['rdkit'])].to_csv('errors.txt')
    df = df[~pd.isnull(df['rdkit'])]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc[0] for desc in Descriptors.descList])

    X = pd.DataFrame([list(calc.CalcDescriptors(mol)) for mol in df['rdkit']],
                     columns=list(calc.GetDescriptorNames()),
                     index=df[0].values)

    X.to_csv('descriptors.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process RDKit descriptors from smiles.txt file.')
    parser.add_argument('filename', metavar='F', type=str,
                        help='the filename of the smiles.txt file')

    args = parser.parse_args()

    calc_descriptors_from_file(args.filename)