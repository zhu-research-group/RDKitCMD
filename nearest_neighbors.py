"""

Takes a source database (SDFile) and finds the nearest neighbor in a target database (SDFile), based of MACCS
keys and the Tanimoto coefficient.  Prints out a matrix of Tanimoto coefficients and a file listing nearest neighbors
with score. SDFile must contain a name filed and is specified after the database name in the command line.

example usage:

python nearest_neighbors.py source.sdf mol_id target.sdf name

"""
from rdkit.Chem import SDMolSupplier, Mol, MACCSkeys
from typing import List
import pandas as pd
from pandas import DataFrame
import argparse
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import jaccard


import os


def create_MACCS_keys(mol: Mol) -> List:
    """ creates MACC keys for a single RDKit Mol """
    return [int(x) for x in MACCSkeys.GenMACCSKeys(mol)]

def create_MACCS_keys_from_mols(mols: SDMolSupplier, name_col: str) -> DataFrame:
    """ creates MACCS keys for each mol in a SDMolSupplier object and index them by the name column """

    names = []
    data = []

    for mol in mols:
        data.append(create_MACCS_keys(mol))
        names.append(mol.GetProp(name_col))

    return DataFrame(data, index=names)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create an SDFile from a text file.')
    parser.add_argument('source_db', metavar='sdb', type=str,
                        help='the SDFile of the source file')
    parser.add_argument('source_name', metavar='sn', type=str,
                        help='the field name of the names in the source db')
    parser.add_argument('target_db', metavar='tdb', type=str,
                        help='the SDFile of the target file')
    parser.add_argument('target_name', metavar='tn', type=str,
                        help='the field name of the names in the target db')

    args = parser.parse_args()

    source_mols = SDMolSupplier(args.source_db)
    target_mols = SDMolSupplier(args.target_db)

    source_name_field = args.source_name
    target_name_field = args.target_name

    source_fps = create_MACCS_keys_from_mols(source_mols, source_name_field)
    target_fps = create_MACCS_keys_from_mols(target_mols, target_name_field)

    similarity_matrix = 1-pairwise_distances(source_fps, target_fps, metric='jaccard')

    nns, similarities = similarity_matrix.argmax(1), similarity_matrix.max(1)

    neighbor_file = open('neighbors.txt', 'w')
    neighbor_file.write('Cmp' + '\t' + 'NN' + '\t' + 'Similarity' + '\n')

    for cmp, similarity, neighbor in zip(source_fps.index, similarities, nns):
        neighbor_file.write(str(cmp) + '\t' + str(target_fps.index[neighbor]) + '\t' + str(similarity) + '\n')

    similarity_matrix = pd.DataFrame(similarity_matrix, index=source_fps.index, columns=target_fps.index)

    neighbor_file.close()





