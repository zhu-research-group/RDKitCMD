from rdkit import Chem
from rdkit.Chem import Descriptors


smile_list = ['c1ccccc1C(=O)O','cc']


def descrptor_1smiles(smiles_str):
    m = Chem.MolFromSmiles(smiles_str)
    descriptor = []
    descriptor.append(Descriptors.HeavyAtomMolWt(m))
    descriptor.append(Descriptors.BalabanJ(m))
    descriptor.append(Descriptors.BertzCT(m))
    descriptor.append(Descriptors.Ipc(m))
    descriptor.append(Descriptors.HallKierAlpha(m))
    descriptor.append(Descriptors.Kappa1(m))
    descriptor.append(Descriptors.Kappa3(m))
    descriptor.append(Descriptors.Chi0(m))
    descriptor.append(Descriptors.Chi1(m))
    descriptor.append(Descriptors.Chi0n(m))
    descriptor.append(Descriptors.Chi4n(m))
    descriptor.append(Descriptors.Chi0v(m))
    descriptor.append(Descriptors.Chi4v(m))
    descriptor.append(Descriptors.MolLogP(m))
    descriptor.append(Descriptors.MolMR(m))
    descriptor.append(Descriptors.MolWt(m))
    descriptor.append(Descriptors.ExactMolWt(m))
    descriptor.append(Descriptors.HeavyAtomCount(m))
    descriptor.append(Descriptors.HeavyAtomMolWt(m))
    descriptor.append(Descriptors.NHOHCount(m))
    descriptor.append(Descriptors.NOCount(m))
    descriptor.append(Descriptors.NumHAcceptors(m))
    descriptor.append(Descriptors.NumHDonors(m))
    descriptor.append(Descriptors.NumHeteroatoms(m))
    descriptor.append(Descriptors.NumRotatableBonds(m))
    descriptor.append(Descriptors.NumValenceElectrons(m))
    descriptor.append(Descriptors.RingCount(m))
    descriptor.append(Descriptors.FractionCSP3(m))
    descriptor.append(Descriptors.NumAromaticHeterocycles(m))
    descriptor.append(Descriptors.NumSaturatedHeterocycles(m))
    descriptor.append(Descriptors.NumAliphaticHeterocycles(m))
    descriptor.append(Descriptors.NumAromaticCarbocycles(m))
    descriptor.append(Descriptors.NumSaturatedCarbocycles(m))
    descriptor.append(Descriptors.NumAliphaticCarbocycles(m))
    print(descriptor)
    return descriptor


descrptor_caculator('c1ccccc1C(=O)O')