# https://www.rdkit.org/
#https://github.com/rdkit/rdkit
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors


import pandas as pd
import numpy as np
from mordred import Calculator, descriptors

# https://zivgitlab.uni-muenster.de/m_kueh11/fp-dm-tool
# It was an excel file, but converted it to csv file 
# It is the energy gap between HOMO-LUMO
url = 'https://raw.githubusercontent.com/gashawmg/molecular-descriptors/main/Orbital_Energies_input_data.csv'
dataset = pd.read_csv(url)
print(dataset.shape)
print(dataset.head())

# There might be one or more valid SMILES that can represent one compound
# Thanks to Pat Walters for this information,checkout his excellent blog: https://www.blogger.com/profile/18223198920629617711
def canonical_smiles(smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles] 
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    return smiles
from rdkit.Chem import Draw
import IPython
from rdkit.Chem.Draw import IPythonConsole

print(Chem.MolFromSmiles('C=CCC'))
print(Chem.MolFromSmiles('CCC=C'))
print(canonical_smiles(['C=CCC']))
print(canonical_smiles(['CCC=C']))


# Canonical SMILES
Canon_SMILES = canonical_smiles(dataset.SMILES)
print(len(Canon_SMILES))

# Put the smiles in the dataframe
dataset['SMILES'] = Canon_SMILES
print (dataset)
# Create a list for duplicate smiles
duplicates_smiles = dataset[dataset['SMILES'].duplicated()]['SMILES'].values
len(duplicates_smiles)

# Create a list for duplicate smiles
dataset[dataset['SMILES'].isin(duplicates_smiles)].sort_values(by=['SMILES'])

dataset_new = dataset.drop_duplicates(subset=['SMILES'])
len(dataset_new)

print (dataset_new)


#Calculate descriptors using RDkit
#General molecular descriptors-about 200 molecular descriptors

def RDkit_descriptors(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles] 
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    
    Mol_descriptors =[]
    for mol in mols:
        # add hydrogens to molecules
        mol=Chem.AddHs(mol)
        # Calculate all 200 descriptors for each molecule
        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors,desc_names 

# Function call
Mol_descriptors,desc_names = RDkit_descriptors(dataset_new['SMILES'])

df_with_200_descriptors = pd.DataFrame(Mol_descriptors,columns=desc_names)
print(df_with_200_descriptors)