# export OE_LICENSE=/home/haotian/.OpenEye/oe_license.txt
# 
from pydocking.docking.chem import *
from pydocking.docking.qvina import prepare_ligand, prepare_target, docking_with_qvina02
from glob import glob
from rdkit import RDLogger
from rdkit.Chem import AllChem
import argparse
from tqdm import tqdm
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
def opt(sdf_file):
    mol = read_sdf(sdf_file, sanitize=True)[0]
    mol.RemoveAllConformers()
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    write_sdf([mol], sdf_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--docking_base', type=str, default='./docking_base/frag2_chemdiv_searched_bondconstraint')
    parser.add_argument('--ori_file', type=str, default='./KRD4.sdf')
    parser.add_argument('--pdb_file', type=str, default='./BRD4_protein.pdb')
    parser.add_argument('--num_sdfs', type=int, default=100)
    args = parser.parse_args()

    protein_pdbqt = prepare_target(args.pdb_file)
    center = sdf2centroid(args.ori_file)
    ori_pdbqt = prepare_ligand(args.ori_file)


    for i in tqdm(range(args.num_sdfs)):
        dock_file = osp.join(args.docking_base, f'{i}.sdf')

        ligand_pdbqt = prepare_ligand(dock_file, verbose=False)
        try:
            dock_out = docking_with_qvina02(protein_pdbqt, ligand_pdbqt, center, verbose=False) 
        except Exception as e:
            print(e)
            print(ligand_pdbqt, 'docking failed')
    