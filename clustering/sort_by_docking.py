import pandas as pd
import os.path as osp
from glob import glob

from rdkit.Chem import Draw
from copy import deepcopy
from PIL import Image
import io

from rdkit import DataStructs
from utils import compute_sims, generalize, find_match, sort_lists_by_first_list, imgs2singlePDF, read_sdf, write_sdf
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--docking_base', type=str, default='./docking_base/frag2_chemdiv_searched')
    parser.add_argument('--query_mols', type=str, default='./tmp/query_mols.sdf')
    parser.add_argument('--ID_marker', type=str, default='IDNUMBER')
    parser.add_argument('--out_pdf', type=str, default=None)
    args = parser.parse_args()

    sdf_qvina_files = glob(osp.join(args.docking_base, '*.sdf'))
    qvina_files = glob(osp.join(args.docking_base, '*qvina.sdf'))
    sdf_files = list(set(sdf_qvina_files) - set(qvina_files))
    query_mols = read_sdf(args.query_mols)

    id_number = []
    show_mols = []
    affins = []
    for i in range(len(sdf_files)):
        try: 
            dock_file = osp.join(args.docking_base, f'{i}.sdf')
            mol = read_sdf(dock_file)[0]
            docked_file = osp.join(args.docking_base, f'{i}_qvina.sdf')
            docked_mol = read_sdf(docked_file)[0]

            id_number.append(mol.GetProp(args.ID_marker))
            show_mols.append(mol)
            affins.append(docked_mol.GetProp('REMARK').strip().split()[2])
        except Exception as e:
            print(e)

    affins, show_mols, id_number= sort_lists_by_first_list([affins, show_mols, id_number], ascending=False)
    print('show mols:', len(show_mols))


    show_mols = deepcopy(show_mols)
    [i.RemoveAllConformers() for i in show_mols] # remove conformers
    number = list(range(len(show_mols)))
    legends = [f'No.{i[0]}; Dock: {i[1]}; ID: {i[2]};' for i in zip(number, affins, id_number)]
    highlight_list = [find_match(mol, query_mols) for mol in show_mols] 

    # group molecules for merging
    group_size = 52  
    show_mols_group = [show_mols[i:i + group_size] for i in range(0, len(show_mols), group_size)]
    highlight_list_group = [highlight_list[i:i + group_size] for i in range(0, len(highlight_list), group_size)]
    legends_group = [legends[i:i + group_size] for i in range(0, len(legends), group_size)]

    imgs = []
    for show_mols, legend, highlights in zip(show_mols_group, legends_group, highlight_list_group):
        print(1)
        print(highlights)
        imgs.append(Draw.MolsToGridImage(show_mols, molsPerRow=4, subImgSize=(500,500), legends=legend, highlightAtomLists=highlights,maxMols=group_size))
    
    if args.out_pdf is None:
        out_pdf = osp.basename(args.docking_base) + '.pdf'
    else:
        out_pdf = args.out_pdf

    imgs2singlePDF(imgs, out_pdf)
    