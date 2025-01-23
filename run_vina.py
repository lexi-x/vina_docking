import argparse
from vina import Vina
import os
import glob

parser = argparse.ArgumentParser( prog='vina_inputs', description='Get inputs for Autodock Vina')
parser.add_argument('ligand_directory',help='directory of ligand files')
parser.add_argument('receptor',help='receptor file')
parser.add_argument('center', type=float, nargs=3, help='center coordinates')
parser.add_argument('box', type=float, nargs=3, help='box coordinates')
parser.add_argument('--exhaustiveness', type=int, default=32, required=False, help='specify exhaustiveness settings')
parser.add_argument('--n_poses', type=int, default=20, required=False, help='specify number of poses')
parser.add_argument('results_dir',default='.', help='specify directory to place docking results')
args = parser.parse_args()

files = glob.glob(f'{args.ligand_directory}/*.pdbqt',recursive=True)

v = Vina(sf_name='vina')

v.set_receptor(args.receptor)

for file in files: 
    print(f"Currently docking: {file}")
    v.set_ligand_from_file(file)
    v.compute_vina_maps(center=args.center, box_size=args.box)

    # Score the current pose
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

    # Dock the ligand
    v.dock(args.exhaustiveness, args.n_poses)
    filepath = file.split('/')
    filename=filepath[-1][:-6]+'_out.pdbqt'
    if not os.path.exists(args.results_dir):
        os.makedirs(args.results_dir)
    v.write_poses(f'{args.results_dir}/{filename}', n_poses=1, overwrite=True)

