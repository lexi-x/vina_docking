import argparse
from vina import Vina
import os
import glob
import logging

def main():
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

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler('my_log.log')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    for file in files: 
        logger.info(f"Docking {file} against {args.receptor} with grid box centered at {args.center} of size {args.box} with exhaustiveness of {args.exhaustiveness} and {args.n_poses} poses.")        
        try:
            print(f"Currently docking: {file}")
            v.set_ligand_from_file(file)
            v.compute_vina_maps(center=args.center, box_size=args.box)
        except:
            error_message = "Error, please input accurate center and box coordinates"
            print(error_message)
            logger.info(error_message)
            continue
        
        # Minimized locally the current pose
        try:
            energy_minimized = v.optimize()
            print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
            v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)
        except:
            error_message = "Error minimizing energy, please re-specify inputs"
            print(error_message)
            logger.info(error_message)
            continue

        # Dock the ligand
        try:
            v.dock(args.exhaustiveness, args.n_poses)
            filepath = file.split('/')
            filename=filepath[-1][:-6]+'_out.pdbqt'
            if not os.path.exists(args.results_dir):
                os.makedirs(args.results_dir)
            v.write_poses(f'{args.results_dir}/{filename}', n_poses=1, overwrite=True)
            logger.info("Docking successful.")           
        except:
            error_message = "Error docking ligand, please re-specify inputs"
            print(error_message)
            logger.info(error_message)
            continue



if __name__ == "__main__":
    main()
