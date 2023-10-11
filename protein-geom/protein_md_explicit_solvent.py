# conda install -c conda-forge openmm
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os.path as osp
from pdbfixer import PDBFixer
import mdtraj as md
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pdb', type=str, required=True)
    parser.add_argument('--out_name', type=str, required=True)
    parser.add_argument('--equil_steps', type=int, default=100)
    parser.add_argument('--total_steps', type=int, default=10000)
    parser.add_argument('--record_steps', type=int, default=None,
                help='if None, record_steps = int(total_steps / 1000)')
    parser.add_argument('--device', type=str, default='CUDA')
    args = parser.parse_args()

    # 0. Solvent the Structure
    fixer = PDBFixer(filename = args.input_pdb)

    # Add missing atoms
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0) 

    # Define the box dimensions (Optional: if you want to define the box size)
    # The box vectors define the dimensions of the box and its angles.
    box_vectors = [unit.Quantity((8, 0, 0), unit.nanometers),
                unit.Quantity((0, 8, 0), unit.nanometers),
                unit.Quantity((0, 0, 8), unit.nanometers)]
    pdb_dir = osp.dirname(args.input_pdb)
    
    # Solvate the system
    # If you defined the box size, add it as an argument: fixer.addSolvent(boxSize=box_vectors)
    fixer.addSolvent(padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)
    # Save the solvated structure
    with open(osp.join(pdb_dir, 'solvated.pdb'), 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # 1. Load the Solved Structure
    pdb = PDBFile(osp.join(pdb_dir, 'solvated.pdb'))
    # 2. Select a Force Field
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers, constraints=HBonds, rigidWater=True,
        ewaldErrorTolerance=0.0005)

    integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    platform = Platform.getPlatformByName(args.device)
    prop = {'CudaPrecision': 'mixed'}

    simulation = Simulation(pdb.topology, system, integrator, platform, prop)
    simulation.context.setPositions(pdb.positions)
    print('Minimizing...')
    simulation.minimizeEnergy()

    if args.equil_steps:
        # Equilibrate
        print('Equilibrating...')
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        simulation.step(args.equil_steps)  # Adjust the number of steps as needed

    if args.record_steps is None:
        record_steps = int(args.total_steps / 1000)
    else:
        record_steps = args.record_steps

    simulation.reporters.append(DCDReporter(osp.join(pdb_dir, 'md_trajectory.dcd'), record_steps))  # Save trajectory
    simulation.reporters.append(StateDataReporter(stdout, record_steps, step=True,
            potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
            speed=True, totalSteps=args.total_steps, separator='\t'))  # Print progress to stdout
    
    print('Running Production...')
    simulation.step(args.total_steps)  # Adjust the number of steps as needed
    print('Done!')

    # Extract and Save the Protein-only Conformations
    traj = md.load(osp.join(pdb_dir, 'md_trajectory.dcd'), top=osp.join(pdb_dir, 'solvated.pdb'))
    traj = traj.superpose(reference=traj[0], atom_indices=traj.topology.select('protein'))
    protein_traj = traj.atom_slice(traj.topology.select('protein'))
    if args.out_name:
        out_file = osp.join(pdb_dir, args.out_name)
    else:
        out_file = osp.join(pdb_dir, 'centered_md_traj.pdb')
        
    protein_traj.save_pdb(out_file)

