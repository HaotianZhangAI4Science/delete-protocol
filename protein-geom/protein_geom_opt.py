import argparse
import openmm as mm
from openmm import app
from openmm.unit import kelvin, picosecond, femtosecond
import os.path as osp
# conda install -c conda-forge openmm

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pdb', type=str, required=True)
    parser.add_argument('--out_name', type=str, required=True)
    args = parser.parse_args()

    # 1. Load the Structure
    pdb = app.PDBFile(args.input_pdb)

    # 2. Select a Force Field
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # 3. Create a System
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
                                    nonbondedCutoff=0.05*mm.unit.nanometer,
                                    constraints=app.HBonds, rigidWater=True)

    # 4. Set Up an Integrator
    integrator = mm.LangevinIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtosecond)

    # 5. Create a Simulation
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    # 6. Perform Energy Minimization
    simulation.minimizeEnergy(maxIterations=1000)  # You can set the maximum number of iterations

    # 7. Extract and Save the Minimized Coordinates
    minimized_positions = simulation.context.getState(getPositions=True).getPositions()
    
    if args.out_name:
        out_file = osp.join(osp.dirname(args.input_pdb), args.out_name)
    else:
        out_file = osp.join(osp.dirname(args.input_pdb), 'geom_opt.pdb')
    
    app.PDBFile.writeFile(simulation.topology, minimized_positions, open(out_file, 'w'))
