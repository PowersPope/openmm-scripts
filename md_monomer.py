#!/usr/bin/env python

# import necessary openmm packages plus tools
import openmm.app as openmm_app
import openmm
from openmm.unit import *

from md_helper import Reporter

# Python Tools
import argparse
from sys import stdout
import os
from collections import defaultdict

### Collect arguments
p = argparse.ArgumentParser(description="Run an Amber MD trajectory to determine macrocycle stability. Assumption is that Chain B is macrocycle.",
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                            )
p.add_argument("--file", type=str, help="Path to Complex (target + cyclic peptide).")
p.add_argument("--output-dir", type=str, default="out/", help="Directory to store output files. Will generate dir, if not yet created.")
p.add_argument("--temperature", type=float, default=300.0, help="Temperature in Kelvin")
p.add_argument("--pressure", type=float, default=1.0, help="Atmospheric Pressure in atm")
p.add_argument("--timestep", type=float, default=2.0, help="Time step in femtoseconds")
p.add_argument("--total-steps", type=int, default=30_000, help="Total number of steps to take in Simulation")
p.add_argument("--postfix-filename", type=str, default="out", help="Specify this to make your output .pdb and .dcd unique")
p.add_argument("--production", action="store_true", help="Specify if you want to run a production NS run, if not only PS are run.")
p.add_argument("--cyclic", action="store_true", help="Specify if you want the peptide to be cyclic.")
p.add_argument("--debug", action="store_true", help="Specify if you want output to STDOUT during Production run for testing.")
p.add_argument("--gpu", action="store_true", help="Specify if you want to run on a GPU")
args = p.parse_args()

# generate our output path if not yet created
if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

#### Setup
# Simulation Settings
pressure = args.pressure*atmosphere
temperature = args.temperature * kelvin
timstep = args.timestep*femtoseconds

# Load in pdb of interest
pdb = openmm_app.PDBFile(args.file)

print("Loading in PDB:", pdb.topology) # Print Topology before the cyclic bond

# Delete Target, as we want to test the monomers stability
modeller = openmm_app.Modeller(pdb.topology, pdb.positions)
modeller.delete([res for res in pdb.topology.residues() if res.chain.index == 0])
    
print("After only peptide selection:", modeller.topology) #print after 


if args.cyclic:
    print("Adding N-C Peptide Bond...")
    residue_all = defaultdict(list)
    for i, c in enumerate(modeller.topology.chains()):
        residue_all[i].extend([r for r in c.residues()])
    # Add our cyclic bond to our modeller object
    residues = residue_all[0]
    N_term = [i for i in residues[0].atoms() if i.name == "N"][0]
    C_term = [i for i in residues[-1].atoms() if i.name == "C"][0]
    modeller.topology.addBond(N_term, C_term)

    print("Topology Modeller:", modeller.topology)
    
# Grab and apply our forcefield amber14
forcefield = openmm_app.ForceField("amber14-all.xml", "amber14/tip4pew.xml")
# Max our solvent explicit and water
modeller.addSolvent(forcefield, model="tip4pew", padding = 1 * nanometer, neutralize = True)

# Generate our system and our simulation
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=openmm_app.PME,
    nonbondedCutoff= 1 * nanometer,
    constraints=openmm_app.HBonds,
    removeCMMotion=False,
)
integrator = openmm.LangevinIntegrator(temperature, 1 / picosecond, 2*femtoseconds)
simulation = openmm_app.Simulation(
    modeller.topology,
    system, 
    integrator, 
    openmm.Platform.getPlatformByName('CUDA') if args.gpu else openmm.Platform.getPlatformByName('CPU'),
)

n_platforms = openmm.Platform.getNumPlatforms()
for i in range(n_platforms):
    print(openmm.Platform.getPlatform(i).getName())

simulation.context.setPositions(modeller.positions)

# Create an instance of our reporter
reporter = Reporter()

# get the before energy
before_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()

# Perform local energy minimization 
print("Minimizing energy...")
simulation.minimizeEnergy(100, reporter=reporter)

simulation.minimizeEnergy() # This can be tested as I origionally had 100
positions = simulation.context.getState(getPositions=True).getPositions()
with open(os.path.join(args.output_dir, f'init-{args.postfix_filename}.pdb'), 'w') as f:
    openmm_app.PDBFile.writeFile(simulation.topology, positions, f)

print("After init topology:", simulation.topology)

print("\n---------------")
print(f"Potential Energy Change After Minimization:\n{before_potential} -> {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
print("---------------\n")

if args.production:
    print("---- Production Level Run for 100ns ----")

    # set the run length
    run_length = 50_000_000 # 50000000 * 2 fs = 100 ns

    # Set up reporter
    simulation.reporters = []
    # size stability
    simulation.reporters.append(openmm_app.DCDReporter(os.path.join(args.output_dir, f"traj-{args.postfix_filename}.dcd"), 5000))
    # Log simulation
    simulation.reporters.append(
        openmm_app.StateDataReporter(
            os.path.join(args.output_dir, f"simulation-{args.postfix_filename}.log"), 5000,
            step = True, potentialEnergy = True, temperature = True,
            progress = True, remainingTime=True, speed=True,
            totalSteps=run_length, time = True,
            separator=","
        )
    )

    if args.debug:
        # Testing purposes
        simulation.reporters.append(
            openmm_app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, elapsedTime=True)
        )

    simulation.step(run_length)

else:
    print("---- Test Level Run for 200 ps ----")
    simulation.reporters = []
    simulation.reporters.append(openmm_app.DCDReporter(os.path.join(args.output_dir, f"traj-{args.postfix_filename}.dcd"), 10))
    simulation.reporters.append(
        openmm_app.StateDataReporter(stdout, 100, step=True, temperature=True, elapsedTime=True)
    )

    simulation.reporters.append(
        openmm_app.StateDataReporter(
            os.path.join(args.output_dir, f"scalars-{args.postfix_filename}.csv"),
            10,
            step=True,
            time=True,
            potentialEnergy=True,
            totalEnergy=True,
            temperature=True,
        )
    )
    simulation.step(args.total_steps)
