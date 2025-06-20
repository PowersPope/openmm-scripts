#!/usr/bin/env python

# import necessary openmm packages plus tools
import openmm.app as openmm_app
import openmm
from openmm.unit import *
import pdbfixer

# custom tools
from md_helper import Reporter, ReporterPullForce, compute_com, get_target_receptor, custom_pull_force

# Python Tools
import argparse
from sys import stdout
from collections import defaultdict
import os
import numpy as np

### Collect arguments
p = argparse.ArgumentParser(description="Run an Amber MD trajectory to determine macrocycle stability.  Assumption is that Chain B is macrocycle.",
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                            )
p.add_argument("--file", type=str, help="Path to Complex (target + cyclic peptide).")
p.add_argument("--output-dir", type=str, default="out/", help="Directory to store output files. Will generate dir, if not yet created.")
p.add_argument("--temperature", type=float, default=300.0, help="Temperature in Kelvin")
p.add_argument("--pressure", type=float, default=1.0, help="Atmospheric Pressure in atm")
p.add_argument("--timestep", type=float, default=2.0, help="Time step in femtoseconds")
p.add_argument("--total-steps", type=int, default=30_000, help="Total number of steps to take in Simulation")
p.add_argument("--postfix-filename", type=str, default="out", help="Specify this to make your output .pdb and .dcd unique")
p.add_argument("--add-hydrogens", action="store_true", help="Add hydrogens (Good to use if you pass --pdbfix)")
p.add_argument("--production", action="store_true", help="Specify if you want to run a production 100ns run, if not only PS and total-steps are run.")
p.add_argument("--debug", action="store_true", help="Specify if you want output to STDOUT during Production run for testing.")
p.add_argument("--cyclic", action="store_true", help="Specify if you want the peptide to be cyclic.")
p.add_argument("--pdbfix", action="store_true", help="Perform PDB cleanup first.")
p.add_argument("--smd", action="store_true", help="Perform Steered Molecular Dynamics (SMD)")
p.add_argument("--target-receptor", type=int, default=0, help="Specify Target Receptor Atom for SDM Reference")
p.add_argument("--gpu", action="store_true", help="Specify if you want to run on a GPU")
args = p.parse_args()

# generate our output path if not yet created
if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

# init a blank file name
outfile_name = ''

### Fix first
if args.pdbfix:
    fixer = pdbfixer.PDBFixer(args.file)
    fixer.findMissingResidues()
    fixer.removeHeterogens(False) # remove all heterogens including water since False
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    outfile_name = os.path.join(args.output_dir, f"{args.file.split('.pdb')[0]}-fixed.pdb")
    openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile_name, 'w'))

#### Setup
# Simulation Settings
pressure = args.pressure*atmosphere
temperature = args.temperature * kelvin
timstep = args.timestep*femtoseconds

# Load in pdb of interest
pdb = openmm_app.PDBFile(outfile_name if args.pdbfix else args.file)

print("Loading in PDB:", pdb.topology) # Print Topology before the cyclic bond

modeller = openmm_app.Modeller(pdb.topology, pdb.positions)

# Make a dict of our residues (so this is done once)
# A residue storing dict where k is the chainid, and v is a list of residue objects
residue_all = defaultdict(list)
for i, c in enumerate(modeller.topology.chains()):
    residue_all[i].extend([r for r in c.residues()])

if args.cyclic:
    print("Adding N-C Peptide Bond...")
    # Add our cyclic bond to our modeller object
    residues = residue_all[1]
    N_term = [i for i in residues[0].atoms() if i.name == "N"][0]
    C_term = [i for i in residues[-1].atoms() if i.name == "C"][0]
    modeller.topology.addBond(N_term, C_term)
    
# Grab and apply our forcefield amber14 and water forcefield
forcefield = openmm_app.ForceField("amber14-all.xml", "amber14/tip4pew.xml")

unmatched = forcefield.getUnmatchedResidues(modeller.topology)
print(unmatched)

# Add hydrogens that are missing if necessary
if args.add_hydrogens:
    modeller.addHydrogens(forcefield)
# Max our solvent explicit and water
modeller.addSolvent(forcefield, model="tip4pew", padding = 1 * nanometer, neutralize = True)

# Generate our system and our simulation
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=openmm_app.PME,
    nonbondedCutoff= 1 * nanometer,
    removeCMMotion=False,
    # constraints=openmm_app.HBonds,
)


if args.smd:
    # assert args.target_receptor!=0,"You need to set the Target Receptor Argument"
    # calculate center of mass target
    com_target = compute_com(
    modeller.positions,
    system, 
    [a for i in residue_all[0] for a in i.atoms()],
    )
    # calculate center of mass peptide
    init_com_peptide = compute_com(
    modeller.positions,
    system, 
    [a for i in residue_all[1] for a in i.atoms()],
    )
    # target_receptor_atom = get_target_receptor(modeller.positions, args.target_receptor, residue_all[0])
    custom_pull_force = custom_pull_force(
        com_target,
        [a.index for i in residue_all[1] for a in i.atoms()], 
        init_com_peptide,
        pull_force_constant = 100.0, # tried 10.0 with glycine
    )
    # custom_pull_force = custom_pull_force(
    #     com_peptide, 
    #     com_target,
    #     residue_target=target_receptor_atom,
    #     pull_force = 1.0,
    #     # spring_constant=100.,
    #     # pulling_velocity=10.,
    # )

    system.addForce(custom_pull_force)

integrator = openmm.LangevinIntegrator(temperature, 1 / picosecond, 2*femtoseconds)
simulation = openmm_app.Simulation(
    modeller.topology,
    system, 
    integrator, 
    # openmm.Platform.getPlatformByName('CUDA') if args.production else openmm.Platform.getPlatformByName('CPU'),
    openmm.Platform.getPlatformByName('CUDA') if args.gpu else openmm.Platform.getPlatformByName('CPU'),
)
    

n_platforms = openmm.Platform.getNumPlatforms()
# for i in range(n_platforms):
#     print(openmm.Platform.getPlatform(i).getName())

# Set the positions of atoms in the simulation
simulation.context.setPositions(modeller.positions)
# Create an instance of our reporter from md_helper.py
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

print("\n---------------")
print(f"Potential Energy Change After Minimization:\n{before_potential} -> {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
print("---------------\n")

if args.production:
    print("---- Production Level Run for 100ns ----")

    # set the run length
    run_length = 50_000_000 # 50,000,000 * 2 fs = 100 ns

    # Set up reporter
    simulation.reporters = []
    simulation.reporters.append(openmm_app.DCDReporter(os.path.join(args.output_dir, f"traj-{args.postfix_filename}.dcd"), 5000))
    # Log simulation
    simulation.reporters.append(
        openmm_app.StateDataReporter(
            os.path.join(args.output_dir, f"simulation-{args.postfix_filename}.csv"), 5000,
            step = True, potentialEnergy = True, temperature = True,
            progress = True, remainingTime=True, speed=True,
            totalSteps=run_length,
            separator=","
        )
    )

    if args.debug:
        # Testing purposes
        simulation.reporters.append(
            openmm_app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, elapsedTime=True)
        )

    if args.smd:
        simulation.step(run_length)

    else:
        simulation.step(run_length)

else:
    print("---- Test Level Run for ps ----")
    simulation.reporters = []
    simulation.reporters.append(openmm_app.DCDReporter(os.path.join(args.output_dir, f"traj-{args.postfix_filename}.dcd"), 10))
    simulation.reporters.append(
        openmm_app.StateDataReporter(stdout, 100, step=True, temperature=True, elapsedTime=True, potentialEnergy=True)
    )

    if args.debug:
        peptide_atoms = [a for i in residue_all[1] for a in i.atoms()]
        print(peptide_atoms)
        simulation.reporters.append(ReporterPullForce(
            os.path.join(args.output_dir,"testPullFile.log"),
            os.path.join(args.output_dir,"test_out.csv"), 
            10, 
            simulation.context, 
            peptide_atoms))


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

    forces = simulation.context.getSystem().getForces()
    custom_force_index, custom_force = len([i for i in forces]), [i for i in forces]

    if args.smd:
        for step in range(args.total_steps):
            simulation.step(1)

            com = compute_com(
            simulation.context.getState(getPositions=True).getPositions(asNumpy=True),
            simulation.context.getSystem(), 
            [a for i in residue_all[1] for a in i.atoms()],
            )
            displacement = com - init_com_peptide
            print(f"Current Step {step} -> {np.sqrt(np.sum(np.power(displacement,2)))} Displacement")

    else:
        simulation.step(args.total_steps)
