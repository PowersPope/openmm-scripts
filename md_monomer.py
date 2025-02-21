#!/usr/bin/env python

# import necessary openmm packages plus tools
import openmm.app as openmm_app
import openmm
from openmm.unit import *
import pdbfixer

from md_helper import Reporter, is_bonded, gen_cycpep_target_modellers

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
p.add_argument("--total-steps-nonproduction", type=int, default=30_000, help="Total number of steps to take in the nonproduction Simulation (Default is 200ps)")
p.add_argument("--total-steps-production", type=int, default=50_000_000, help="Total number of steps to take in the production Simulation (Default is 100ns = 50,000,000 * 2fs)")
p.add_argument("--production-report-every-n", type=int, default=5000, help="Every designated time step it will take a frame and report it.")
p.add_argument("--postfix-filename", type=str, default="out", help="Specify this to make your output .pdb and .dcd unique")
p.add_argument("--production", action="store_true", help="Specify if you want to run a production NS run, if not only PS are run.")
p.add_argument("--cyclic", action="store_true", help="Specify if you want the peptide to be cyclic.")
p.add_argument("--debug", action="store_true", help="Specify if you want output to STDOUT during Production run for testing.")
p.add_argument("--gpu", action="store_true", help="Specify if you want to run on a GPU")
p.add_argument("--run-without-changes", action="store_true", help="Specify if you want to simply run something without deleting the target")
p.add_argument("--peptide-chain", type=int, default=0, help="Chain (int) that is our peptide of interest (if applying --cyclic)")
p.add_argument("--pdbfix", action="store_true", help="Perform PDB cleanup first.")
args = p.parse_args()

# generate our output path if not yet created
if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

# init a blank file name
outfile_name = ''

id_dict = dict(
    zip(
        [0,1,2,3,4,5,6,7,8,9,10],
        ['A','B','C','D','E','F','G','H','I','J']
    )
)

### Fix first
if args.pdbfix:
    fixer = pdbfixer.PDBFixer(args.file)
    fixer.findMissingResidues()
    print(fixer.missingResidues)
    if not args.cyclic:
        for chain in fixer.topology.chains():
            if chain.id == id_dict[args.peptide_chain]:
                lastIndexInChain = [i for i, res in enumerate(chain.residues())][-1]
                fixer.missingResidues[(chain.index, lastIndexInChain + 1)] = ["NME"]
                fixer.missingResidues[(chain.index, 0)] = ["ACE"]
    print(fixer.missingResidues)
    # fixer.removeHeterogens(True) # remove all heterogens including water since False
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # fixer.addMissingHydrogens(7.0)
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

if not args.run_without_changes:
    # Delete Target, as we want to test the monomers stability
    # modeller = openmm_app.Modeller(pdb.topology, pdb.positions)
    # if len([c for c in modeller.topology.chains()]) > 1:
    #     modeller.delete([res for res in pdb.topology.residues() if res.chain.index != args.peptide_chain])
    modeller, peptide_modeller = gen_cycpep_target_modellers(pdb, args.peptide_chain)
        
    print("After only peptide selection:", modeller.topology) #print after 
elif args.run_without_changes and args.cyclic:
    modeller, peptide_modeller = gen_cycpep_target_modellers(pdb, args.peptide_chain)
    print("TARGET MODELLER:", modeller.topology)
    print("PEPTIDE MODELLER:", peptide_modeller.topology)
else:
    modeller = openmm_app.Modeller(pdb.topology, pdb.positions)

true_pep_idx = args.peptide_chain if args.run_without_changes else 0

if args.cyclic:
    print("Adding N-C Peptide Bond...")
    delete_atoms = ["OXT"]
    # residue_all = defaultdict(list)
    # for i, c in enumerate(modeller.topology.chains()):
        # residue_all[i].extend([r for r in c.residues()])
    residues = [r for r in peptide_modeller.topology.residues()]
    # Add our cyclic bond to our modeller object
    # residues = residue_all[true_pep_idx]
    print("Topology Peptide Modeller:", peptide_modeller.topology)
    d_idx = [i for i in residues[-1].atoms() if i.name in delete_atoms]
    N_term = [i for i in residues[0].atoms() if i.name == "N"][0]
    C_term = [i for i in residues[-1].atoms() if i.name == "C"][0]
    peptide_modeller.topology.addBond(N_term, C_term)
    peptide_modeller.delete(d_idx)
    print("Topology Modeller:", peptide_modeller.topology)


# print("BEFORE SOLVENT:", modeller.topology)

# print([i for i in modeller.topology.residues()])
    
# Grab and apply our forcefield amber14
forcefield = openmm_app.ForceField("amber14-all.xml", "amber14/tip4pew.xml")
    # Might add an explicit delete call here where it deletes only the N hydrogens, so that I don't have to make this conditional
if args.cyclic and args.run_without_changes:
    modeller.addHydrogens(forcefield)

modeller.add(peptide_modeller.topology, peptide_modeller.positions)
print("Modeller FINAL:", modeller.topology)
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
    run_length = args.total_steps_production # default:50000000 * 2 fs = 100 ns

    # Set up reporter
    simulation.reporters = []
    # size stability
    simulation.reporters.append(openmm_app.DCDReporter(os.path.join(args.output_dir, f"traj-{args.postfix_filename}.dcd"), args.production_report_every_n))
    # Log simulation
    simulation.reporters.append(
        openmm_app.StateDataReporter(
            os.path.join(args.output_dir, f"simulation-{args.postfix_filename}.log"), args.production_report_every_n,
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
    simulation.step(args.total_steps_nonproduction)
