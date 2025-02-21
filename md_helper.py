#!/usr/bin/env python

## Import
import openmm.app as openmm_app
import openmm
from openmm.unit import *

import numpy as np
import os


## Functions
def gen_cycpep_target_modellers(input, peptide_chain: int):
    """
    Generate two Modeller objects that are Target and Peptide

    PARAMS
    ------
    :input: OpenMM topology object.
    :peptide_chain: Chain number for our peptide

    RETURNS
    -------
    :target_modeller: Modeller Toplogy that has only our Non Peptide Chain/s
    :peptide_modeller: Modeller Toplogy that has only our Peptide Chain
    """
    # make a modeller of our target topology
    target_modeller = openmm.app.Modeller(input.topology, input.positions)
    peptide_modeller = openmm.app.Modeller(input.topology, input.positions)

    # define holding lists
    del_pep, del_target = list(), list()
    for i, c in enumerate(target_modeller.topology.chains()):
        for r in c.residues():
            if i != peptide_chain:
                del_target.extend(r.atoms())
            else:
                del_pep.extend(r.atoms())
    print(del_pep)
    target_modeller.delete(del_pep)
    peptide_modeller.delete(del_target)
    return target_modeller, peptide_modeller

def is_bonded(topology, atom1, atom2):
    """
    Checks if a bond exists between two given atoms in an OpenMM topology. 
    
    PARAMS
    ------
    :topology: OpenMM topology object. 
    :atom1: The first atom.
    :atom2: The second atom.
    
    RETURNS
    -------
    :bool: True if a bond exists between atom1 and atom2, False otherwise. 
    """
    for bond in topology.bonds():
        if ((bond[0].name == atom1.name) and (bond[0].residue.index == atom1.residue.index) and (bond[1].name == atom2.name) and(bond[1].residue.index == atom2.residue.index)) or \
            ((bond[1].name == atom1.name) and (bond[1].residue.index == atom1.residue.index) and (bond[0].name == atom2.name) and(bond[0].residue.index == atom2.residue.index)):
            print(bond, atom1, bond[0].name == atom1.name, bond[0].residue.index == atom1.residue.index)
            return True
    return False 

def compute_com(positions, system, atom_idx) -> np.array:
    """Compute the Center of Mass for our given atom_idx

    PARAMS
    ------
    :positions: The positions of our atoms in our topology object
    :system: Our system object
    :atom_idx: An array of atom indices over which we want to compute 
        the center of mass

    RETURNS
    -------
    com: np.array
        Center of mass x, y, z coordinates
    """
    # Extract masses
    masses = np.array([system.getParticleMass(i.index) for i in atom_idx])
    # compute com
    com = sum([m * positions[i.index] for i,m in zip(atom_idx, masses)]) / masses.sum()
    return com

def net_force(simulation, force_index: int, atom_idx: np.array) -> np.float32:
    """Sum the particular forces expressed in atom_idx within the simulation

    PARAMS
    ------
    :simulation: The entire simulation object with our protein and water and forces
    :force_index: The index for the particular force that we want to select.
    :atom_idx: Array of atoms where we want our forces summed over

    RETURNS
    -------
    :total_force: summed force value of atom_idx along the first axis
    """
    # extract forces
    print(simulation.context.getSystem().getForce(force_index))
    print(dir(simulation.context.getSystem().getForce(force_index)))
    return np.array([simulation.context.getSystem().getForce(force_index).getParticleParameters(i.index) for i in atom_idx], axis=0)

def get_target_receptor(positions, residue: int, residue_list: list) -> list:
    """Extract out the alphaCarbon atom of our residue of interest

    PARAMS
    ------
    :positions: All target positions
    :residue: Our target residue of interest (int)
    :residue_list: Target chain list of residue objects

    RETURNS
    -------
    :target_atom: The alpha carbon atom of our resiude of interest xyz
    """
    #  Extract out the resiude of interest
    target_res = residue_list[residue-1]
    print(target_res)
    target_atom = [positions[a.index] for a in target_res.atoms() if a.name == "CA"][0]
    print(target_atom)
    return target_atom

def custom_pull_force(
    com_target,
    peptide_atoms,
    peptide_com,
    pull_force_constant: float = 10.0,
    ):
    """Add an external pulling force that increases with time.

    PARAMS
    ------
    :com_target: The center of mass of the receptor.
    :peptide_atoms: All atoms that are within the peptide
    :pull_force_constant: Our constant pull force
    """
    # Set units onto our pull force 
    pull_force_constant *= kilocalories_per_mole / 1.0 * angstroms

    # We need to set our pulling point 5 angstroms away
    # pull_diff = peptide_com - com_target
    # direction = pull_diff / np.linalg.norm(pull_diff)
    # fx, fy, fz = direction

    # Apply force in only the x direction
    fx = 0.0 
    fy = 0.0
    fz = -1.0 

    # specify our force
    pullforce = openmm.CustomCentroidBondForce(1, "pull_force_constant * (x1*fx + y1*fy + z1*fz)")
    pullforce.addGlobalParameter("pull_force_constant", pull_force_constant)
    pullforce.addGlobalParameter("fx", fx)
    pullforce.addGlobalParameter("fy", fy)
    pullforce.addGlobalParameter("fz", fz)
    pullforce.addGroup(peptide_atoms)
    pullforce.addBond([0])

    return pullforce



def add_pull_force(
    com_atom,
    com_target,
    residue_target: int,
    pull_force: float = 1.0,
    ):
    """Add an external pulling force that increases with time.

    PARAMS
    ------
    :com_target: The center of mass of the receptor.
    :com_atom: All atoms that are within the peptide
    :residue_target: X,Y,Z of residue on target that interacts with our anchor
    :spring_constant: Spring constant value for our force (1.0)
    :pulling_velocity: Scaling of the velocity
    """
    # create forces 
    pull_force *= kilocalories_per_mole / angstroms

    # We need to set our pulling point 5 angstroms away
    pull_diff = com_atom - com_target
    direction = pull_diff / np.linalg.norm(pull_diff)

    # Set distance 
    distance = 15 * angstroms
    pull_point = residue_target + distance * direction

    fx, fy, fz = pull_force * direction
    # print(direction)
    # print(fx, fy, fz)

    # print(f"Difference of 10 angstrom {residue_target} -> {pull_point}")
    # print(f"Peptide COM {com_atom} and Target COM {com_target}")
    # print(f"Point Distance away from peptide: {np.sqrt(np.sum(np.power(com_atom - pull_point,2)))}")
    # print(f"Point Distance away from Target Residue: {np.sqrt(np.sum(np.power(residue_target - pull_point,2)))}")

    # Custom external force
    # force_expression = "0.5 * k * ((x1 - (x0 + v*t))^2 + (y1 - (y0 + v*t))^2 + (z1 - (z0 + v*t))^2)"
    # force_expression = "0.5 * k * ((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)"
    force_expression = "fx*(x1 - x0) + fy*(y1 - y0) + fz*(z1 - z0)"
    pulling_force = openmm.CustomExternalForce(force_expression)
    # pulling_force.addGlobalParameter("k", spring_constant)
    pulling_force.addGlobalParameter("fx", fx)
    pulling_force.addGlobalParameter("fy", fy)
    pulling_force.addGlobalParameter("fz", fz)
    pulling_force.addGlobalParameter("x0", com_atom[0])
    pulling_force.addGlobalParameter("y0", com_atom[1])
    pulling_force.addGlobalParameter("z0", com_atom[2])


    pulling_force.addGlobalParameter("x1", pull_point[0])
    pulling_force.addGlobalParameter("y1", pull_point[1])
    pulling_force.addGlobalParameter("z1", pull_point[2])
    # pulling_force.addGlobalParameter("v", pulling_velocity)
    # Add time
    # pulling_force.addGlobalParameter("t", 0.0)
    # for a in peptide_atoms:
    #     pulling_force.addParticle(a.index, modeller.positions[a.index])

    # Add our fource to the system
    # system.addForce(pulling_force)
    return pulling_force


### Classes
class ReporterPullForce(openmm_app.StateDataReporter):
    """Custom Data Reporter for OpenMM, where if we are applying SDM then we
        want the Pull Force to be reported as well"""
    def __init__(self, pullforcefile, file, reportInterval, context, pull_idx, **kwargs):
        """Init this report similarly to StateDataReporter, but with additional arguments
            for reporting our Pull Force

            PARAMS
            ------
            :file: filename to write to
            :reportInterval: The step interval to extract information
            :context: Openmm simulation context
            :pull_idx: atom indices where the pull force originates from
        """
        super().__init__(file, reportInterval, **kwargs)
        self.pullForceFile = pullforcefile
        self.context = context
        self.pull_idx = pull_idx

    def report(self, simulation, state):
        """Get normal information and then add in extra pull information"""
        super().report(simulation, state)

        # report the pulling force
        forces = self.context.getSystem().getForces()
        custom_force_index, custom_force = len([i for i in forces]), [i for i in forces]

        # compute the center of mass and the forces associated with it
        com = compute_com(simulation, self.pull_idx)
        com_forces = net_force(simulation, custom_force_index-1, self.pull_idx)
        force_mag = np.linalg.norm(com_forces)
        print("COM:", com)
        print("COM FORCES NET:", com_forces)
        print("COM FORCE MAGNITUDE:", force_mag)

        if os.path.exists(self.pullForceFile):
            with open(self.pullForceFile, 'a') as outfile:
                print(f"{simulation.currentStep}\t{force_mag:.3f}\t{com}", file=outfile)
        else:
            with open(self.pullForceFile, 'w') as outfile:
                print("Step\tPull Force\tCOM", file=outfile)
                print(f"{simulation.currentStep}\t{force_mag:.3f}\t{com}", file=outfile)
        return None





class Reporter(openmm.MinimizationReporter):
    """Sublcass that is used to report the progress of local minimization"""
    interval = 10 # report interval
    energies = [] # array to record progress
    def report(self, iteration, x, grad, args):
        # print current system energy to screen 
        if iteration % self.interval == 0:
            print(iteration, args['system energy'])

        # save energy at each iteration to an array we can use later
        self.energies.append(args['system energy'])

        return False
