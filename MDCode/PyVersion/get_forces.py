import numpy as np
from lj_pot import lj_pot

def get_forces(phi_cutoff, epsilon, sigma, sim):
    """
    Compute forces and potentials for all atoms.

    Parameters:
    phi_cutoff (float): Cutoff distance for potential
    epsilon (float): LJ epsilon
    sigma (float): LJ sigma
    sim (MDSimulation): Simulation object
    """
    sim.fx = np.zeros(sim.n_atoms)
    sim.fy = np.zeros(sim.n_atoms)
    sim.phi = np.zeros(sim.n_atoms)

    for i in range(sim.n_atoms):
        for j in range(sim.n_atoms):
            if i == j:
                continue
            dx = sim.x[i] - sim.x[j]
            dy = sim.y[i] - sim.y[j]
            r = np.sqrt(dx**2 + dy**2)

            if r > phi_cutoff:
                continue

            a_phi, dphi_dr = lj_pot(r, epsilon, sigma)

            ang = np.arctan2(dy, dx)
            dfx = -dphi_dr * np.cos(ang)
            dfy = -dphi_dr * np.sin(ang)

            sim.phi[i] += a_phi
            sim.fx[i] += dfx
            sim.fy[i] += dfy