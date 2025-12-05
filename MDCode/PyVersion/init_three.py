import numpy as np
from add_list_atomic import add_list_atomic

def init_three(sim):
    """
    Initialize the simulation with three atoms.

    Parameters:
    sim (MDSimulation): Simulation object
    """
    sim.do_plot = True
    sim.dt = 15e-15
    sim.t_stop = 100000 * sim.dt
    sim.init_dist = 0.0
    sim.method = 'VE'  # VE: Verlet, FD: Forward Difference

    sim.mass0 = 14 * sim.C['am']  # Silicon
    sim.mass1 = 5 * sim.C['am']   # Argon

    sim.atom_spacing = 0.5430710e-9
    sim.lj_sigma = sim.atom_spacing / 2**(1/6)
    sim.lj_epsilon = 1e-21

    sim.phi_cutoff = 3 * sim.atom_spacing * 1.1

    sim.t_target = 300

    x0 = sim.atom_spacing * np.array([-0.8, 0, 0.8])
    y0 = sim.atom_spacing * np.array([0, -1.2, 0])
    v = 0
    vx0 = np.array([0, 0, 0])
    vy0 = np.array([0, v/2, v/2])
    types = np.array([0, 1, 0])
    add_list_atomic(x0, y0, vx0, vy0, types, 0, 0, sim)

    size = 2 * sim.atom_spacing
    sim.limits = [-size, size, -size, size]
    sim.pl_delt = 5 * sim.dt
    sim.marker_size = 10
    sim.plot_file = 'Block.gif'
    sim.do_plot_image = False
    sim.plot_size = [100, 100, 1049, 1049]
    sim.plot_pos_only = True

    sim.scale_v = 0.2e-11
    sim.scale_f = 20