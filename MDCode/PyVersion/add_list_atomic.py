import numpy as np

def add_list_atomic(x0, y0, vx0, vy0, types, init_dist, temp, sim):
    """
    Add atoms to the simulation.

    Parameters:
    x0, y0 (array): Initial positions
    vx0, vy0 (array): Initial velocities
    types (array): Atom types (0 or 1)
    init_dist (float): Initial distribution factor
    temp (float): Temperature for thermal initialization
    sim (MDSimulation): Simulation object
    """
    num_atoms = len(x0)
    xp = np.array(x0)
    yp = np.array(y0)

    # Add random displacement
    sim.x = np.concatenate([sim.x, xp + (np.random.rand(num_atoms) - 0.5) * sim.atom_spacing * init_dist])
    sim.y = np.concatenate([sim.y, yp + (np.random.rand(num_atoms) - 0.5) * sim.atom_spacing * init_dist])

    # Masses
    mass = np.ones(num_atoms) * sim.mass0
    mass[np.array(types) == 1] = sim.mass1
    sim.mass = np.concatenate([sim.mass, mass])

    sim.atom_type = np.concatenate([sim.atom_type, np.array(types)])

    if temp == 0:
        vx = np.zeros(num_atoms)
        vy = np.zeros(num_atoms)
    else:
        std = np.sqrt(sim.C['kb'] * temp / mass)
        vx = std * np.random.randn(num_atoms)
        vy = std * np.random.randn(num_atoms)

    # Center velocities and add initial
    vx = vx - np.mean(vx) + np.array(vx0)
    vy = vy - np.mean(vy) + np.array(vy0)

    sim.vx = np.concatenate([sim.vx, vx])
    sim.vy = np.concatenate([sim.vy, vy])

    sim.n_atoms += num_atoms