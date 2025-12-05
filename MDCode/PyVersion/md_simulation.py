import numpy as np
import matplotlib.pyplot as plt
from get_forces import get_forces
from init_three import init_three

class MDSimulation:
    def __init__(self):
        # Constants
        self.C = {
            'q_0': 1.60217653e-19,  # electron charge
            'hb': 1.054571596e-34,  # Dirac constant
            'h': 1.054571596e-34 * 2 * np.pi,  # Planck constant
            'm_0': 9.10938215e-31,  # electron mass
            'kb': 1.3806504e-23,    # Boltzmann constant
            'eps_0': 8.854187817e-12,  # vacuum permittivity
            'mu_0': 1.2566370614e-6,  # vacuum permeability
            'c': 299792458,         # speed of light
            'g': 9.80665,           # gravity
            'am': 1.66053892e-27    # atomic mass unit
        }

        # Initialize arrays
        self.x = np.array([])
        self.y = np.array([])
        self.vx = np.array([])
        self.vy = np.array([])
        self.fx = np.array([])
        self.fy = np.array([])
        self.phi = np.array([])
        self.mass = np.array([])
        self.atom_type = np.array([])
        self.n_atoms = 0

        # Other variables
        self.time = []
        self.phi_tot = []
        self.ke_tot = []
        self.temps = []
        self.t0 = []
        self.t1 = []

        # Plotting
        self.do_plot = True
        self.do_plot_image = False
        self.plot_count = 0
        self.scale_v = 0
        self.scale_f = 0
        self.marker_size = 12
        self.fig = None
        self.axes = None

    def run_simulation(self):
        # Initialization
        init_three(self)

        # Set up plotting
        if self.do_plot:
            self.fig, self.axes = plt.subplots(1, 2, figsize=(12, 4))
            plt.ion()  # Interactive mode

        # Initial setup
        max_x = np.max(self.x) * 1.5
        min_x = np.min(self.x) * 1.5
        max_y = np.max(self.y) * 1.5
        min_y = np.min(self.y) * 1.5

        self.fx = np.zeros(self.n_atoms)
        self.fy = np.zeros(self.n_atoms)
        self.phi = np.zeros(self.n_atoms)

        get_forces(self.phi_cutoff, self.lj_epsilon, self.lj_sigma, self)

        self.phi0 = np.sum(self.phi) / self.n_atoms / 2

        v2 = self.vx**2 + self.vy**2
        pty0 = self.atom_type == 0
        pty1 = self.atom_type == 1
        n_atoms0 = np.sum(pty0)
        n_atoms1 = np.sum(pty1)

        ke0 = np.mean(v2[pty0]) * self.mass0 * 0.5 if n_atoms0 > 0 else 0
        ke1 = np.mean(v2[pty1]) * self.mass1 * 0.5 if n_atoms1 > 0 else 0
        ke_tot = (ke0 * n_atoms0 + ke1 * n_atoms1)

        t = 0
        c = 0
        self.time.append(0)
        self.phi_tot.append(np.sum(self.phi) / 2)
        self.ke_tot.append(ke_tot)
        self.temps.append(ke_tot / self.n_atoms / self.C['kb'])
        self.t0.append(ke0 / self.C['kb'] if n_atoms0 > 0 else 0)
        self.t1.append(ke1 / self.C['kb'] if n_atoms1 > 0 else 0)

        # Initial positions for Verlet
        xp = self.x - self.dt * self.vx
        xpp = self.x - 2 * self.dt * self.vx
        yp = self.y - self.dt * self.vy
        ypp = self.y - 2 * self.dt * self.vy

        plt0 = self.pl_delt

        while t < self.t_stop:
            get_forces(self.phi_cutoff, self.lj_epsilon, self.lj_sigma, self)

            if self.method == 'VE':
                self.x = -xpp + 2 * xp + self.dt**2 / self.mass * self.fx
                self.y = -ypp + 2 * yp + self.dt**2 / self.mass * self.fy
                self.vx = (self.x - xpp) / (2 * self.dt)
                self.vy = (self.y - ypp) / (2 * self.dt)
            elif self.method == 'FD':
                dvx = self.fx * self.dt / self.mass
                self.vx += dvx
                dx = self.vx * self.dt + self.fx * self.dt**2 / 2 / self.mass
                self.x = xp + dx

                dvy = self.fy * self.dt / self.mass
                self.vy += dvy
                dy = self.vy * self.dt + self.fy * self.dt**2 / 2 / self.mass
                self.y = yp + dy

            xpp = xp
            ypp = yp
            xp = self.x
            yp = self.y

            c += 1
            t += self.dt
            self.time.append(t)

            self.phi_tot.append(np.sum(self.phi) / 2)
            v2 = self.vx**2 + self.vy**2
            ke0 = np.mean(v2[pty0]) * self.mass0 * 0.5 if n_atoms0 > 0 else 0
            ke1 = np.mean(v2[pty1]) * self.mass1 * 0.5 if n_atoms1 > 0 else 0
            ke_tot = (ke0 * n_atoms0 + ke1 * n_atoms1)
            self.ke_tot.append(ke_tot)
            self.temps.append(ke_tot / self.n_atoms / self.C['kb'])
            self.t0.append(ke0 / self.C['kb'])
            self.t1.append(ke1 / self.C['kb'])

            if t > plt0:
                print(f'time: {t:g} ({t / self.t_stop * 100:.2f} %)')

                if self.do_plot:
                    self.plot_vars(c)

                plt0 += self.pl_delt
                plt.pause(0.00001)

        if self.do_plot_image:
            # Save GIF if needed
            pass

    def plot_vars(self, c):
        # Clear previous plots
        for ax in self.axes.flat:
            ax.clear()

        # Positions and forces
        self.axes[0].scatter(self.x[self.atom_type == 0], self.y[self.atom_type == 0], c='b', s=self.marker_size)
        self.axes[0].scatter(self.x[self.atom_type == 1], self.y[self.atom_type == 1], c='g', s=self.marker_size)
        self.axes[0].quiver(self.x, self.y, self.fx * self.scale_f, self.fy * self.scale_f, color='m', width=0.005)
        self.axes[0].set_xlim(self.limits[0], self.limits[1])
        self.axes[0].set_ylim(self.limits[2], self.limits[3])
        self.axes[0].set_title('Atoms and Forces')
        self.axes[0].set_xlabel('X')
        self.axes[0].set_ylabel('Y')

        # Energies
        self.axes[1].plot(self.time, self.ke_tot, 'b', linewidth=2, label='KE')
        self.axes[1].plot(self.time, np.array(self.phi_tot) + self.lj_epsilon * self.n_atoms, 'g', linewidth=2, label='Phi + Epsilon')
        self.axes[1].plot(self.time, np.array(self.phi_tot) + np.array(self.ke_tot) + self.lj_epsilon * self.n_atoms, 'k', linewidth=2, label='Total')
        self.axes[1].set_xlabel('Time')
        self.axes[1].legend()

        plt.tight_layout()
        plt.draw()

if __name__ == "__main__":
    sim = MDSimulation()
    sim.run_simulation()