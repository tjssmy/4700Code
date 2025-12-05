# Python Version of MD Simulation

This is a Python translation of the MATLAB molecular dynamics simulation code.

## Files
- `md_simulation.py`: Main simulation class and runner.
- `lj_pot.py`: Lennard-Jones potential and derivative.
- `get_forces.py`: Force calculation.
- `add_list_atomic.py`: Add atoms to simulation.
- `init_three.py`: Initialize with three atoms.
- `requirements.txt`: Dependencies.

## Usage
1. Install dependencies: `pip install -r requirements.txt`
2. Run: `python md_simulation.py`

## Notes
- Translated from MATLAB to Python using NumPy and Matplotlib.
- Uses Verlet integration by default.
- Plots positions, forces, and energies in real-time.