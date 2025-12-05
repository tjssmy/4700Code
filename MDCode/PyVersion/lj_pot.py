def lj_pot(r, epsilon, sigma):
    """
    Compute Lennard-Jones potential and its derivative.

    Parameters:
    r (float): Distance between atoms
    epsilon (float): Depth of the potential well
    sigma (float): Finite distance at which the inter-particle potential is zero

    Returns:
    tuple: (phi, dphi_dr) where phi is the potential and dphi_dr is the derivative
    """
    phi = 4 * epsilon * (sigma**12 * r**(-12) - sigma**6 * r**(-6))
    dphi_dr = 4 * epsilon * (sigma**12 * (-12) * r**(-13) - sigma**6 * (-6) * r**(-7))
    return phi, dphi_dr