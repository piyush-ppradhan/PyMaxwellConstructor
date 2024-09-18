import numpy as np
from scipy.integrate import quad

from src.brent_def import Brent
from src.eos import *


# from src.utils import generate_intervals
class MaxwellIntegrator:
    """
    Similar to the Brent class, this class stores all the necessary information for the Maxwell integral.
    The actual implementation for the integration function is from scipy.integrate.quad

    brent_: Brent
        Brent class object
    bounds: list, float
        Bounds used for integration
    eos: EOS
        Equation of state being used
    error_tol: float
        Tolerance used for error computation
    max_iter: int
        Default: 1000
    ndivisions: int
        Default: 1000
    """

    def __init__(self, bounds, eos, error_tol, max_iter=1000, ndivisions=1000):
        self.bounds = bounds
        self.eos = eos
        self.error_tol = error_tol
        self.max_iter = max_iter
        self.ndivisions = ndivisions
        self.brent_ = Brent(bounds, ndivisions, eos.EOS)
        self.dp_drho = self.eos.dp_drho
        self.d2p_drho2 = self.eos.d2p_drho2

        self.brent_.analyze_curve()
        self.brent_.find_all_roots()

        self.rho_min = np.min(self.brent_.roots)
        self.rho_max = np.max(self.brent_.roots)

        self.p_min = self.eos.EOS(self.rho_min)
        self.p_max = self.eos.EOS(self.rho_max)
