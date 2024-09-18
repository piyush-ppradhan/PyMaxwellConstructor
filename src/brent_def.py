import numpy as np
from scipy.optimize import NoConvergence, brent


class Brent:
    """
    Brent solver class. Note that this does not implement brent solver itself but rather stores, all the necessary information that will be used.
    Actual brent solver implementation is from scipy.optimize.brent

    Parameters:
        bounds: tuple
            Tuple of values where bounds[0] < bounds[1]
        ndivisions: int
            Number of divisions of bounds for root search
        f: lambda
            Function to be used for root check
    """

    def __init__(self, bounds, ndivisions, f):
        self.bounds = bounds
        self.ndivisions = ndivisions
        self.f = f

    def analyze_curve(self):
        signchange = 0
        sign_change_min = []
        sign_change_max = []
        x = np.linspace(self.bounds[0], self.bounds[1], self.ndivisions)
        for i in range(np.shape(x)[0] - 1):
            if self.f(x[i]) * self.f(x[i + 1]) <= 0:
                signchange += 1
                sign_change_min.append(x[i])
                sign_change_max.append(x[i + 1])
        self.nroots = signchange
        if self.nroots == 0:
            raise ValueError("No roots found for the given function")
        self.root_bounds = np.zeros((self.nroots, 2))
        self.root_bounds[:, 0] = sign_change_min
        self.root_bounds[:, 1] = sign_change_max

    def find_all_roots(self):
        self.roots = np.zeros((self.nroots,))
        for i in range(self.nroots):
            self.roots[i] = brent(
                self.f, brack=(self.root_bounds[i, 0], self.root_bounds[i, 1])
            )
