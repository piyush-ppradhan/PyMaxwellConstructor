import numpy as np
from scipy.integrate import quad
from scipy.optimize import NoConvergence, brent

from src.brent_def import *


class MaxwellConstructor:

    def __init__(self, eos, ndivisions=1001, tol=1e-6):
        self.eos = eos
        self.tol = tol
        self.ndivisions = ndivisions

    def get_sat_density(self, T_sat, show_output=False):
        if T_sat > self.eos.T_c:
            ValueError("Temperature is greater than critical temperature.")
        p_sat_found = False

        f = lambda rho: self.eos.EOS(rho, T_sat)
        rho_first = self.tol
        rho_end = 1.0 / self.eos.b - 100 * self.tol
        self.brent_ = Brent([rho_first, rho_end], self.ndivisions, f)
        self.brent_.analyze_curve()
        self.brent_.find_all_roots()

        rho_min = self.brent_.roots[1]
        rho_max = self.brent_.roots[0]
        p_min = self.eos.EOS(rho_min, T_sat)
        p_max = self.eos.EOS(rho_max, T_sat)

        # Choosing the pressure interval where the p_sat value can be expected to lie in.
        p_first = max(100 * self.tol, p_min)
        p_last = p_max - 10 * self.tol
        p_sat_guess = np.linspace(p_first, p_last, self.ndivisions)

        for p in p_sat_guess:
            fl = lambda rho: self.eos.EOS(rho, T_sat) - p
            fg = lambda rho: self.eos.EOS(rho, T_sat) - p
            try:
                rhog_sat = brent(fg, brack=[rho_min, rho_end])
                rhol_sat = brent(fl, brack=[rho_first, rho_max])
            except NoConvergence:
                continue
            integrand = lambda rho: p - self.eos.EOS(rho, T_sat)
            area = quad(integrand, rhog_sat, rhol_sat)
            if abs(area[0]) <= self.tol:
                p_sat_found = True
                if show_output:
                    print(f"Saturated pressure value: {p}")
                    print(f"rho_g: {rhog_sat}")
                    print(f"rho_l: {rhol_sat}")
                    return rhog_sat, rhol_sat
        if not p_sat_found:
            print("Saturated pressure for given values cannot be calculated")
