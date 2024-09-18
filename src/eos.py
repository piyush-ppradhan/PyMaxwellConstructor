"""
Definition for all Equation of States
"""


class EOS:

    def __init__(self):
        pass


class CarnahanStarling(EOS):

    def __init__(self, **kwargs):
        self.a = kwargs.get("a")
        self.b = kwargs.get("b")
        self.R = kwargs.get("R")
        self.set_critical_parameters()

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        if value > 0:
            self._a = value

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        if value > 0:
            self._b = value

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if value > 0:
            self._R = value

    def EOS(self, rho, T):
        x = 0.25 * self.b * rho
        return rho * T * (1 + x + x**2 - x**3) / (1 - x) ** 3 - self.a * rho**2

    def dp_drho(self, rho, T):
        BB = self.b / 4
        x = BB * rho

        return (
            T * (1.0 + x + x**2 - x**3) / (1.0 - x) ** 3
            + x / BB * T * (BB + 2.0 * BB * x - 3.0 * BB * x**2) / (1.0 - x) ** 3
            + 3 * x * T * (1.0 + x + (x**2) - x**3) / (1.0 - x) ** 4
            - 2.0 * self.a * x / BB
        )

    def d2p_drho2(self, rho, T):
        BB = self.b / 4
        x = BB * rho

        return (
            2.0 * T * (BB + 2.0 * BB * x - 3.0 * BB * x**2) / (1.0 - x) ** 3
            + 6.0 * T * (1.0 + x + x**2 - x**3) / (1.0 - x) ** 4 * BB
            + x / BB * T * (2.0 * BB**2 - 6.0 * BB**2 * x) / (1.0 - x) ** 3
            + 6.0 * x * T * (BB + 2 * BB * x - 3.0 * BB * x**2) / (1.0 - x) ** 4
            + 12.0 * x * BB * T * (1.0 + x + x**2 - x**3) / (1.0 - x) ** 5
            - 2.0 * self.a
        )

    def set_critical_parameters(self):
        r = 0.52177553676981578055344660
        ca = 0.49638805772940987326931762
        cb = 0.18729456694673304069903652

        self.T_c = cb * self.a / (ca * self.b)
        self.p_c = cb / self.b * self.T_c
        self.rho_c = r / self.b


class PengRobinson(EOS):

    def __init__(self, **kwargs):
        self.a = kwargs.get("a")
        self.b = kwargs.get("b")
        self.omega = kwargs.get("omega")
        self.R = kwargs.get("R")
        self.set_critical_parameters()

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        if value > 0:
            self._a = value

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        if value > 0:
            self._b = value

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, value):
        if value > 0:
            self._omega = value

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if value > 0:
            self._R = value

    def alpha(self, T):
        return (
            1.0
            + (0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2)
            * (1.0 - (T / self.T_c) ** 0.5)
        ) ** 2

    def EOS(self, rho, T):
        alpha = self.alpha(T)
        return (rho * self.R * T) / (1.0 - self.b * rho) - (
            self.a * alpha * rho**2
        ) / (1.0 + 2 * self.b * rho - (self.b * rho) ** 2)

    def dp_drho(self, rho, T):
        alpha = self.alpha(T)
        a = self.a
        b = self.b
        return (
            T / (1.0 - b * rho)
            + rho * T / (1.0 - b * rho) ** 2 * b
            - 2.0 * alpha * a * rho / (1.0 + 2.0 * b * rho - b**2 * rho**2)
            + alpha
            * a
            * rho**2
            / (1.0 + 2.0 * b * rho - b**2 * rho**2) ** 2
            * (2.0 * b - 2.0 * b**2 * rho)
        )

    def d2p_drho2(self, rho, T):
        a = self.a
        b = self.b
        alpha = self.alpha(T)
        return (
            2.0 * T / (1.0 - b * rho) ** 2 * b
            + 2.0 * rho * T / (1.0 - b * rho) ** 3 * b**2
            - 2.0 * alpha * a / (1.0 + 2.0 * b * rho - b**2 * rho**2)
            + 4.0
            * alpha
            * a
            * rho
            / (1.0 + 2.0 * b * rho - b**2 * rho**2) ** 2
            * (2.0 * b - 2.0 * b**2 * rho)
            - 2.0
            * alpha
            * a
            * rho**2
            / (1.0 + 2.0 * b * rho - b**2 * rho**2) ** 3
            * (2.0 * b - 2.0 * b**2 * rho) ** 2
            - 2.0
            * alpha
            * a
            * rho**2
            / (1.0 + 2.0 * b * rho - b**2 * rho**2) ** 2
            * b**2
        )

    def set_critical_parameters(self):
        rc = 0.25307658654159946227082744
        c1 = 0.07779607390388845597184472
        c2 = 0.45723552892138218938346024

        self.p_c = self.a / self.b**2 * c1**2 / c2
        self.T_c = c1 * self.a / (c2 * self.b)
        self.rho_c = rc / self.b


class RedlichKwong(EOS):

    def __init__(self, **kwargs):
        self.a = kwargs.get("a")
        self.b = kwargs.get("b")
        self.R = kwargs.get("R")
        self.set_critical_parameters()

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        if value > 0:
            self._a = value

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        if value > 0:
            self._b = value

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if value > 0:
            self._R = value

    def EOS(self, rho, T):
        return (rho * self.R * T) / (1.0 - self.b * rho) - (self.a * rho**2) / (
            (T**0.5) * (1.0 + self.b * rho)
        )

    def dp_drho(self, rho, T):
        a = self.a
        b = self.b
        return (
            T / (1.0 - b * rho)
            + rho * T * b / (1.0 - b * rho) ** 2
            - 2.0 * a * rho * T ** (-0.5) / (1.0 + b * rho)
            + a * rho**2 * b * T ** (-0.5) / (1.0 + b * rho) ** 2
        )

    def d2p_drho2(self, rho, T):
        a = self.a
        b = self.b
        return (
            2.0 * T * b / (1.0 - b * rho) ** 2
            + 2.0 * rho * T * b**2 / (1.0 - b * rho) ** 3
            - 2.0 * a * T ** (-0.5) / (1.0 + b * rho)
            + 4.0 * a * rho * b * T ** (-0.5) / (1.0 + b * rho) ** 2
            - 2.0 * a * rho**2 * b**2 * T ** (-0.5) / (1.0 + b * rho) ** 3
        )

    def set_critical_parameters(self):
        rc = 0.25992104989487316476721061
        c1 = 0.08664034996495772158907019
        c2 = 0.42748023354034140439099065

        self.T_c = (c1 * self.a / (c2 * self.b)) ** (2.0 / 3.0)
        self.p_c = self.T_c * c1 / self.b
        self.rho_c = rc / self.b


class RedlichKwongSoave(EOS):

    def __init__(self, **kwargs):
        self.a = kwargs.get("a")
        self.b = kwargs.get("b")
        self.omega = kwargs.get("omega")
        self.R = kwargs.get("R")
        self.set_critical_parameters()

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        if value > 0:
            self._a = value

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        if value > 0:
            self._b = value

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, value):
        if value > 0:
            self._omega = value

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if value > 0:
            self._R = value

    def alpha(self, T):
        return (
            1.0
            + (0.48 + 1.574 * self.omega - 0.176 * self.omega**2)
            * (1.0 - (T / self.T_c) ** 0.5)
        ) ** 2

    def EOS(self, rho, T):
        return (rho * self.R * T) / (1.0 - self.b * rho) - (
            self.a * self.alpha * rho**2
        ) / (1.0 + self.b * rho)

    def dp_drho(self, rho, T):
        alpha = self.alpha(T)
        a = self.a
        b = self.b
        return (
            T / (1 - b * rho)
            + rho * T / (1 - b * rho) ** 2 * b
            - 2 * alpha * a * rho / (1 + b * rho)
            + alpha * a * rho**2 / (1 + b * rho) ** 2 * b
        )

    def d2p_drho2(self, rho, T):
        alpha = self.alpha(T)
        a = self.a
        b = self.b
        return (
            2 * T / (1 - b * rho) ** 2 * b
            + 2 * rho * T / (1 - b * rho) ** 3 * b**2
            - 2 * alpha * a / (1 + b * rho)
            + 4 * alpha * a * rho / (1 + b * rho) ** 2 * b
            - 2 * alpha * a * rho**2 / (1 + b * rho) ** 3 * b**2
        )

    def set_critical_parameters(self):
        rc = 0.25992104989487316476721061
        ca = 0.42748023354034140439099065
        cb = 0.08664034996495772158907019

        self.T_c = cb * self.a / (ca * self.b)
        self.p_c = self.T_c * cb / self.b
        self.rho_c = rc / self.b


class VanderWall(EOS):

    def __init__(self, **kwargs):
        self.a = kwargs.get("a")
        self.b = kwargs.get("b")
        self.R = kwargs.get("R")
        self.set_critical_parameters()

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        if value > 0:
            self._a = value

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        if value > 0:
            self._b = value

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if value > 0:
            self._R = value

    def EOS(self, rho, T):
        return (rho * self.R * T) / (1.0 - self.b * rho) - self.a * rho**2

    def dp_drho(self, rho, T):
        a = self.a
        b = self.b
        return T / (1.0 - b * rho) + rho * T * b / (1.0 - b * rho) ** 2 - 2.0 * a * rho

    def dp2_drho2(self, rho, T):
        a = self.a
        b = self.b
        return (
            2.0 * T / (1.0 - b * rho) ** 2 * b
            + 2.0 * rho * T / (1.0 - b * rho) ** 3 * b**2
            - 2.0 * a
        )

    def set_critical_parameters(self):
        self.p_c = self.a / (27.0 * self.b**2)
        self.T_c = 8.0 * self.a / (27.0 * self.b)
        self.rho_c = 1.0 / (3.0 * self.b)
