from src.base import MaxwellConstructor
from src.eos import CarnahanStarling

a = 1.0
b = 4.0
R = 1.0

kwargs = {"a": a, "b": b, "R": R}

eos = CarnahanStarling(**kwargs)

Tc = (a / b) * (0.18727 / 0.4963) * (1 / R)

mc = MaxwellConstructor(eos)
mc.get_sat_density(T_sat=0.75 * Tc, show_output=True)
