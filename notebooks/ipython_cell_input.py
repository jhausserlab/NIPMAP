import random as r
import math as m

# PI approximation using monte-carlo simulation.
def pi_approximator(total_darts, radius):
    inside = 0
    for i in range(0, total_darts):
        x2 = r.uniform(-radius, radius)
        y2 = r.uniform(-radius, radius)
        if m.sqrt(x2**2 + y2**2) < radius:
            inside += 1
    return (float(inside) / total_darts) * 4
