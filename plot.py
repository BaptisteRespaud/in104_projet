import numpy as np
import matplotlib.pyplot as plt

pi_approximations = np.loadtxt("pi_iterations")

plt.plot(pi_approximations)
plt.xlabel("Iteration")
plt.ylabel("Approximation of π")
plt.savefig("plot.png")   
