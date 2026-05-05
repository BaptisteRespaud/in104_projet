import numpy as np
import matplotlib.pyplot as plt

"""
pi_approximations = np.loadtxt("pi_iterations")

plt.plot(pi_approximations)
plt.xlabel("Iteration")
plt.ylabel("Approximation of π")
plt.savefig("plot.png")  
"""

avg_sqr_mag_L6 = np.loadtxt("iterations/sqr_avg_mag_vsT_L6.txt")
avg_sqr_mag_L8 = np.loadtxt("iterations/sqr_avg_mag_vsT_L8.txt")
avg_sqr_mag_L10 = np.loadtxt("iterations/sqr_avg_mag_vsT_L10.txt")
avg_sqr_mag_L12 = np.loadtxt("iterations/sqr_avg_mag_vsT_L12.txt")

T0 = 1.5
dT = 0.1

temps = T0 + dT * np.arange(len(avg_sqr_mag_L6))

plt.plot(temps, avg_sqr_mag_L6, marker='o', linestyle='-', color='blue', label="L = 6")
plt.plot(temps, avg_sqr_mag_L8, marker='o', linestyle='-', color='orange', label="L = 8")
plt.plot(temps, avg_sqr_mag_L10, marker='o', linestyle='-', color='green', label="L = 10")
plt.plot(temps, avg_sqr_mag_L12, marker='o', linestyle='-', color='red', label="L = 12")

plt.xlabel("Temperature T")
plt.ylabel(r"$\langle m^2 \rangle$")
plt.title("Squared Magnetization vs Temperature")
plt.grid(True)

plt.legend()  
plt.savefig("plot_avg_sqr_mag.png")

