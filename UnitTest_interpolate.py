import numpy as np
import matplotlib.pyplot as plt

def interpolate(z1, z2, z3, z4, x):
    """
    Third Order Interpolation function
    Reference: Section 2.5.7.1 of GSA's "Ionospheric Correction
    Algorithm for Galileo Single Frequency Users"
    """
    if abs(2 * x) < 10 ** -10:
        return z2

    delta = 2 * x - 1
    g1 = z3 + z2
    g2 = z3 - z2
    g3 = z4 + z1
    g4 = (z4 - z1) / 3.0

    a0 = 9 * g1 - g3
    a1 = 9 * g2 - g4
    a2 = g3 - g1
    a3 = g4 - g2

    return 1 / 16.0 * (a0 + a1 * delta + a2 * delta ** 2 + a3 * delta ** 3)

def lagrange(z1,z2,z3,z4, x):

    a0 = -z1/6.0 + 0.5 * z2 -0.5 *z3 + z4/6.0
    a1 = 0.5*z1 - z2 + 0.5*z3
    a2 = -z1/3.0 -0.5*z2 + z3 -z4/6.0
    a3 = z2

    return a3 + a2 * x + a1*x**2 + a0 * x**3

z1 = 2
z2 = 0
z3 = 1
z4 = 4

x = np.linspace(-1,2,100)
y = np.empty(100)
z = np.empty(100)
for i in range(100):
    y[i] = interpolate(z1,z2,z3,z4, x[i])
    z[i] = lagrange(z1,z2,z3,z4, x[i])

plt.plot(x,y)
plt.plot(x,z)
plt.grid()
plt.show()
