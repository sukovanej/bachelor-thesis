from numpy import sqrt, cos, sin, linspace, exp
import scipy.integrate as integrate
import matplotlib.pyplot as plt

hbar, c, m = 1, 1, 1
TIME = 100
eps = 1e-6

def function(x, t):
    def real(p):
        k = - 1 / hbar * c * (sqrt(p ** 2 + (m * c)**2)) * t
        l = 1 / hbar * p * x
        return cos(k) * cos(l) - sin(k) * sin(l) * exp(- eps * p**2)

    def imaginary(p):
        k = - 1 / hbar * c * (sqrt(p ** 2 + (m * c)**2)) * t
        l = 1 / hbar * p * x
        return sin(k) * cos(l) + cos(k) * sin(l) * exp(- eps * p**2)

    return real, imaginary

x_list = linspace(0, 3 * TIME, 500)
y_list = []

for x in x_list:
    real, imaginary = function(x, TIME)

    r = integrate.quad(real, -c/10, c/10)
    i = integrate.quad(imaginary, -c/10, c/10)
    y_list.append(sqrt(r[0] ** 2 + i[0] ** 2) / 2)

plt.xlabel("x-coordinate for $t = {}c$".format(TIME))
plt.ylabel("propagation probability")
plt.plot(x_list, y_list)
plt.plot([TIME, TIME], [0, max(y_list) * 1.2])
plt.gca().set_ylim([0, max(y_list) * 1.2])
#plt.savefig("free_relativistic_particle.png")
plt.show()
