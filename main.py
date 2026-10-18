"""
langmuir linearization project
noah lovett
started oct. 18 2026
"""

import numpy as np
import matplotlib.pyplot as plt

methane = open("methane.txt")
words = methane.read().split()
values = []
for word in words:
    try:
        values.append(float(word))
    except ValueError:
        pass
values = values[2:]
methane.close()
pressures = np.array([]) #bar
concentrations = np.array([]) #mmol/g
for i in range(len(values)):
    if i % 2 == 0:
        pressures = np.append(pressures, values[i])
    else:
        concentrations = np.append(concentrations, values[i])

linearized_pressures = pressures**(-1) # 1/bar
linearized_concentrations = concentrations**(-1) # g /mmol
sum_lp = sum(linearized_pressures)
sum_lc = sum(linearized_concentrations)
n = len(linearized_pressures)

#linear regression plot
m = (n*sum(linearized_concentrations*linearized_pressures) - sum_lc*sum_lp) / (n*sum(linearized_pressures**2) - sum_lp**2)
b = np.mean(linearized_concentrations) - m*(np.mean(linearized_pressures))
x = np.array([])
for i in range(120):
    x = np.append(x, i)

#m = 1/k*Nm; b = 1/Nm
Nm = b**-1
print("Nm = " + str(Nm))
k = (m*Nm)**-1
print("k = " + str(k))

fig, ax = plt.subplots()
ax.plot(x, m*x+b, color='red', label=('Best fit: y = ' + str(round(m, 5)) + "x + " + str(round(b, 5))))
ax.scatter(linearized_pressures, linearized_concentrations, color='blue', marker='o', label='Data')
ax.set_xlabel('1/Pressure [1/bar]')
ax.set_ylabel('1/Concentration [g/mmol]')
ax.legend()

x = x[:5] # pressures up to 5 bar
#langmuir adsorption isotherm:
N = Nm*(k*x)/(1+k*x)
a = 1
N1 = Nm*(k*(1+k*a) - (k**2)*a) / (1+k*a)**2
N2 = -2*Nm*(k**2)*(1+k*a)**-3
N3 = 6*Nm*(k**3)*(1+k*a)**-4
N4 = -24*Nm*(k**4)*(1+k*a)**-5

taylor_series = Nm*(k*a)/(1+k*a) + N1*(x-a) + N2*((x-a)**2)/2 + N3*((x-a)**3)/6 + N4*((x-a)**4)/24

fig, ax2 = plt.subplots()
ax2.plot(x, N, color='orange', label='Langmuir Isotherm Adsorption (T = 293K)')
ax2.scatter(pressures[:10], concentrations[:10], color='blue', marker='o', label='Data')
ax2.legend()
ax2.set_xlabel('Pressure [bar]')
ax2.set_ylabel('Concentration [mmol/g]')

fig, ax3 = plt.subplots()
ax3.scatter(pressures[:10], concentrations[:10], color='blue', marker='o', label='Data')
ax3.plot(x, taylor_series, color='green', label='Taylor Series')
ax3.legend()
ax3.set_xlabel('Pressure [bar]')
ax3.set_ylabel('Concentration [mmol/g]')

plt.show()