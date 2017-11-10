"""
This is a handcrafted plot of 14N/15N vs E_upper for 

"""

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

Line = namedtuple("Line", ['Ju', 'Freq', 'Eu'])

hcn_lines = [
    Line(1 , 86.33992140, 4.14369),
    Line(2 , 172.67785120, 12.43099),
    Line(3 , 259.01179760, 24.86166),
    Line(4 , 345.33976930, 41.43546),
    Line(5 , 431.65977480, 62.15200),
    Line(6 , 517.96982100, 87.01076),
    Line(7 , 604.26791400, 116.01122),
    Line(8 , 690.55207900, 149.15270),
    Line(9 , 776.82031970, 186.43453),
    Line(10, 863.07063400, 227.85560),
    Line(11, 949.30103900, 273.41526),
    Line(12, 1035.50955010, 323.11210),
    Line(13, 1121.69417700, 376.94531)
]

observed_Jus = [1,3,4,6,7,8,9]

observed_hcn_lines = [x for x in hcn_lines if x.Ju in observed_Jus]
E_uppers = np.array([x.Eu for x in observed_hcn_lines])

Ratio_14_15 = namedtuple("Ratio", ['Ju', 'value', 'error'])

fourteen_fifteen_nitrogen_ratios = [
    Ratio_14_15 (1 , 390 , 70),
    Ratio_14_15 (3 , 240 , 60),
    Ratio_14_15 (4 , 410 , 100),
    Ratio_14_15 (6 , 120 , 10),
    Ratio_14_15 (7 , 150 , 40),
    Ratio_14_15 (8 , 150 , 90),
    Ratio_14_15 (9 , 210 , 70)
]

ratios = np.array([x.value for x in fourteen_fifteen_nitrogen_ratios])
errors = np.array([x.error for x in fourteen_fifteen_nitrogen_ratios])

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111)

ax.errorbar(E_uppers, ratios, yerr=errors, fmt='ro', ms=4)
ax.set_xlabel("Upper state energy (K)", fontsize=16)
ax.set_ylabel(r"$\frac{^{14}\rm{N}}{^{15}\rm{N}}$", fontsize=24, rotation='horizontal', labelpad=20)
ax.set_title(r"$^{14}$N/$^{15}$N in IRAS 16293 from single-dish observations")

ax.plot([-100,200], [388,388], ':', color='C2', lw=1, scalex=False)
ax.text(55, 390, "ISM standard: 388", color='C2', fontsize=12)
ax.plot([-100,200], [270,270], ':', color='C1', lw=1, scalex=False)
ax.text(55, 272, "Solar system standard: 270", color='C1', fontsize=12)

plt.show()

fig.savefig("COST_plot_1415ratio_vs_Eu.pdf", bbox_inches='tight')
fig.savefig("COST_plot_1415ratio_vs_Eu.png", bbox_inches='tight')