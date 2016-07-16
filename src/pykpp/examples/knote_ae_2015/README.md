# Typical Summer Day

This example application is based on the Knote et al. (2015) simulation of a typical summer day with CB05. The example here uses the "Summer" PBL and temperature profile. At present, the example uses the a relatively "clean" initial environment (NOx(t=0) = 1ppb). 

Where possible, the inputs are as described in the publication, but has been adapted for pykpp. The model uses initial conditions, temperature, PBL and deposition rates as described in the publications. The emissions were produced from the NEI using a similar process as described in the paper, but may be somewhat different.

Figures showing the inputs and outputs are distributed with the example and will be remade if you run the example (python summerclean.py). Figure 1 shows the input diurnal profiles of temperature, pressure, NOx emissions, and photolysis. Figure 2 shows the output diurnal profiles of ozone, hydroxyl radical and hydroperoxy radical concentrations. Figure 3 shows the output diurnal profiles of NO, NO2 and the sum.

![Figure 1](physical.pdf)

![Figure 2](chemical_ozone.pdf)

![Figure 3](chemical_nox.pdf)

[1] Knote, C., Tuccella, P., Curci, G., Emmons, L., Orlando, J. J., Madronich, S., Baró, R., Jiménez-Guerrero, P., Luecken, D., Hogrefe, C., Forkel, R., Werhahn, J., Hirtl, M., Pérez, J. L., San José, R., Giordano, L., Brunner, D., Yahya, K., Zhang, Y., Influence of the choice of gas-phase mechanism on predictions of key gaseous pollutants during the AQMEII phase-2 intercomparison, Atmospheric Environment, Volume 115, August 2015, Pages 553-568, ISSN 1352-2310, http://dx.doi.org/10.1016/j.atmosenv.2014.11.066.