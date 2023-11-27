"""
Title:
Measurement of the quantum-confined Stark effect in InAs/In(Ga)As quantum dots
with p-doped quantum dot barriers.

Authors:
- Dr. Nicolás Abadía, Institute for Compound Semiconductors and Cardiff
University.
- Gina Arnau Torner, Princeton University.

Journal:
- Journal: Optics Express
- Date: 2022
- Volume: 30
- Issue: 11
- Page range: 17730-17738
- DOI: 10.1364/OE.455491 - https://doi.org/10.1364/OE.455491

Abstract:
The quantum-confined Stark effect in InAs/In(Ga)As quantum dots (QDs) using
non-intentionally doped and p-doped QD barriers was investigated to compare
their performance for use in optical modulators. The measurements indicate
that the doped QD barriers lead to a better figure of merit (FoM), defined as
the ratio of the change in absorption Δα for a reverse bias voltage swing to
the loss at 1 V α(1 V), FoM=Δα/α (1 V). The improved performance is due to the
absence of the ground-state absorption peak and an additional component to the
Stark shift. Measurements indicate that p-doping the QD barriers can lead to
more than a 3x increase in FoM modulator performance between temperatures of
−73 °C to 100 °C when compared with the stack with NID QD barriers.

Description:
This Python script corresponds to the implementation of the research described
in the paper titled "Measurement of the quantum-confined Stark effect in
InAs/In(Ga)As quantum dots with p-doped quantum dot barriers." This code plots
the Figure of Merit (FoM) vs. wavelength (λ) for a 5V-swing at different
temperatures of the doped and non-intentionally-doped stacks in [1].

Usage:
Run the code to plot Fig. 5(a) in [1]

License:
Published by Optica Publishing Group under the terms of the Creative Commons
Attribution 4.0 License. Further distribution of this work must maintain
attribution to the author(s) and the published article's title, journal
citation, and DOI.

Citation:
[1] Joe Mahoney, Mingchu Tang, Huiyun Liu, and Nicolás Abadía, "Measurement of
the quantum-confined Stark effect in InAs/In(Ga)As quantum dots with p-doped
quantum dot barriers," Opt. Express 30, 17730-17738 (2022) -
https://doi.org/10.1364/OE.455491

Acknowledgments:
We acknowledge Gina Arnau Torner for polishing the code.

Contact:
Email: abadian@cardiff.ac.uk

Note:
Thank you for using our Python code! If you encounter any typos, bugs, or have
suggestions for improvements, we encourage you to reach out to us. Your
feedback is valuable in enhancing the quality of our code.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# The data is stored in eV, and the plots are in nm. This constant is used to
# convert from energy to wavelength using the formula λ = (h·c)/E [nm].
CONVERSION_E_WAVELENGTH = 1e9 * 6.63e-34 * 3e8 / (1.6e-19)

# This routine sets a relative path for Unix-like or Windows systems.


def dir_create(Dir):
    # '/' on Unix-like systems, '\' on Windows systems.
    D = Dir.replace(os.sep, '/')
    return D


# Directory containing the data of the p-doped sample.
p_doped_sample = dir_create('./data/doped')
# Directory containing the data of the non-intentionally-doped sample.
non_intentionally_doped_sample = dir_create('./data/undoped')

# This part extracts the data of the p-doped sample.
# Note: An M is used in the name of the variables to represent the minus sign
# in temperatures below 0 °C. E.g., -73 °C is M73.

# Data at -73 ºC for the p-doped sample is loaded.
file_name_1V_TM73C_doped = \
    'UCL_355_R4_1681_B1_S2S3_MICROSTAT_N200K_1V/loss.txt'
path_file_name_1V_TM73C_doped = os.path.join(
    p_doped_sample, file_name_1V_TM73C_doped)
file_name_5V_TM73C_doped = \
    'UCL_355_R4_1681_B1_S2S3_MICROSTAT_N200K_5V_take_2/loss.txt'
path_file_name_5V_TM73C_doped = os.path.join(
    p_doped_sample, file_name_5V_TM73C_doped)
data_wavelength_vs_loss_1V_TM73C_doped = np.loadtxt(
    path_file_name_1V_TM73C_doped, delimiter=',', dtype=complex)
data_wavelength_vs_loss_5V_TM73C_doped = np.loadtxt(
    path_file_name_5V_TM73C_doped, delimiter=',', dtype=complex)
FoM_TM73C_doped = (data_wavelength_vs_loss_5V_TM73C_doped[:300, 1] -
                   data_wavelength_vs_loss_1V_TM73C_doped[:300, 1]-0.8) / \
                  (data_wavelength_vs_loss_1V_TM73C_doped[:300, 1]+0.8)
max_FoM_TM73C_doped = np.max(FoM_TM73C_doped)
max_E_1V_TM73C_doped = np.where(FoM_TM73C_doped == np.max(FoM_TM73C_doped))[0]
max_wavelength_TM73C_doped = CONVERSION_E_WAVELENGTH / \
    data_wavelength_vs_loss_5V_TM73C_doped[max_E_1V_TM73C_doped, 0]

# Data at 21 ºC for the p-doped sample is loaded. The FoM is calculated.
file_name_1V_T21C_doped = 'UCL_1681_355_R4_B1_S2S3_1V/loss.txt'
path_file_name_1V_T21C_doped = os.path.join(
    p_doped_sample, file_name_1V_T21C_doped)
file_name_5V_T21C_doped = \
    'UCL_355_1681_B1_S2S3_21_22nd_Sep_5V_21_Deg/loss.txt'
path_file_name_5V_T21C_doped = os.path.join(
    p_doped_sample, file_name_5V_T21C_doped)
data_wavelength_vs_loss_1V_T21C_doped = np.loadtxt(
    path_file_name_1V_T21C_doped, delimiter=',', dtype=complex)
data_wavelength_vs_loss_5V_T21C_doped = np.loadtxt(
    path_file_name_5V_T21C_doped, delimiter=',', dtype=complex)
FoM_T21C_doped = (-data_wavelength_vs_loss_1V_T21C_doped[:290, 1] +
                  data_wavelength_vs_loss_5V_T21C_doped[:290, 1]) / \
    (data_wavelength_vs_loss_1V_T21C_doped[:290, 1])
max_FoM_T21C_doped = np.max(FoM_T21C_doped)
max_E_1V_T21C_doped = np.where(FoM_T21C_doped == np.max(FoM_T21C_doped))[0]
max_wavelength_T21C_doped = CONVERSION_E_WAVELENGTH / \
    data_wavelength_vs_loss_5V_T21C_doped[max_E_1V_T21C_doped, 0]

# Data at 50 ºC for the p-doped sample is loaded. The FoM is calculated.
file_name_1V_T50C_doped = 'UCL_1681_355_R4_B1_S2S3_50_DEG_1V/loss.txt'
path_file_name_1V_T50C_doped = os.path.join(
    p_doped_sample, file_name_1V_T50C_doped)
file_name_5V_T50C_doped = 'UCL_1681_355_R4_B1_S2S3_50_DEG_5V/loss.txt'
path_file_name_5V_T50C_doped = os.path.join(
    p_doped_sample, file_name_5V_T50C_doped)
data_wavelength_vs_loss_1V_T50C_doped = np.loadtxt(
    path_file_name_1V_T50C_doped, delimiter=',', dtype=complex)
data_wavelength_vs_loss_5V_T50C_doped = np.loadtxt(
    path_file_name_5V_T50C_doped, delimiter=',', dtype=complex)
FoM_T50C_doped = (-data_wavelength_vs_loss_1V_T50C_doped[:300, 1] +
                  data_wavelength_vs_loss_5V_T50C_doped[:300, 1]) / \
    (data_wavelength_vs_loss_1V_T50C_doped[:300, 1])
max_FoM_T50C_doped = np.max(FoM_T50C_doped)
max_E_1V_T21C = np.where(FoM_T50C_doped == np.max(FoM_T50C_doped))[0]
max_wavelength_T50C_doped = CONVERSION_E_WAVELENGTH / \
    data_wavelength_vs_loss_5V_T50C_doped[max_E_1V_T21C, 0]

# Data at 75 ºC for the p-doped sample is loaded. The FoM is calculated.
file_name_1V_T75C_doped = 'UCL_1681_355_R4_B1_S2S3_75_DEG_1V/loss.txt'
path_file_name_1V_T75C_doped = os.path.join(
    p_doped_sample, file_name_1V_T75C_doped)
file_name_5V_T75C_doped = 'UCL_1681_355_R4_B1_S2S3_75_DEG_5V/loss.txt'
path_file_name_5V_T75C_doped = os.path.join(
    p_doped_sample, file_name_5V_T75C_doped)
data_wavelength_vs_loss_1V_T75C_doped = np.loadtxt(
    path_file_name_1V_T75C_doped, delimiter=',', dtype=complex)
data_wavelength_vs_loss_5V_T75C_doped = np.loadtxt(
    path_file_name_5V_T75C_doped, delimiter=',', dtype=complex)
FoM_T75C_doped = (-data_wavelength_vs_loss_1V_T75C_doped[:290, 1] +
                  data_wavelength_vs_loss_5V_T75C_doped[:290, 1]
                  + 0.2) / (data_wavelength_vs_loss_1V_T75C_doped[:290, 1]-0.2)
max_FoM_T75C_doped = np.max(FoM_T75C_doped)
max_E_1V_T75C_doped = np.where(FoM_T75C_doped == np.max(FoM_T75C_doped))[0]
max_wavelength_T75C_doped = CONVERSION_E_WAVELENGTH / \
    data_wavelength_vs_loss_5V_T75C_doped[max_E_1V_T75C_doped, 0]

# Data at 100 ºC for the p-doped sample is loaded. The FoM is calculated.
file_name_1V_T100C_doped = \
    'UCL_355_1681_B1_S2S3_21_22nd_Sep_1V_100_Deg/loss.txt'
path_file_name_1V_T100C_doped = os.path.join(
    p_doped_sample, file_name_1V_T100C_doped)
file_name_5V_T100C_doped = \
    'UCL_355_B1_R4_1681_00_DEG_5V/loss.txt'
path_file_name_5V_T100C_doped = os.path.join(
    p_doped_sample, file_name_5V_T100C_doped)
data_wavelength_vs_loss_1V_T100C_doped = np.loadtxt(
    path_file_name_1V_T100C_doped, delimiter=',', dtype=complex)
data_wavelength_vs_loss_5V_T100C_doped = np.loadtxt(
    path_file_name_5V_T100C_doped, delimiter=',', dtype=complex)
FoM_T100C_doped = (-data_wavelength_vs_loss_1V_T100C_doped[:300, 1]-0.4 +
                   data_wavelength_vs_loss_5V_T100C_doped[:300, 1]) / \
    (data_wavelength_vs_loss_1V_T100C_doped[:300, 1]+0.4)
max_FoM_T100C_doped = np.max(FoM_T100C_doped)
max_E_1V_T100C_doped = np.where(FoM_T100C_doped == np.max(FoM_T100C_doped))[0]
max_wavelength_T100C_doped = CONVERSION_E_WAVELENGTH / \
    data_wavelength_vs_loss_5V_T100C_doped[max_E_1V_T100C_doped, 0]

# This part plots the FoM of the sample with the p-doped stack at different
# temperatures.
plt.figure(1)

plt.plot(CONVERSION_E_WAVELENGTH /
         data_wavelength_vs_loss_5V_TM73C_doped[:300, 0],
         FoM_TM73C_doped, 'b', label='PD @ -73$^\circ$C', color='black')
plt.plot(CONVERSION_E_WAVELENGTH /
         data_wavelength_vs_loss_5V_T21C_doped[:290, 0],
         FoM_T21C_doped, 'b', label='PD @ 21$^\circ$C', color='blue')
plt.plot(CONVERSION_E_WAVELENGTH /
         data_wavelength_vs_loss_5V_T50C_doped[:300, 0],
         FoM_T50C_doped, 'g', label='PD @ 50$^\circ$C', color='g')
plt.plot(CONVERSION_E_WAVELENGTH /
         data_wavelength_vs_loss_5V_T75C_doped[:290, 0], FoM_T75C_doped,
         'y', label='PD @ 75$^\circ$C')
plt.plot(CONVERSION_E_WAVELENGTH /
         data_wavelength_vs_loss_5V_T100C_doped[:300, 0], FoM_T100C_doped,
         'r', label='PD @ 100$^\circ$C')

# This part extracts the data of the non-intentionally-doped sample.
# Note: An M is used in the name of the variables to represent the minus sign
# in temperatures below 0 °C. E.g., -73 °C is M73.

# Data at -73 °C for the non-intentionally-doped sample is loaded. The
# extinction ratio (ER) and the insertion loss (IL) are calculated.
E_TM73C = 'D1_5V_FOM/E_200.txt'
path_E_TM73C = os.path.join(non_intentionally_doped_sample, E_TM73C)
data_E_TM73C = np.loadtxt(path_E_TM73C, delimiter=',', dtype=complex)
ER_TM73C = 'D1_5V_FOM/ER_200.txt'
path_ER_TM73C = os.path.join(non_intentionally_doped_sample, ER_TM73C)
data_ER_TM73C = np.loadtxt(path_ER_TM73C, delimiter=',', dtype=complex)
IL_TM73C = 'D1_5V_FOM/IL_200.txt'
path_IL_TM73C = os.path.join(non_intentionally_doped_sample, IL_TM73C)
data_IL_TM73C = np.loadtxt(path_IL_TM73C, delimiter=',', dtype=complex)

# Data at 21 °C for the non-intentionally-doped sample is loaded. The
# extinction ratio (ER) and the insertion loss (IL) are calculated.
E_T21C = 'D1_5V_FOM/E_21.txt'
path_E_T21C = os.path.join(non_intentionally_doped_sample, E_T21C)
data_E_T21C = np.loadtxt(path_E_T21C, delimiter=',', dtype=complex)
ER_T21C = 'D1_5V_FOM/ER_21.txt'
path_ER_T21C = os.path.join(non_intentionally_doped_sample, ER_T21C)
data_ER_T21C = np.loadtxt(path_ER_T21C, delimiter=',', dtype=complex)
IL_T21C = 'D1_5V_FOM/IL_21.txt'
path_IL_T21C = os.path.join(non_intentionally_doped_sample, IL_T21C)
data_IL_T21C = np.loadtxt(path_IL_T21C, delimiter=',', dtype=complex)

# Data at 50 °C for the non-intentionally-doped sample is loaded. The
# extinction ratio (ER) and the insertion loss (IL) are calculated.
E_T50C = 'D1_5V_FOM/E_50.txt'
path_E_T50C = os.path.join(non_intentionally_doped_sample, E_T50C)
data_E_T50C = np.loadtxt(path_E_T50C, delimiter=',', dtype=complex)
ER_T50C = 'D1_5V_FOM/ER_50.txt'
path_ER_T50C = os.path.join(non_intentionally_doped_sample, ER_T50C)
data_ER_T50C = np.loadtxt(path_ER_T50C, delimiter=',', dtype=complex)
IL_T50C = 'D1_5V_FOM/IL_50.txt'
path_IL_T50C = os.path.join(non_intentionally_doped_sample, IL_T50C)
data_IL_T50C = np.loadtxt(path_IL_T50C, delimiter=',', dtype=complex)

# Data at 75 °C for the non-intentionally-doped sample is loaded. The
# extinction ratio (ER) and the insertion loss (IL) are calculated.
E_T75C = 'D1_5V_FOM/E_75.txt'
path_E_T75C = os.path.join(non_intentionally_doped_sample, E_T75C)
data_E_T75C = np.loadtxt(path_E_T75C, delimiter=',', dtype=complex)
ER_T75C = 'D1_5V_FOM/ER_75.txt'
path_ER_T75C = os.path.join(non_intentionally_doped_sample, ER_T75C)
data_ER_T75C = np.loadtxt(path_ER_T75C, delimiter=',', dtype=complex)
IL_T75C = 'D1_5V_FOM/IL_75.txt'
path_IL_T75C = os.path.join(non_intentionally_doped_sample, IL_T75C)
data_IL_T75C = np.loadtxt(path_IL_T75C, delimiter=',', dtype=complex)

# Data at 100 °C for the non-intentionally-doped sample is loaded. The
# extinction ratio (ER) and the insertion loss (IL) are calculated.
E_T100C = 'D1_5V_FOM/E_100.txt'
path_E_T100C = os.path.join(non_intentionally_doped_sample, E_T100C)
data_E_T100C = np.loadtxt(path_E_T100C, delimiter=',', dtype=complex)
ER_T100C = 'D1_5V_FOM/ER_100.txt'
path_ER_T100C = os.path.join(non_intentionally_doped_sample, ER_T100C)
data_ER_T100C = np.loadtxt(path_ER_T100C, delimiter=',', dtype=complex)
IL_T100C = 'D1_5V_FOM/IL_100.txt'
path_IL_T100C = os.path.join(non_intentionally_doped_sample, IL_T100C)
data_IL_T100C = np.loadtxt(path_IL_T100C, delimiter=',', dtype=complex)


# This part plots the FoM of the sample with the non-intentionally-doped stack
# at different temperatures.
plt.plot(CONVERSION_E_WAVELENGTH/data_E_TM73C[:255], data_ER_TM73C[:255] /
         data_IL_TM73C[:255], 'k--', label='NID @ -73$^\circ$C')
plt.plot(CONVERSION_E_WAVELENGTH/data_E_T21C[:255], data_ER_T21C[:255] /
         data_IL_T21C[:255], 'b--', label='NID @ 21$^\circ$C')
plt.plot(CONVERSION_E_WAVELENGTH/data_E_T50C[:235], data_ER_T50C[:235] /
         data_IL_T50C[:235], 'g--', label='NID @ 50$^\circ$C')
plt.plot(CONVERSION_E_WAVELENGTH/data_E_T75C[:300], data_ER_T75C[:300] /
         data_IL_T75C[:300], 'y--', label='NID @ 75$^\circ$C')
plt.plot(CONVERSION_E_WAVELENGTH/data_E_T100C[:330], data_ER_T100C[:330] /
         data_IL_T100C[:330], 'r--', label='NID @ 100$^\circ$C')
plt.xlim(1250, 1420)
plt.ylim(-0.06, 9)
plt.xlabel('$\lambda$ [nm]')
plt.ylabel(r'FoM = $\Delta$$\alpha\cdot\alpha^{-1}(1V)$')
plt.rcParams.update({'font.size': 17})
plt.legend(prop={'size': 6.5}, loc=1)
# This variable is used to mark a black dotted line around 2.4 dBmm^{-1}
reference_FoM = np.ones(len(data_wavelength_vs_loss_5V_T100C_doped[:, 0]))
reference_FoM *= 2.4
plt.plot(CONVERSION_E_WAVELENGTH /
         data_wavelength_vs_loss_5V_T100C_doped[:, 0], reference_FoM, 'k--')
plt.annotate('a) 4V-swing',
             xy=(1255, 8.2))
plt.annotate('FoM = 2.4',
             xy=(1350, 2.5))
plt.savefig('Fig_5a.png', dpi=1080, bbox_inches='tight')

"""
Note:
If you find this code useful for your work and choose to incorporate it, we
kindly request that you cite our research paper:

Citation:
[1] Joe Mahoney, Mingchu Tang, Huiyun Liu, and Nicolás Abadía, "Measurement of
the quantum-confined Stark effect in InAs/In(Ga)As quantum dots with p-doped
quantum dot barriers," Opt. Express 30, 17730-17738 (2022) -
https://doi.org/10.1364/OE.455491
"""
