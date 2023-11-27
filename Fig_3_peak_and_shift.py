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
the peak absorption of the quantum-confined Stark effect vs. voltage for
different temperatures and the wavelength of the peak absorption vs. voltage.

Usage:
Run the code to plot Fig. 3 in [1]

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
file_name_peak_TM73C_doped = 'doped_peak_200K.txt'
path_file_name_peak_TM73C_doped = os.path.join(
    p_doped_sample, file_name_peak_TM73C_doped)
file_name_shift_TM73C_doped = 'doped_shift_200K.txt'
path_file_name_shift_TM73C_doped = os.path.join(
    p_doped_sample, file_name_shift_TM73C_doped)
data_peak_TM73C_doped = np.loadtxt(
    path_file_name_peak_TM73C_doped, delimiter=',', dtype=complex)
data_shift_TM73C_doped = np.loadtxt(
    path_file_name_shift_TM73C_doped, delimiter=',', dtype=complex)

# Data at 21 ºC for the p-doped sample is loaded.
file_name_peak_T21C_doped = 'doped_peak_294K.txt'
path_file_name_peak_T21C_doped = os.path.join(
    p_doped_sample, file_name_peak_T21C_doped)
file_name_shift_T21C_doped = 'doped_shift_294K.txt'
path_file_name_shift_T21C_doped = os.path.join(
    p_doped_sample, file_name_shift_T21C_doped)
data_peak_T21C_doped = np.loadtxt(
    path_file_name_peak_T21C_doped, delimiter=',', dtype=complex)
data_shift_T21C_doped = np.loadtxt(
    path_file_name_shift_T21C_doped, delimiter=',', dtype=complex)

# Data at 100 ºC for the p-doped sample is loaded.
file_name_peak_T100C_doped = 'doped_peak_373K.txt'
path_file_name_peak_T100C_doped = os.path.join(
    p_doped_sample, file_name_peak_T100C_doped)
file_name_shift_T100C_doped = 'doped_shift_373K.txt'
path_file_name_shift_T100C_doped = os.path.join(
    p_doped_sample, file_name_shift_T100C_doped)
data_peak_T100C_doped = np.loadtxt(
    path_file_name_peak_T100C_doped, delimiter=',', dtype=complex)
data_shift_T100C_doped = np.loadtxt(
    path_file_name_shift_T100C_doped, delimiter=',', dtype=complex)

# Voltages
V1 = [0, 1, 2, 3, 4, 5]
V2 = [0, 1, 3, 5, 7, 9, 10]
V3 = [0, 1, 3, 5, 7, 10]

# This part plots the peak absorption vs. voltage of the sample with the
# p-doped stack at different temperatures.
plt.figure(0)
plt.plot(V1, -data_peak_TM73C_doped, 'bo', label='PD @ -73$^\circ$C')
plt.plot(V2, -data_peak_T21C_doped, 'ko', label='PD @ 21$^\circ$C')
plt.plot(V3, -data_peak_T100C_doped, 'ro', label='PD @ 100$^\circ$C')
plt.xlabel('Voltage [V]')
plt.ylabel('Peak Absorption [cm$^-$$^1$]')
plt.rcParams.update({'font.size': 17})

# This part plots the wavelength of the peak absorption vs. voltage of the
# sample with the p-doped stack at different temperatures.
plt.figure(1)
plt.plot(V1, CONVERSION_E_WAVELENGTH/data_shift_TM73C_doped,
         'bo', label='PD @ -73$^\circ$C')
plt.plot(V2, CONVERSION_E_WAVELENGTH/data_shift_T21C_doped,
         'ko', label='PD @ 21$^\circ$C')
plt.plot(V3, CONVERSION_E_WAVELENGTH/data_shift_T100C_doped,
         'ro', label='PD @ 100$^\circ$C')
plt.xlabel('Voltage [V]')
plt.ylabel('Peak Absorption Energy [eV]')
plt.legend(prop={'size': 8})

# This part extracts the data of the non-intentionally-doped sample.
# Note: An M is used in the name of the variables to represent the minus sign
# in temperatures below 0 °C. E.g., -73 °C is M73.

# Data at -73 ºC for the non-intentionally-doped sample is loaded.
file_name_peak_TM73C_NID = 'undoped_peak_200K.txt'
path_file_name_peak_TM73C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_peak_TM73C_NID)
file_name_shift_TM73C_NID = 'undoped_shift_200K.txt'
path_file_name_shift_TM73C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_shift_TM73C_NID)
data_peak_TM73C_NID = np.loadtxt(
    path_file_name_peak_TM73C_NID, delimiter=',', dtype=complex)
data_shift_TM73C_NID = np.loadtxt(
    path_file_name_shift_TM73C_NID, delimiter=',', dtype=complex)

# Data at 21 ºC for the non-intentionally-doped sample is loaded.
file_name_peak_T21C_NID = 'undoped_peak_294K.txt'
path_file_name_peak_T21C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_peak_T21C_NID)
file_name_shift_T21C_NID = 'undoped_shift_294K.txt'
path_file_name_shift_T21C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_shift_T21C_NID)
data_peak_T21C_NID = np.loadtxt(
    path_file_name_peak_T21C_NID, delimiter=',', dtype=complex)
data_shift_T21C_NID = np.loadtxt(
    path_file_name_shift_T21C_NID, delimiter=',', dtype=complex)

# Data at 100 ºC for the non-intentionally-doped sample is loaded.
file_name_peak_T100C_NID = 'undoped_peak_373K.txt'
path_file_name_peak_T100C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_peak_T100C_NID)
file_name_shift_T100C_NID = 'undoped_shift_373K.txt'
path_file_name_shift_T100C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_shift_T100C_NID)
data_peak_T100C_NID = np.loadtxt(
    path_file_name_peak_T100C_NID, delimiter=',', dtype=complex)
data_shift_T100C_NID = np.loadtxt(
    path_file_name_shift_T100C_NID, delimiter=',', dtype=complex)

# Voltages
V1 = [0, 1, 2, 3, 4, 5]
V2 = [0, 1, 3, 5, 7, 9, 10]

# This part plots the peak absorption vs. voltage of the sample with the
# non-intentionally-doped stack at different temperatures.
plt.figure(0)
plt.plot(V1, -data_peak_TM73C_NID, 'bs', label='NID @ -73$^\circ$C')
plt.plot(V1, -data_peak_T21C_NID, 'ks', label='NID @ 21$^\circ$C')
plt.plot(V1, -data_peak_T100C_NID, 'rs', label='NID @ 100$^\circ$C')
plt.xlabel('Voltage [V]')
plt.ylabel('Peak Absorption [cm$^-$$^1$]')
plt.legend(prop={'size': 9})
plt.annotate('a)',
             xy=(0, 60))
plt.savefig('Fig_3a.png', dpi=1080, bbox_inches='tight')

# This part plots the wavelength of the peak absorption vs. voltage of the
# sample with the non-intentionally-doped stack at different temperatures.
plt.figure(1)
plt.plot(V1, CONVERSION_E_WAVELENGTH/data_shift_TM73C_NID,
         'bs', label='NID @ -73$^\circ$C')
plt.plot(V1, CONVERSION_E_WAVELENGTH/data_shift_T21C_NID,
         'ks', label='NID @ 21$^\circ$C')
plt.plot(V1, CONVERSION_E_WAVELENGTH/data_shift_T100C_NID,
         'rs', label='NID @ 100$^\circ$C')
plt.annotate('b)',
             xy=(0, 1290))
plt.xlabel('Voltage [V]')
plt.ylabel('Wavelength [nm]')
plt.legend(prop={'size': 9}, loc='upper right', bbox_to_anchor=(1, 0.55))
plt.savefig('Fig_3b.png', dpi=1080, bbox_inches='tight')

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
