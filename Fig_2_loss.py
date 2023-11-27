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
the absorption (α) and change in absorption (Δα) for a 9V-swing at different
temperatures of the doped and non-intentionally-doped stacks in [1].

Usage:
Run the code to plot Fig. 2 in [1]

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
    D = Dir.replace(os.sep, "/")
    return D


# Directory containing the data of the p-doped sample.
p_doped_sample = dir_create("./data/doped")
# Directory containing the data of the non-intentionally-doped sample.
non_intentionally_doped_sample = dir_create("./data/undoped")

# This part extracts the data of the p-doped sample.
# Note: An M is used in the name of the variables to represent the minus sign
# in temperatures below 0 °C. E.g., -73 °C is M73.

# Data at -73 ºC for the p-doped sample is loaded.
file_name_1V_TM73C_doped = \
    "UCL_355_R4_1681_B1_S2S3_MICROSTAT_N200K_1V/loss.txt"
path_file_name_1V_TM73C_doped = os.path.join(
    p_doped_sample, file_name_1V_TM73C_doped)
file_name_10V_TM73C_doped = \
    "UCL_355_R4_1681_B1_S2S3_MICROSTAT_N200K_5V_take_2/loss.txt"
path_file_name_10V_TM73C_doped = os.path.join(
    p_doped_sample, file_name_10V_TM73C_doped)
data_wavelength_vs_loss_1V_TM73C_doped = np.loadtxt(
    path_file_name_1V_TM73C_doped, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_TM73C_doped = np.loadtxt(
    path_file_name_10V_TM73C_doped, delimiter=",", dtype=complex
)

# Data at 21 ºC for the p-doped sample is loaded.
file_name_1V_T21C_doped = "UCL_1681_355_R4_B1_S2S3_1V/loss.txt"
path_file_name_1V_T21C_doped = os.path.join(
    p_doped_sample, file_name_1V_T21C_doped)
file_name_10V_T21C_doped = "UCL_1681_355_R4_B1_S2S3_10V/loss.txt"
path_file_name_10V_T21C_doped = os.path.join(
    p_doped_sample, file_name_10V_T21C_doped)
data_wavelength_vs_loss_1V_T21C_doped = np.loadtxt(
    path_file_name_1V_T21C_doped, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T21C_doped = np.loadtxt(
    path_file_name_10V_T21C_doped, delimiter=",", dtype=complex
)

# Data at 50 ºC for the p-doped sample is loaded.
file_name_1V_T50C_doped = "UCL_1681_355_R4_B1_S2S3_50_DEG_1V/loss.txt"
path_file_name_1V_T50C_doped = os.path.join(
    p_doped_sample, file_name_1V_T50C_doped)
file_name_10V_T50C_doped = "UCL_1681_355_R4_B1_S2S3_50_DEG_10V/loss.txt"
path_file_name_10V_T50C_doped = os.path.join(
    p_doped_sample, file_name_10V_T50C_doped)
data_wavelength_vs_loss_1V_T50C_doped = np.loadtxt(
    path_file_name_1V_T50C_doped, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T50C_doped = np.loadtxt(
    path_file_name_10V_T50C_doped, delimiter=",", dtype=complex
)

# Data at 75 ºC for the p-doped sample is loaded.
file_name_1V_T75C_doped = "UCL_1681_355_R4_B1_S2S3_75_DEG_1V/loss.txt"
path_file_name_1V_T75C_doped = os.path.join(
    p_doped_sample, file_name_1V_T75C_doped)
file_name_10V_T75C_doped = "UCL_1681_355_R4_B1_S2S3_75_DEG_10V/loss.txt"
path_file_name_10V_T75C_doped = os.path.join(
    p_doped_sample, file_name_10V_T75C_doped)
data_wavelength_vs_loss_1V_T75C_doped = np.loadtxt(
    path_file_name_1V_T75C_doped, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T75C_doped = np.loadtxt(
    path_file_name_10V_T75C_doped, delimiter=",", dtype=complex
)

# Data at 100 ºC for the p-doped sample is loaded.
file_name_1V_T100C_doped = \
    "UCL_355_1681_B1_S2S3_21_22nd_Sep_1V_100_Deg/loss.txt"
path_file_name_1V_T100C_doped = os.path.join(
    p_doped_sample, file_name_1V_T100C_doped)
file_name_10V_T100C_doped = \
    "UCL_355_1681_B1_S2S3_21_22nd_Sep_10V_100_Deg/loss.txt"
path_file_name_10V_T100C_doped = os.path.join(
    p_doped_sample, file_name_10V_T100C_doped)
data_wavelength_vs_loss_1V_T100C_doped = np.loadtxt(
    path_file_name_1V_T100C_doped, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T100C_doped = np.loadtxt(
    path_file_name_10V_T100C_doped, delimiter=",", dtype=complex
)

# This part plots the absorption (α) and change in absorption (Δα) for a
# 9V-swing of the sample with the p-doped stack at different temperatures.
plt.figure(1)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T21C_doped[:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T21C_doped[:, 1],
    "b",
    label=r"$\alpha(1V)$ @ 21$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T21C_doped[:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T21C_doped[:, 1]
        - data_wavelength_vs_loss_1V_T21C_doped[:, 1]
    ),
    "b--",
    label=r"$\Delta\alpha(9V-swing)$ @ 21$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T50C_doped[:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T50C_doped[:, 1],
    "g",
    label=r"$\alpha(1V)$ @ 50$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T50C_doped[:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T50C_doped[:, 1]
        - data_wavelength_vs_loss_1V_T50C_doped[:, 1]
    ),
    "g--",
    label=r"$\Delta\alpha(9V-swing)$ @ 50$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T75C_doped[:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T75C_doped[:, 1],
    "y",
    label=r"$\alpha(1V)$ @ 75$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T75C_doped[:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T75C_doped[:, 1]
        - data_wavelength_vs_loss_1V_T75C_doped[:, 1]
    ),
    "y--",
    label=r"$\Delta\alpha(9V-swing)$ @ 75$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T100C_doped[:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T100C_doped[:, 1],
    "r",
    label=r"$\alpha(1V)$ @ 100$^\circ$C",
)
plt.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T100C_doped[:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T100C_doped[:, 1]
        - data_wavelength_vs_loss_1V_T100C_doped[:, 1]
    ),
    "r--",
    label=r"$\Delta\alpha(9V-swing)$ @ 100$^\circ$C",
)
plt.xlabel("$\lambda$ [nm]")
plt.ylabel(r"$\Delta\alpha$ and $\alpha(1V)$ [dBmm$^-$$^1$]")
plt.legend(prop={"size": 9}, loc=1)
plt.rcParams.update({"font.size": 17})
# This variable is used to mark a black dotted line around 12 dBmm^{-1}
a = np.ones(len(data_wavelength_vs_loss_10V_T100C_doped[:, 0]))
a *= 12
plt.plot(
    CONVERSION_E_WAVELENGTH /
    data_wavelength_vs_loss_10V_T100C_doped[:, 0], a, "k--"
)
plt.annotate(r"12dBmm$^-$$^1$", xy=(1330, 13))
plt.annotate("a) PD", xy=(1225, 41))
plt.xlim(1200, 1400)
plt.ylim(-1, 45)
plt.legend(prop={"size": 6.5}, loc=1)
plt.rcParams.update({"font.size": 19})
plt.savefig("Fig_2a.png", dpi=1080, bbox_inches="tight")

# This part extracts the data of the non-intentionally-doped sample.
# Note: An M is used in the name of the variables to represent the minus sign
# in temperatures below 0 °C. E.g., -73 °C is M73.

# Data at 21 °C for the non-intentionally-doped sample is loaded.
file_name_1V_T21C_NID = "UCL_355_R2_1661_A2_S2S3_21_DEG_1V_20th_OCT/loss.txt"
path_file_name_1V_T21C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_1V_T21C_NID
)
file_name_10V_T21C_NID = "UCL_355_R2_1661_A2_21_DEG_S2S3_10V_20th_OCT/loss.txt"
path_file_name_10V_T21C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_10V_T21C_NID
)
data_wavelength_vs_loss_1V_T21C_NID = np.loadtxt(
    path_file_name_1V_T21C_NID, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T21C_NID = np.loadtxt(
    path_file_name_10V_T21C_NID, delimiter=",", dtype=complex
)

# Data at 50 °C for the non-intentionally-doped sample is loaded.
file_name_1V_T50C_NID = "loss_180mA_R2_A_50_1V.txt"
path_file_name_1V_T50C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_1V_T50C_NID
)
file_name_10V_T50C_NID = "loss_180mA_R2_A_50_10V_2.txt"
path_file_name_10V_T50C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_10V_T50C_NID
)
data_wavelength_vs_loss_1V_T50C_NID = np.loadtxt(
    path_file_name_1V_T50C_NID, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T50C_NID = np.loadtxt(
    path_file_name_10V_T50C_NID, delimiter=",", dtype=complex
)
file_name_E_T50C = "Energy_50_R2_A2.txt"
path_file_name_E_T50C = os.path.join(
    non_intentionally_doped_sample, file_name_E_T50C)
data_wavelength_vs_loss_1V_and_10V_T50C_NID = np.loadtxt(
    path_file_name_E_T50C, delimiter=",", dtype=complex
)

# Data at 75 °C for the non-intentionally-doped sample is loaded.
file_name_1V_T75C_NID = "75_retake_6_1V/loss.txt"
path_file_name_1V_T75C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_1V_T75C_NID
)
file_name_10V_T75C_NID = "UCL_355_R2_1661_A2_S2S3_75_DEG_10V_take_2/loss.txt"
path_file_name_10V_T75C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_10V_T75C_NID
)
data_wavelength_vs_loss_1V_T75C_NID = np.loadtxt(
    path_file_name_1V_T75C_NID, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T75C_NID = np.loadtxt(
    path_file_name_10V_T75C_NID, delimiter=",", dtype=complex
)

# Data at 100 °C for the non-intentionally-doped sample is loaded.
file_name_1V_T100C_NID = "UCL_2355_R2_1661_A2_100_DEG_S2S3_1V_take_2/loss.txt"
path_file_name_1V_T100C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_1V_T100C_NID
)
file_name_10V_T100C_NID = "UCL_2355_R2_1661_A2_100_DEG_S2S3_10V/loss.txt"
path_file_name_10V_T100C_NID = os.path.join(
    non_intentionally_doped_sample, file_name_10V_T100C_NID
)
data_wavelength_vs_loss_1V_T100C_NID = np.loadtxt(
    path_file_name_1V_T100C_NID, delimiter=",", dtype=complex
)
data_wavelength_vs_loss_10V_T100C_NID = np.loadtxt(
    path_file_name_10V_T100C_NID, delimiter=",", dtype=complex
)

# This part plots the absorption (α) and change in absorption (Δα) for a
# 9V-swing of the sample with the non-intentionally-doped stack at different
# temperatures.
fig1 = plt.figure(2)
ax = fig1.add_subplot(111)
k = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T21C_NID[:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T21C_NID[:, 1],
    "b",
    label=r"$\alpha(1V)$ @ 21$^\circ$C",
)
l = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T21C_NID[:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T21C_NID[:, 1]
        - data_wavelength_vs_loss_1V_T21C_NID[:, 1]
    ),
    "b--",
    label=r"$\Delta\alpha(9V-swing)$ @ 21$^\circ$C",
)
m = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_and_10V_T50C_NID,
    -0.434 * data_wavelength_vs_loss_1V_T50C_NID,
    "g",
    label=r"$\alpha(1V)$ @ 50$^\circ$C",
)
n = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_and_10V_T50C_NID,
    -0.434
    * (data_wavelength_vs_loss_10V_T50C_NID -
       data_wavelength_vs_loss_1V_T50C_NID),
    "g--",
    label=r"$\Delta\alpha(9V-swing)$ @ 50$^\circ$C",
)
o = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T75C_NID[:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T75C_NID[:, 1],
    "y",
    label=r"$\alpha(1V)$ @ 75$^\circ$C",
)
p = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T75C_NID[:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T75C_NID[:, 1]
        - data_wavelength_vs_loss_1V_T75C_NID[:, 1]
    ),
    "y--",
    label=r"$\Delta\alpha(9V-swing)$ @ 75$^\circ$C",
)
q = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_1V_T100C_NID[50:, 0],
    -0.434 * data_wavelength_vs_loss_1V_T100C_NID[50:, 1],
    "r",
    label=r"$\alpha(1V)$ @ 100$^\circ$C",
)
r = ax.plot(
    CONVERSION_E_WAVELENGTH / data_wavelength_vs_loss_10V_T100C_NID[50:, 0],
    -0.434
    * (
        data_wavelength_vs_loss_10V_T100C_NID[50:, 1]
        - data_wavelength_vs_loss_1V_T100C_NID[49:, 1]
    ),
    "r--",
    label=r"$\Delta\alpha(9V-swing)$ @ 100$^\circ$C",
)
plt.xlim(1200, 1440)
plt.ylim(-1, 50)
ax.set_xlabel("$\lambda$ [nm]")
ax.set_ylabel("$\alpha$ [dBmm$^-$$^1$]")
lns = k + l + m + n + o + p + q + r
labs = [s.get_label() for s in lns]
ax.legend(lns, labs, prop={"size": 9}, loc=1, ncol=2)
ax.set_xlim(1200, 1400)
ax.set_ylim(-1, 45)
ax.set_ylabel(r"$\Delta\alpha$ and $\alpha(1V)$ [dBmm$^-$$^1$]")
# This variable is used to mark a black dotted line around 12 dBmm^{-1}
reference_absorption = np.ones(
    len(data_wavelength_vs_loss_10V_T100C_NID[:, 0]))
reference_absorption *= 12
plt.plot(
    CONVERSION_E_WAVELENGTH /
    data_wavelength_vs_loss_10V_T100C_NID[:, 0], a, "k--"
)
plt.ylabel(r"$\Delta\alpha$ and $\alpha(1V)$ [dBmm$^-$$^1$]")
plt.legend(prop={"size": 7}, loc=1)
plt.annotate(r"12dBmm$^-$$^1$", xy=(1330, 13))
plt.annotate("b) NID", xy=(1225, 41))
plt.savefig("Fig_2b.png", dpi=1080, bbox_inches="tight")

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
