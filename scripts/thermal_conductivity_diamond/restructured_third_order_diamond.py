import os

from PhoGap import PhononsQuippyFiniteDisplacements
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath

path_to_potential = "/CECI/home/ucl/modl/jgeorge/Potentials"

# Please adapt the following lines according to your folder structure.
potential_filenames = [str(os.path.join(path_to_potential, "gp_iter6_sparse9k.xml")),
                       str(os.path.join(path_to_potential, "Fig_5/GAP-18_plus_SCs_f0.01/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_5/GAP-18_plus_SCs_f0.001/gp_iter6C.xml"))]






                   

# This will construct all paths for the potentials, please adapt to the potentials that you would like to use
potential_names = ["GAP18",
                   "Fig_5_GAP-18_plus_SCs_f0.01",
                   "Fig_5_GAP-18_plus_SCs_f0.001"]
structure_mp_149 = Structure.from_file("POSCAR_diamond")

# get conventional cell
kpath = HighSymmKpath(structure_mp_149, symprec=0.01)
structure_mp_149 = kpath.conventional

# displacement for phono3py
DISPLACEMENT_DISTANCE = 0.03

structures = [structure_mp_149]
structure_names = ["mp-149"]
# supercell matrix
supercells_matrix = [[[2, 0, 0], [0, 2, 0], [0, 0, 2]]]

# reference data from DFT computation
reference_data_mean = [
    5748.268,
    4529.575,
    3617.478,
    2924.504,
    2393.040,
    1982.287,
    1662.283,
    1410.767,
    1211.176,
    1051.193,
    921.646,
    815.682,
    728.153,
    655.168,
    593.762,
    541.655,
    497.085,
    458.672,
    425.331,
    396.199,
    370.587,
    347.938,
    327.800,
    309.802,
    293.640,
    279.061,
    265.855,
    253.845,
    242.882,
    232.840,
    223.611,
    215.103,
    207.236,
    199.943,
    193.164,
    186.847,
    180.947,
    175.424,
    170.244,
    165.376,
    160.792,
    156.469,
    152.384,
    148.518,
    144.854,
    141.376,
    138.071,
    134.925,
    131.927,
    129.068,
    126.336,
    123.725,
    121.225,
    118.830,
    116.533,
    114.328,
    112.209,
    110.172,
    108.212,
    106.324,
    104.504,
    102.749,
    101.055,
    99.418,
    97.837,
    96.308,
    94.828,
    93.395,
    92.007,
    90.662,
    89.357,
    88.091,
    86.862,
    85.669,
    84.509,
    83.382,
    82.286,
    81.219,
    80.181,
    79.170,
    78.186,
    77.227,
    76.292,
    75.380,
    74.490,
    73.623,
    72.776,
    71.948,
    71.141,
    70.352,
    69.580,
    68.826,
    68.089,
    67.368,
    66.663,
    65.973,
    65.297,
    64.636,
    63.988,
    63.353,
    62.732,
    62.123,
    61.526,
    60.940,
    60.366,
    59.803,
    59.251,
    58.709,
    58.177,
    57.655,
    57.143,
    56.640,
    56.146,
    55.660,
    55.184,
    54.715,
    54.255,
    53.802,
    53.357,
    52.920,
    52.490,
    52.067,
    51.651,
    51.241,
    50.839,
    50.442,
    50.052,
    49.668,
    49.290,
    48.918,
    48.552,
    48.191,
    47.835,
    47.485,
    47.140,
    46.801,
    46.466,
    46.136,
    45.810,
    45.490,
    45.174,
    44.862,
    44.555,
    44.252,
    43.953,
    43.658,
    43.367,
    43.080,
    42.797,
    42.518,
    42.242,
    41.970,
    41.702,
    41.437,
    41.175,
    40.917,
    40.662,
    40.410,
    40.161,
    39.916,
    39.673,
    39.434,
    39.197,
    38.963,
    38.732,
    38.504,
    38.279,
    38.056,
    37.836,
    37.618,
    37.403,
    37.190,
    36.980,
    36.772,
    36.567,
    36.364,
    36.163,
    35.964,
    35.768,
    35.573,
    35.381,
    35.191,
    35.003,
    34.817,
    34.633,
    34.451,
    34.271,
    34.093,
    33.917,
    33.742,
    33.569
]

reference_data_T = [
    50.0,
    55.0,
    60.0,
    65.0,
    70.0,
    75.0,
    80.0,
    85.0,
    90.0,
    95.0,
    100.0,
    105.0,
    110.0,
    115.0,
    120.0,
    125.0,
    130.0,
    135.0,
    140.0,
    145.0,
    150.0,
    155.0,
    160.0,
    165.0,
    170.0,
    175.0,
    180.0,
    185.0,
    190.0,
    195.0,
    200.0,
    205.0,
    210.0,
    215.0,
    220.0,
    225.0,
    230.0,
    235.0,
    240.0,
    245.0,
    250.0,
    255.0,
    260.0,
    265.0,
    270.0,
    275.0,
    280.0,
    285.0,
    290.0,
    295.0,
    300.0,
    305.0,
    310.0,
    315.0,
    320.0,
    325.0,
    330.0,
    335.0,
    340.0,
    345.0,
    350.0,
    355.0,
    360.0,
    365.0,
    370.0,
    375.0,
    380.0,
    385.0,
    390.0,
    395.0,
    400.0,
    405.0,
    410.0,
    415.0,
    420.0,
    425.0,
    430.0,
    435.0,
    440.0,
    445.0,
    450.0,
    455.0,
    460.0,
    465.0,
    470.0,
    475.0,
    480.0,
    485.0,
    490.0,
    495.0,
    500.0,
    505.0,
    510.0,
    515.0,
    520.0,
    525.0,
    530.0,
    535.0,
    540.0,
    545.0,
    550.0,
    555.0,
    560.0,
    565.0,
    570.0,
    575.0,
    580.0,
    585.0,
    590.0,
    595.0,
    600.0,
    605.0,
    610.0,
    615.0,
    620.0,
    625.0,
    630.0,
    635.0,
    640.0,
    645.0,
    650.0,
    655.0,
    660.0,
    665.0,
    670.0,
    675.0,
    680.0,
    685.0,
    690.0,
    695.0,
    700.0,
    705.0,
    710.0,
    715.0,
    720.0,
    725.0,
    730.0,
    735.0,
    740.0,
    745.0,
    750.0,
    755.0,
    760.0,
    765.0,
    770.0,
    775.0,
    780.0,
    785.0,
    790.0,
    795.0,
    800.0,
    805.0,
    810.0,
    815.0,
    820.0,
    825.0,
    830.0,
    835.0,
    840.0,
    845.0,
    850.0,
    855.0,
    860.0,
    865.0,
    870.0,
    875.0,
    880.0,
    885.0,
    890.0,
    895.0,
    900.0,
    905.0,
    910.0,
    915.0,
    920.0,
    925.0,
    930.0,
    935.0,
    940.0,
    945.0,
    950.0,
    955.0,
    960.0,
    965.0,
    970.0,
    975.0,
    980.0,
    985.0,
    990.0,
    995.0,
    1000.0]

# this should not have an influence on the overall calculation
npoints_band = 51
kpoint_density = 12000

# kpoint density for phono3py calculation
kpoint_density_phono3py = 12000
foldername = "Results"
with open(os.path.join(foldername, "results.txt"), 'w') as f:
    f.write("Pot Structure RMS RMS300 Kappa300 Kappa600 Kappa900\n")
start_from_files = False

for ipot, pot in enumerate(potential_filenames[0:], 0):
    try:
        from quippy import *

        potential = quippy.potential.Potential("IP GAP", param_filename=pot)
    except ModuleNotFoundError:
        start_from_files = True

    for istructure, structure in enumerate(structures):
        print(structure_names[istructure])

        runner = PhononsQuippyFiniteDisplacements(structure, potential, smat=supercells_matrix[istructure],
                                                  path_parameters=os.path.join(foldername, structure_names[
                                                      istructure] + "_" + str(
                                                      potential_names[ipot]) + '_gap_phonopy.yaml'),
                                                  displacementdistance=DISPLACEMENT_DISTANCE,
                                                  npoints_band=npoints_band, kpoint_density=kpoint_density,
                                                  kpoint_density_phono3py=kpoint_density_phono3py,
                                                  set_thermal_conductivity=True, work_with_primitive=False)

        runner.run_all()
        runner.save_kappa_plot(mean_xx_yy_zz=True, xx=False, yy=False, zz=False, yz=False, xz=False, xy=False,
                               filename=os.path.join(foldername, structure_names[istructure] + "_" + str(
                                   potential_names[ipot]) + '_kappa.eps'), mean_exp_data_T=reference_data_T,
                               mean_exp_data=reference_data_mean)
        rms = runner.calculate_rms(exp_data_x=reference_data_T, exp_data_y=reference_data_mean)

        rms_xyz = runner.calculate_rms_xyz(exp_data_T=reference_data_T, exp_data_xx=reference_data_mean,
                                           exp_data_yy=reference_data_mean, exp_data_zz=reference_data_mean)
        rms_300 = runner.calculate_rms(exp_data_x=[300], exp_data_y=[126.336])

        kappa_300 = runner.kappa_at_temperature(temperature=300, whichvalue='xx')
        kappa_600 = runner.kappa_at_temperature(temperature=600, whichvalue='xx')
        kappa_900 = runner.kappa_at_temperature(temperature=900, whichvalue='xx')

        with open(os.path.join(foldername, "results.txt"), 'a') as f:
            f.write(str(potential_names[ipot]) + " " + str(structure_names[istructure]) + " " + str(rms) + " " + str(
                rms_300) + " " + str(kappa_300) + " " + str(kappa_600) + " " + str(kappa_900) + "\n")
        os.rename("kappa-m111111.hdf5", os.path.join(foldername, potential_names[ipot] + '_kappa-m111111.hdf5'))
