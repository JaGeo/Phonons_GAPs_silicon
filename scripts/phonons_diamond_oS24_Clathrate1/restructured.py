import os

from PhoGap import PhonopyFiniteDisplacements, PhononsQuippyFiniteDisplacements, ComparePhononBS
from pymatgen.core.structure import Structure

# Please adapt the following lines according to your folder structure.
path_to_potential = "/CECI/home/ucl/modl/jgeorge/Potentials"

# This will construct all paths for the potentials, please adapt to the potentials that you would like to use
potential_filenames = [str(os.path.join(path_to_potential, "gp_iter6_sparse9k.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/only_disp_nsp_1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/only_disp_nsp_3000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/only_disp_nsp_5000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/only_disp_nsp_7000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/only_disp_nsp_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_MP_nsp_1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_MP_nsp_3000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_MP_nsp_5000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_MP_nsp_7000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_MP_nsp_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_wd_nsp_1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_wd_nsp_3000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_wd_nsp_5000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_wd_nsp_7000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/02_MP_series/scratch_wd_nsp_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/03_fas_series/only_disp_9k_0.0001/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/03_fas_series/only_disp_9k_0.0010/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/03_fas_series/only_disp_9k_0.0100/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/03_fas_series/only_disp_9k_0.1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential,
                                        "07_fit/03_fas_series/scratch_wd_nsp9k_0.0001/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential,
                                        "07_fit/03_fas_series/scratch_wd_nsp9k_0.0010/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential,
                                        "07_fit/03_fas_series/scratch_wd_nsp9k_0.0100/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential,
                                        "07_fit/03_fas_series/scratch_wd_nsp9k_0.1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/04_GAP-18C/nsp_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/04_GAP-18C/nsp_9000_onlydisp/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/04_GAP-18C/nsp_9000_3000/fas_0.0010/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/04_GAP-18C/nsp_9000_3000/fas_0.0100/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential,
                                        "07_fit/06_GAP-18C_wrnd/nsp_9000_3000_0.0100/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "07_fit/06_GAP-18C_wrnd/nsp_9000_3000_0.0010/gp_iter6C.xml"))
                       ]
# This will construct the names of the potentials
potential_names = ["GAP18",
                   "02_MP_series_only_disp_nsp_1000",
                   "02_MP_series_only_disp_nsp_3000",
                   "02_MP_series_only_disp_nsp_5000",
                   "02_MP_series_only_disp_nsp_7000",
                   "02_MP_series_only_disp_nsp_9000",
                   "02_MP_series_scratch_MP_nsp_1000",
                   "02_MP_series_scratch_MP_nsp_3000",
                   "02_MP_series_scratch_MP_nsp_5000",
                   "02_MP_series_scratch_MP_nsp_7000",
                   "02_MP_series_scratch_MP_nsp_9000",
                   "02_MP_series_scratch_wd_nsp_1000",
                   "02_MP_series_scratch_wd_nsp_3000",
                   "02_MP_series_scratch_wd_nsp_5000",
                   "02_MP_series_scratch_wd_nsp_7000",
                   "02_MP_series_scratch_wd_nsp_9000",
                   "03_fas_series_only_disp_9k_0.0001",
                   "03_fas_series_only_disp_9k_0.0010",
                   "03_fas_series_only_disp_9k_0.0100",
                   "03_fas_series_only_disp_9k_0.1000",
                   "03_fas_series_scratch_wd_nsp9k_0.0001",
                   "03_fas_series_scratch_wd_nsp9k_0.0010",
                   "03_fas_series_scratch_wd_nsp9k_0.0100",
                   "03_fas_series_scratch_wd_nsp9k_0.1000",
                   "04_GAP-18C_nsp_9000",
                   "04_GAP-18C_nsp_9000_onlydisp",
                   "04_GAP-18C_nsp_9000_3000_fas_0.0010",
                   "04_GAP-18C_nsp_9000_3000_fas_0.0100",
                   "06_GAP-18C_wrnd_nsp_9000_3000_0.0100",
                   "06_GAP-18C_wrnd_nsp_9000_3000_0.0010"
                   ]


def get_structures_names():
    """
    Function to create all structures
    Returns:

    """
    structure_mp_149 = Structure.from_file("POSCAR_diamond")

    structure_oS24 = Structure.from_spacegroup(sg=63,
                                               lattice=[[3.82236, 0.0, 0.0], [0.0, 10.7007, 0.0], [0.0, 0.0, 12.6258]],
                                               species=["Si", "Si", "Si"],
                                               coords=[[0., 0.243, 0.5551], [0.0, 0.5705, 0.3412],
                                                       [0.0, 0.0284, 0.5903]])
    structure_clathrate_1 = Structure.from_spacegroup(sg=223,
                                                      lattice=[[10.355, 0, 0], [0, 10.355, 0], [0, 0, 10.355]],
                                                      species=["Si", "Si", "Si"],
                                                      coords=[[0, 0.25, 0.5], [0.1837, 0.1837, 0.1837],
                                                              [0, 0.1172, 0.3077]])

    structures = [structure_mp_149, structure_oS24, structure_clathrate_1]
    structure_names = ["mp-149", "oS24", "Clathrate-1"]
    return structures, structure_names


structures, structure_names = get_structures_names()

# will determine which castep calculations are used for comparison
castep_calc_ids = [22, 9, 33]

# these are the supercell matrices that are used for the comparison
supercells_matrix = [[[4, 0, 0], [0, 4, 0], [0, 0, 4]], [[3, 0, 0], [0, 3, 0], [0, 0, 2]],
                     [[2, 0, 0], [0, 2, 0], [0, 0, 2]]]

# number of points for each path in phonon band structure
npoints_band = 51
kpoint_density = 12000

# finite displacement for phonon calculation, here only 0.01 is tested
DISTANCES = [0.01]

for DISPLACEMENT_DISTANCE in DISTANCES:
    with open("Results/results_" + str(DISPLACEMENT_DISTANCE) + ".txt", 'w') as f:
        f.write("Pot Structure RMS\n")

start_from_files = False
for ipot, pot in enumerate(potential_filenames[0:], 0):
    try:
        from quippy import *

        potential = quippy.potential.Potential("IP GAP", param_filename=pot)
    except ModuleNotFoundError:
        start_from_files = True
    foldername = "Results"
    ylim = [[0, 16], [0, 16], [0, 16], [0, 16]]

    for istructure, structure in enumerate(structures):
        for DISPLACEMENT_DISTANCE in DISTANCES:
            print(structure_names[istructure])

            runner = PhononsQuippyFiniteDisplacements(structure, potential, smat=supercells_matrix[istructure],
                                                      path_parameters=os.path.join(foldername, structure_names[
                                                          istructure] + "_" + str(
                                                          potential_names[ipot]) + '_gap_phonopy.yaml'),
                                                      displacementdistance=DISPLACEMENT_DISTANCE,
                                                      npoints_band=npoints_band, kpoint_density=kpoint_density)

            runner_castep = PhonopyFiniteDisplacements(
                filename=os.path.join("../DFT_Benchmark_CASTEP", str(castep_calc_ids[istructure]) + '_phonopy.yaml'),
                kpoint_density=kpoint_density)

            runner.run_all()

            runner.save_plot_band_structure(os.path.join(foldername, structure_names[
                istructure] + "_" + str(potential_names[ipot]) + "_" + str(
                DISPLACEMENT_DISTANCE) + '_' + "_gap_band_structure.eps"), ylim=ylim[istructure])
            runner.save_plot_dos(os.path.join(foldername, structure_names[
                istructure] + "_" + str(potential_names[ipot]) + "_" + str(
                DISPLACEMENT_DISTANCE) + '_' + "_gap_dos_structure.eps"))

            runner_castep.run_all()

            runner_castep.save_plot_band_structure(os.path.join(foldername, structure_names[
                istructure] + "_castep_band_structure.eps"), ylim=ylim[istructure])
            runner_castep.save_plot_dos(os.path.join(foldername, structure_names[
                istructure] + "_castep_dos_structure.eps"))

            try:
                comparison = ComparePhononBS(runner_castep.phonon_band_structure_pymatgen,
                                             runner.phonon_band_structure_pymatgen)
                rms = comparison.rms_overall()
                # To check if frequencies in phonopy are correctly sorted
                rms2 = comparison.rms_overall_second_definition()
                comparison.rms_kdep_plot(whichkpath=1, filename=os.path.join(foldername,
                                                                             structure_names[istructure] + '_' + str(
                                                                                 potential_names[ipot]) + "_" + str(
                                                                                 DISPLACEMENT_DISTANCE) + '_rms_phonons.eps'))
                comparison.compare_plot(filename=os.path.join(foldername, structure_names[istructure] + '_' + str(
                    potential_names[ipot]) + "_" + str(
                    DISPLACEMENT_DISTANCE) + '_comparison_phonons.eps'))
                MARE = comparison.mean_absolute_error()

            except:
                rms = "None"
                rms2 = "None"
                MARE = "None"
            with open("Results/results_" + str(DISPLACEMENT_DISTANCE) + ".txt", 'a') as f:
                f.write(str(potential_names[ipot]) + ' ' + str(structure_names[istructure]) + ' ' + str(rms) + '\n')
