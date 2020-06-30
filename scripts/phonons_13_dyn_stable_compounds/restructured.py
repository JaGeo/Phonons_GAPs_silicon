import os

from PhoGap import PhonopyFiniteDisplacements, PhononsQuippyFiniteDisplacements, ComparePhononBS

# Please adapt the following lines according to your folder structure.
path_to_potential = "/CECI/home/ucl/modl/jgeorge/Potentials"

# This will construct all paths for the potentials, please adapt to the potentials that you would like to use
potential_filenames = [str(os.path.join(path_to_potential, "gp_iter6_sparse9k.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_random/M_1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_random/M_3000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_random/M_5000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_random/M_7000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_random/M_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_individual/M_1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_individual/M_3000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_individual/M_5000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_individual/M_7000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/only_individual/M_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/combined/M_1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/combined/M_3000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/combined/M_5000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/combined/M_7000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_2/combined/M_9000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/only_individual/f_0.1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/only_individual/f_0.0100/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/only_individual/f_0.0010/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/only_individual/f_0.0001/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/combined/f_0.1000/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/combined/f_0.0100/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/combined/f_0.0010/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_4/combined/f_0.0001/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_5/GAP-18_plus_SCs_f0.01/gp_iter6C.xml")),
                       str(os.path.join(path_to_potential, "Fig_5/GAP-18_plus_SCs_f0.001/gp_iter6C.xml"))
                       ]

# This will construct the names of the potentials
potential_names = ["GAP18",
                    "Fig_2_only_random_M_1000",
                    "Fig_2_only_random_M_3000",
                    "Fig_2_only_random_M_5000",
                    "Fig_2_only_random_M_7000",
                    "Fig_2_only_random_M_9000",
                    "Fig_2_only_individual_M_1000",
                    "Fig_2_only_individual_M_3000",
                    "Fig_2_only_individual_M_5000",
                    "Fig_2_only_individual_M_7000",
                    "Fig_2_only_individual_M_9000",
                    "Fig_2_combined_M_1000",
                    "Fig_2_combined_M_3000",
                    "Fig_2_combined_M_5000",
                    "Fig_2_combined_M_7000",
                    "Fig_2_combined_M_9000",
                    "Fig_4_only_individual_f_0.1000",
                    "Fig_4_only_individual_f_0.0100",
                    "Fig_4_only_individual_f_0.0010",
                    "Fig_4_only_individual_f_0.0001",
                    "Fig_4_combined_f_0.1000",
                    "Fig_4_combined_f_0.0100",
                    "Fig_4_combined_f_0.0010",
                    "Fig_4_combined_f_0.0001",
                    "Fig_5_GAP-18_plus_SCs_f0.01",
                    "Fig_5_GAP-18_plus_SCs_f0.001"]
                   
                        

                 

# number of points per path in phonon band structure
npoints_band = 51
kpoint_density = 12000

# finite displacement for phonon calculation
DISTANCES = [0.01]

for DISPLACEMENT_DISTANCE in DISTANCES:
    with open("Results/results_" + str(DISPLACEMENT_DISTANCE) + ".txt", 'w') as f:
        f.write("Pot Structure mpid RMS imagmodes(pot) imagmodes(dft) \n")

start_from_files = False
for ipot, pot in enumerate(potential_filenames[0:], 0):
    try:
        from quippy import *

        potential = quippy.potential.Potential("IP GAP", param_filename=pot)
    except ModuleNotFoundError:
        # Just for testing if quippy is not installed
        from FakeCalculator import FakeCalc

        potential = FakeCalc
        print("Will use fake calculator")
    foldername = "Results"
    # mpids of the structures that are analysed
    mpids = ["mp-1072544", "mp-1079297", "mp-1095269", "mp-1200830", "mp-1203790", "mp-149", "mp-16220", "mp-165",
             "mp-168", "mp-571520", "mp-971661", "mp-971662", "mp-999200"]

    for iistructure, istructure in enumerate([5, 6, 9, 15, 18, 22, 23, 24, 25, 28, 32, 33, 34]):
        # 5: mp-1072544
        # 6: mp-1079297
        # 9: mp-1095269
        # 15: mp-1200830
        # 18: mp-1203790
        # 22: mp-149
        # 23: mp-16220
        # 24: mp-165
        # 25: mp-168
        # 28: mp-571520
        # 32: mp-971661
        # 33: mp-971662
        # 34: mp-999200
        for DISPLACEMENT_DISTANCE in DISTANCES:

            runner_castep = PhonopyFiniteDisplacements(
                filename=os.path.join("../DFT_Benchmark_CASTEP", str(istructure) + '_phonopy.yaml'),
                kpoint_density=kpoint_density)
            runner_castep.run_all()

            smat = runner_castep.phonon.get_supercell_matrix().tolist()
            runner = PhononsQuippyFiniteDisplacements(runner_castep.optimized_structure, potential, smat=smat,
                                                      path_parameters=os.path.join(foldername,
                                                                                   str(istructure) + "_" + str(
                                                                                       potential_names[
                                                                                           ipot]) + '_gap_phonopy.yaml'),
                                                      displacementdistance=DISPLACEMENT_DISTANCE,
                                                      npoints_band=npoints_band, kpoint_density=kpoint_density)

            runner.run_all()

            runner.save_plot_band_structure(
                os.path.join(foldername, str(istructure) + "_" + str(potential_names[ipot]) + "_" + str(
                    DISPLACEMENT_DISTANCE) + "_gap_band_structure.eps"), ylim=[-1, 16])
            runner.save_plot_dos(
                os.path.join(foldername, str(istructure) + "_" + str(potential_names[ipot]) + "_" + str(
                    DISPLACEMENT_DISTANCE) + "_gap_dos_structure.eps"))

            runner_castep.save_plot_band_structure(
                os.path.join(foldername, str(istructure) + "_castep_band_structure.eps"), ylim=[-1, 16])
            runner_castep.save_plot_dos(os.path.join(foldername, str(istructure) + "_castep_dos_structure.eps"))

            try:
                comparison = ComparePhononBS(runner.phonon_band_structure_pymatgen,
                                             runner_castep.phonon_band_structure_pymatgen)
                rms = comparison.rms_overall()
                comparison.rms_kdep_plot(whichkpath=2, filename=os.path.join(foldername,
                                                                             str(istructure) + '_' + str(
                                                                                 potential_names[ipot]) + "_" + str(
                                                                                 DISPLACEMENT_DISTANCE) + '_''_rms_phonons.eps'))
            except:
                rms = "NONE"

            with open("Results/results_" + str(DISPLACEMENT_DISTANCE) + ".txt", 'a') as f:
                f.write(str(potential_names[ipot]) + ' ' + str(istructure) + ' ' + str(mpids[iistructure]) + ' ' + str(
                    rms) + ' ' + str(
                    runner.has_imag_modes(0.1)) + ' ' + str(runner_castep.has_imag_modes(0.1)) + '\n')
