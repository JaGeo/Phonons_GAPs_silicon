import unittest
from unittest.mock import MagicMock

from pymatgen.ext.matproj import MPRester
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos

from FakeCalculator import FakeCalc
from PhoGap import *

try:
    from phono3py import Phono3py
except:
    pass

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
testfile_dir = os.path.join(module_dir, "test_files")


class Test_HelperFunctions(unittest.TestCase):
    def setUp(self) -> None:
        with MPRester() as a:
            self.structure = a.get_structure_by_material_id('mp-149')
            self.phonon_band = a.get_phonon_bandstructure_by_material_id("mp-149")

    def test_get_prim(self):
        prim_structure = get_prim(structure=self.structure)
        self.assertEqual(prim_structure.num_sites, 2)

    def test_get_kpath(self):
        prim_structure = get_prim(structure=self.structure)
        k1, k2 = get_kpath(prim_structure)
        # check if each key is in k1["kpoints"]
        # print(k2)
        for key, qpoint in self.phonon_band.labels_dict.items():
            found = False
            if key in k1["kpoints"]:
                found = True
            self.assertTrue(found)

    def test_get_ase_atoms(self):
        atoms = get_ase_atoms(self.structure)
        self.assertAlmostEqual(self.structure.lattice.volume, atoms.get_volume())
        lengths_and_angles = atoms.get_cell_lengths_and_angles()
        for iabc, abc in enumerate(self.structure.lattice.abc):
            self.assertAlmostEqual(lengths_and_angles[iabc], abc)
        for iabc, alphabetagamma in enumerate(self.structure.lattice.angles):
            self.assertAlmostEqual(lengths_and_angles[3 + iabc], alphabetagamma)

    def test_get_structure(self):
        atoms = get_ase_atoms(self.structure)
        structure = get_structure(atoms)
        for iabc, abc in enumerate(self.structure.lattice.abc):
            self.assertAlmostEqual(structure.lattice.abc[iabc], abc)


class TestCompareZPE(unittest.TestCase):

    def setUp(self) -> None:
        phonopy1 = MagicMock()
        phonopy1.get_zero_point_energy = MagicMock(return_value=1)
        phonopy1.get_free_energy = MagicMock(return_value=2)

        # mock get_free_energy and get_zero_point_energy
        values = {0: 1, 100: 0}

        def side_effect(temperature):
            return values[temperature]

        phonopy2 = MagicMock()
        phonopy2.get_zero_point_energy = MagicMock(return_value=2)
        phonopy2.get_free_energy = side_effect
        self.compare2 = CompareZPE(phonopy1, phonopy2)

    def test_rms_zpe(self):
        self.assertAlmostEqual(self.compare2.get_rms_zpe(), 1.0)

    def test_rms_free_energy(self):
        self.assertAlmostEqual(self.compare2.get_rms_free_energy(startt=0, stopt=100, stept=100),
                               np.sqrt((2 ** 2 + 1 ** 2) / 2))


class TestComparePhononBs(unittest.TestCase):

    def setUp(self) -> None:
        self.phonopyfinite = PhonopyFiniteDisplacements(filename=os.path.join(testfile_dir, "phonopy.yaml"))
        self.phonopyfinite.run_all()
        self.phonopyfinite2 = PhonopyFiniteDisplacements(filename=os.path.join(testfile_dir, "phonopy.yaml"))
        self.phonopyfinite2.run_all()
        self.phonopyfinite3 = PhonopyFiniteDisplacements(
            filename=os.path.join(testfile_dir, "phonopy_manipulated.yaml"))
        self.phonopyfinite3.run_all()

        self.compare = ComparePhononBS(self.phonopyfinite.phonon_band_structure_pymatgen,
                                       self.phonopyfinite2.phonon_band_structure_pymatgen)

        self.compare2 = ComparePhononBS(self.phonopyfinite.phonon_band_structure_pymatgen,
                                        self.phonopyfinite3.phonon_band_structure_pymatgen)
        self.band1_dict = self.phonopyfinite2.phonon_band_structure_pymatgen.as_dict()
        new_bands = []
        for iband, band in enumerate(self.band1_dict["bands"]):
            new_bands.append([])
            for frequency in band:
                new_bands[iband].append(frequency + 0.1)
        self.band1_dict["bands"] = new_bands
        self.bandnew = PhononBandStructureSymmLine.from_dict(self.band1_dict)
        self.compare3 = ComparePhononBS(self.phonopyfinite.phonon_band_structure_pymatgen,
                                        self.bandnew)

    def test_rms_overall(self):
        self.assertAlmostEqual(self.compare.rms_overall(), 0)
        self.assertNotAlmostEqual(self.compare2.rms_overall(), 0)
        self.assertAlmostEqual(self.compare3.rms_overall(), 0.1)

    def test_rms_overall_2(self):
        self.assertAlmostEqual(self.compare.rms_overall_second_definition(), 0)
        self.assertNotAlmostEqual(self.compare2.rms_overall_second_definition(), 0)
        self.assertAlmostEqual(self.compare3.rms_overall_second_definition(), 0.1)

    def test_kdep_rms(self):
        for el in self.compare.rms_kdep():
            self.assertAlmostEqual(el, 0.0)
        for el in self.compare3.rms_kdep():
            self.assertAlmostEqual(el, 0.1)

    def test_kdep_rms_plot(self):
        self.compare.rms_kdep_plot(1)


class TestPhonopyFiniteDisplacements(unittest.TestCase):

    def setUp(self) -> None:
        self.phonopyfinite = PhonopyFiniteDisplacements(filename=os.path.join(testfile_dir, "phonopy.yaml"))
        self.phonopyfinite.run_all()
        self.phonon = self.phonopyfinite.read_yaml(os.path.join(testfile_dir, "phonopy.yaml"))
        self.structure = get_pmg_structure(self.phonon.get_unitcell())

    def test_read_yaml(self):
        # TODO: rethink tests
        phonon = self.phonopyfinite.read_yaml(os.path.join(testfile_dir, "phonopy.yaml"))
        self.assertTrue(Phonopy, type(phonon))

    def test_get_primitive_cell(self):
        # creates supercell and then tests if it arrives at same primitive cell!
        structure = get_pmg_structure(self.phonon.get_supercell())
        structure_to_compare = get_pmg_structure(self.phonon.get_unitcell())
        structure_prim = self.phonopyfinite._get_primitive_cell(structure)
        structure_prim_to_compare = self.phonopyfinite._get_primitive_cell(structure_to_compare)
        for ilat1, lat1 in enumerate(structure_prim_to_compare.lattice.abc):
            self.assertAlmostEqual(lat1, structure_prim.lattice.abc[ilat1])

    def test_get_ase_from_pmg(self):
        ase = self.phonopyfinite._get_ase_from_pmg(self.structure)
        for ilength, length in enumerate(ase.cell.lengths()):
            self.assertAlmostEqual(length, self.structure.lattice.abc[ilength])

    def test_get_pmg_from_ase(self):
        ase = self.phonopyfinite._get_ase_from_pmg(self.structure)
        structure = self.phonopyfinite._get_pmg_from_ase(ase)
        for ilength, length in enumerate(ase.cell.lengths()):
            self.assertAlmostEqual(length, structure.lattice.abc[ilength])

    def test_get_bandstructure_calc(self):
        bs = self.phonopyfinite._get_bandstructure_calc(structure=self.phonopyfinite.optimized_structure,
                                                        phonon=self.phonopyfinite.phonon)
        self.assertEqual(type(bs), PhononBandStructureSymmLine)

    def test_get_dos_calc(self):
        dos = self.phonopyfinite._get_dos_calc(structure=self.phonopyfinite.optimized_structure,
                                               phonon=self.phonopyfinite.phonon, kpoint_density=100)
        self.assertEqual(type(dos), PhononDos)

    def test_save_plot_band_structure(self):
        with tempfile.TemporaryDirectory() as d:
            filename = os.path.join(d, 'band.eps')
            self.phonopyfinite.save_plot_band_structure(filename=filename)
            self.assertTrue(os.path.exists(filename))

    def test_save_plot_dos(self):
        with tempfile.TemporaryDirectory() as d:
            filename = os.path.join(d, 'band.eps')
            self.phonopyfinite.save_plot_dos(filename)
            self.assertTrue(os.path.exists(filename))
        # os.rmdir(filename)

    def test_get_imaginary_modes(self):
        self.assertFalse(self.phonopyfinite.has_imag_modes(0.1))

    def test_get_zero_point(self):
        self.phonopyfinite.phonon.run_thermal_properties(0, 100, 10)
        # phonopy result in kJ/mol-c, result from pymatgen in J/mol (already per formula unit)
        self.assertAlmostEqual(self.phonopyfinite.phonon.get_thermal_properties_dict()['free_energy'][0] / 2.0,
                               self.phonopyfinite.get_zero_point_energy() / 1000.0, 1)

    def test_get_free_energy(self):
        self.phonopyfinite.phonon.run_thermal_properties(100, 100, 10)
        # phonopy result in kJ/mol-c, result from pymatgen in J/mol (already per formula unit)
        self.assertAlmostEqual(self.phonopyfinite.phonon.get_thermal_properties_dict()['free_energy'][0] / 2.0,
                               self.phonopyfinite.get_free_energy(temperature=100.0) / 1000.0, 1)


# TODO: test this class! make sure it is using the correct supercells
class TestPhononsQuippyFiniteDisplacements(unittest.TestCase):
    def setUp(self) -> None:
        self.structure = Structure.from_file(os.path.join(testfile_dir, "POSCAR"))
        self.structure_conv = Structure.from_file(os.path.join(testfile_dir, "POSCAR.conv"))
        self.phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                              kpoint_density=100, kpoint_density_phono3py=100)

    def test_run_all(self):
        self.phononsquippy.run_all()
        for ilength, length in enumerate(self.phononsquippy.optimized_structure.lattice.abc):
            self.assertNotAlmostEqual(length, self.structure.lattice.abc[ilength])
        self.phononsquippy.energy_optimized_structure
        self.phononsquippy.optimized_structure
        self.phononsquippy.phonon_band_structure_pymatgen
        self.phononsquippy.phonon_dos_pymatgen
        self.phononsquippy.phonon

    def test_parameters_kpoint_density(self):
        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=10, npoints_band=1)
        phononsquippy.run_all()

        phononsquippy2 = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                          kpoint_density=100, npoints_band=1)
        phononsquippy2.run_all()

        notequal = False
        for density1, density2 in zip(phononsquippy.phonon_dos_pymatgen.densities,
                                      phononsquippy2.phonon_dos_pymatgen.densities):
            if abs(density1 - density2) > 1e-5:
                notequal = True
        self.assertTrue(notequal)

    def test_parameters_npoints_band(self):
        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=10, npoints_band=10)
        phononsquippy.run_all()
        phononsquippy2 = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                          kpoint_density=100, npoints_band=1)
        phononsquippy2.run_all()
        notequal = False
        for qpoint1, qpoint2 in zip(phononsquippy.phonon_band_structure_pymatgen.qpoints,
                                    phononsquippy2.phonon_band_structure_pymatgen.qpoints):
            if (qpoint1.frac_coords[0] != qpoint2.frac_coords[0]) or (
                    qpoint1.frac_coords[1] != qpoint2.frac_coords[1]) or qpoint1.frac_coords[2] != qpoint2.frac_coords[
                2]:
                notequal = True
        self.assertTrue(notequal)

    def test_parameters_work_with_primitive(self):

        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure_conv, potential=FakeCalc,
                                                         kpoint_density=10, npoints_band=1, work_with_primitive=False,
                                                         set_phonons=False)
        phononsquippy.run_all()
        self.assertEqual(phononsquippy.optimized_structure.num_sites, self.structure_conv.num_sites)

    def test_parameters_kpoint_density_phono3py(self):
        try:
            if os.path.exists("kappa-m111.hdf5"):
                ValueError("A kappa-m111.hdf5 exists. Please remove it before you restart the test")
            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=10, kpoint_density_phono3py=10,
                                                             set_thermal_conductivity=True,
                                                             smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                             max_distance_third_order=1.5)
            phononsquippy.run_all()
            self.assertTrue(os.path.exists("kappa-m111.hdf5"))
            os.remove("kappa-m111.hdf5")



        except ValueError:
            ValueError("A kappa-m111.hdf5 exists. Please remove it before you restart the test")
        except:
            print("phono3py not implemented")

    def test_parameters_displacementdistancephono3py(self):

        try:

            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=10, kpoint_density_phono3py=10,
                                                             set_thermal_conductivity=True,
                                                             smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                             max_distance_third_order=4,
                                                             displacementdistancephono3py=0.01)

            phononsquippy.run_all()
            phononsquippy2 = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                              kpoint_density=10, kpoint_density_phono3py=10,
                                                              set_thermal_conductivity=True,
                                                              smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                              max_distance_third_order=1.5,
                                                              displacementdistancephono3py=0.1)
            phononsquippy2.run_all()

            phono3py1 = phononsquippy.phono3py
            phono3py2 = phononsquippy2.phono3py

            # test if at least one coordinate in each structure is influenced by the displacement setting

            notthesame = False

            for scell1, scell2 in zip(phono3py1.get_supercells_with_displacements(),
                                      phono3py2.get_supercells_with_displacements()):
                structure1 = get_pmg_structure(scell1)
                structure2 = get_pmg_structure(scell2)
                for atoms1, atoms2 in zip(structure1, structure2):
                    for xyz1, xyz2 in zip(atoms1.coords, atoms2.coords):
                        if abs(xyz1 - xyz2) > 1e-3:
                            notthesame = True
            self.assertTrue(notthesame)

        except:
            print("Phono3py not available")

    def test_parameters_smat(self):
        smat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=10, kpoint_density_phono3py=10,
                                                         set_thermal_conductivity=False, smat=smat)
        phononsquippy.run_all()

        for firstarray, secondarray in zip(phononsquippy.phonon.get_supercell_matrix(), smat):
            for element, element2 in zip(firstarray, secondarray):
                self.assertEqual(element, element2)

        smat = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=10, kpoint_density_phono3py=10,
                                                         set_thermal_conductivity=False, smat=smat)
        phononsquippy.run_all()

        for firstarray, secondarray in zip(phononsquippy.phonon.get_supercell_matrix(), smat):
            for element, element2 in zip(firstarray, secondarray):
                self.assertEqual(element, element2)

    def test_parameters_temperature_range_kappa(self):
        try:

            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=10, kpoint_density_phono3py=10,
                                                             set_thermal_conductivity=True,
                                                             smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                             max_distance_third_order=1.5, temperature_range_kappa=[0])
            phononsquippy.run_all()
            print(phononsquippy.kappa)
            self.assertEqual(len(phononsquippy.kappa), 1)
        except:
            print("phono3py not available on your system")

    def test_max_distance_third_order(self):
        try:

            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=10, kpoint_density_phono3py=10,
                                                             set_thermal_conductivity=True,
                                                             smat=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                                                             max_distance_third_order=1.5, temperature_range_kappa=[0])

            phononsquippy.run_all()
            phono3py1 = phononsquippy.phono3py
            self.assertTrue(None in phono3py1.get_supercells_with_displacements())


        except:
            print("phono3py not imported")

    def test_kappa_at_temperature(self):
        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=100, kpoint_density_phono3py=100)
        phononsquippy.temperature_range_kappa = [100, 200]
        phononsquippy.kappa_xx = [1, 0]
        phononsquippy.kappa_yy = [0, 2]
        phononsquippy.kappa_zz = [3, 0]
        phononsquippy.kappa_xz = [4, 0]
        phononsquippy.kappa_xy = [5, 0]
        phononsquippy.kappa_yz = [6, 0]

        self.assertEqual(phononsquippy.kappa_at_temperature(100, 'xx'), 1)
        self.assertEqual(phononsquippy.kappa_at_temperature(200, 'yy'), 2)
        self.assertEqual(phononsquippy.kappa_at_temperature(100, 'zz'), 3)
        self.assertEqual(phononsquippy.kappa_at_temperature(100, 'xy'), 5)
        self.assertEqual(phononsquippy.kappa_at_temperature(100, 'yz'), 6)
        self.assertEqual(phononsquippy.kappa_at_temperature(100, 'xz'), 4)

    def test_calculate_rms(self):

        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=100, kpoint_density_phono3py=100)
        phononsquippy.temperature_range_kappa = [100, 200]
        phononsquippy.kappa_mean = [1, 0]
        self.assertAlmostEqual(phononsquippy.calculate_rms(exp_data_x=[100, 200], exp_data_y=[1, 0]), 0.0)
        self.assertAlmostEqual(phononsquippy.calculate_rms(exp_data_x=[100, 200], exp_data_y=[2, 1]), 1.0)
        self.assertAlmostEqual(phononsquippy.calculate_rms(exp_data_x=[100], exp_data_y=[1]), 0.0)
        self.assertAlmostEqual(phononsquippy.calculate_rms(exp_data_x=[200], exp_data_y=[1]), 1.0)

    def test_calculate_rms_xyz(self):
        phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                         kpoint_density=100, kpoint_density_phono3py=100)
        phononsquippy.temperature_range_kappa = [100, 200]
        phononsquippy.kappa_xx = [1, 0]
        phononsquippy.kappa_yy = [1, 0]
        phononsquippy.kappa_zz = [1, 0]

        self.assertAlmostEqual(
            phononsquippy.calculate_rms_xyz(exp_data_T=[100, 200], exp_data_xx=[1, 0], exp_data_yy=[1, 0],
                                            exp_data_zz=[1, 0]), 0.0)
        self.assertAlmostEqual(
            phononsquippy.calculate_rms_xyz(exp_data_T=[100, 200], exp_data_xx=[2, 1], exp_data_yy=[2, 1],
                                            exp_data_zz=[2, 1]), 1.0)
        self.assertAlmostEqual(
            phononsquippy.calculate_rms_xyz(exp_data_T=[100], exp_data_xx=[1], exp_data_yy=[1], exp_data_zz=[1]), 0.0)
        self.assertAlmostEqual(
            phononsquippy.calculate_rms_xyz(exp_data_T=[200], exp_data_xx=[1], exp_data_yy=[1], exp_data_zz=[1]), 1.0)

    def test_save_kappa_plot(self):
        with tempfile.TemporaryDirectory() as d:
            filename = os.path.join(d, 'kappa.eps')
            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=100, kpoint_density_phono3py=100)
            phononsquippy.temperature_range_kappa = [100, 200]
            phononsquippy.kappa_mean = [0, 1]
            phononsquippy.save_kappa_plot(mean_xx_yy_zz=True, xx=False, yy=False, zz=False, xy=False, xz=False,
                                          yz=False, filename=filename)
            self.assertTrue(os.path.exists(filename))

    def test_get_thermal_conductivity_matrix(self):
        try:
            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=100, kpoint_density_phono3py=100,
                                                             set_thermal_conductivity=True,
                                                             smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                             max_distance_third_order=1.5)
            phononsquippy.run_all()

            self.assertEqual(type(phononsquippy._get_thermal_conductivity_matrix([100])), np.ndarray)



        except:
            print("phono3py not implemented")

    def test_get_phono3pyobject_phono3py(self):
        try:
            phononsquippy = PhononsQuippyFiniteDisplacements(structure=self.structure, potential=FakeCalc,
                                                             kpoint_density=100, kpoint_density_phono3py=100,
                                                             set_thermal_conductivity=True,
                                                             smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                             max_distance_third_order=1.5)
            phononsquippy.run_all()
            self.assertEqual(phononsquippy.get_phono3pyobject_phono3py(), Phono3py)
        except:
            print("phono3py not implemented")

    def test_get_optimized_cell(self):
        struct_opt = self.phononsquippy._get_optimized_cell(self.structure, FakeCalc)
        for ilength, length in enumerate(struct_opt.lattice.abc):
            self.assertNotAlmostEqual(length, self.structure.lattice.abc[ilength])

    def test_get_potential_energy(self):
        energy = self.phononsquippy._get_potential_energy(self.structure, FakeCalc)
        # just make sure it does something
        self.assertAlmostEqual(energy, -0.06831434931360365)

    def test_get_phononobject_phonopy(self):
        tempfilename = os.path.join(tempfile.gettempdir(), "result.yaml")
        smat = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        phonopyobject = self.phononsquippy._get_phononobject_phonopy(self.structure, FakeCalc,
                                                                     smat=smat,
                                                                     displacement_distance=0.01, save_parameters=True,
                                                                     path=tempfilename)
        self.assertEqual(type(phonopyobject), Phonopy)

        new_phonopy = load(tempfilename)
        for iel1, el1 in enumerate(new_phonopy.get_supercell_matrix()):
            for iel2, el2 in enumerate(el1):
                self.assertEqual(smat[iel1][iel2], el2)

        phonopyobject2 = self.phononsquippy._get_phononobject_phonopy(self.structure, FakeCalc,
                                                                      smat=smat,
                                                                      displacement_distance=0.02, save_parameters=True,
                                                                      path=tempfilename)
        self.assertAlmostEqual(new_phonopy.get_displacements()[0][2] * 2.0, phonopyobject2.get_displacements()[0][2], 3)


if __name__ == '__main__':
    unittest.main()
