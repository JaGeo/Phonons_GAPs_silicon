import copy
import math
import os
import tempfile
import warnings

import matplotlib as mpl
import numpy as np
from ase import Atoms
from ase.constraints import UnitCellFilter, StrainFilter
from ase.optimize import BFGS
from phonopy import Phonopy, load
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.units import VaspToTHz
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_ph_dos
from pymatgen.io.phonopy import get_phonopy_structure, get_ph_bs_symm_line, get_pmg_structure
from pymatgen.io.vasp import Kpoints
from pymatgen.phonon.plotter import PhononBSPlotter, PhononDosPlotter
from pymatgen.symmetry.bandstructure import HighSymmKpath

try:
    from quippy import *
except:
    pass
try:
    from phono3py import Phono3py
    from phono3py.phonon3.fc3 import show_drift_fc3
except:
    pass

mpl.use('Agg')


def get_prim(structure: Structure) -> Structure:
    """
    get a primitive structure
    Args:
        structure: Structure object

    Returns: Structure object

    """
    kpath = HighSymmKpath(structure, symprec=0.01)
    structure = kpath.prim
    return structure


def get_kpath(structure: Structure):
    """
    get high-symmetry points in k-space
    Args:
        structure: Structure Object

    Returns:

    """
    kpath = HighSymmKpath(structure, symprec=0.01)
    kpath_save = kpath.kpath
    labels = copy.deepcopy(kpath_save["path"])
    path = copy.deepcopy(kpath_save["path"])

    for ilabelset, labelset in enumerate(labels):
        for ilabel, label in enumerate(labelset):
            path[ilabelset][ilabel] = kpath_save["kpoints"][label]
    return kpath_save, path


def get_ase_atoms(structure):
    pymatgentoase = AseAtomsAdaptor()
    atoms_new = pymatgentoase.get_atoms(structure)
    return atoms_new


def get_structure(atoms):
    pymatgentoase = AseAtomsAdaptor()
    structure = pymatgentoase.get_structure(atoms)
    return structure


class CompareZPE:
    """
    class to compare thermal properties
    """

    def __init__(self, phoncalc1, phoncalc2):
        self.phoncalc1 = phoncalc1
        self.phoncalc2 = phoncalc2

    def get_rms_zpe(self):
        zpe1 = self.phoncalc1.get_zero_point_energy()
        zpe2 = self.phoncalc2.get_zero_point_energy()
        return np.sqrt((zpe1 - zpe2) ** 2)

    def get_rms_free_energy(self, startt=0, stopt=1000, stept=10):
        list1 = []
        list2 = []
        for t in range(startt, stopt + stept, stept):
            list1.append(self.phoncalc1.get_free_energy(temperature=t))
            list2.append(self.phoncalc2.get_free_energy(temperature=t))
        diff = np.array(list1) - np.array(list2)
        return np.sqrt(np.mean(diff ** 2))


class ComparePhononBS:
    """
    Class that will compare phonon band structure objects from pymatgen
    rms values?
    """

    def __init__(self, bs1, bs2):
        self.bs1 = bs1
        self.bs2 = bs2
        self.bands1 = self.bs1.bands
        self.bands2 = self.bs2.bands

    def rms_kdep(self):
        diff = self.bands1 - self.bands2

        diff = np.transpose(diff)
        kpointdep = [np.sqrt(np.mean(diff[i] ** 2)) for i in range(len(diff))]
        return kpointdep

    def rms_kdep_plot(self, whichkpath=1, filename="rms.eps", format="eps"):
        rms = self.rms_kdep()

        if whichkpath == 1:
            plotter = PhononBSPlotter(bs=self.bs1)
        elif whichkpath == 2:
            plotter = PhononBSPlotter(bs=self.bs2)

        distances = []
        for element in plotter.bs_plot_data()["distances"]:
            distances.extend(element)
        import matplotlib.pyplot as plt
        plt.close("all")
        plt.plot(distances, rms)
        plt.xticks(ticks=plotter.bs_plot_data()["ticks"]["distance"], labels=plotter.bs_plot_data()["ticks"]["label"])
        plt.xlabel("Wave vector")
        plt.ylabel("Phonons RMS (THz)")
        plt.savefig(filename, format=format)

    #TODO: implement test!
    def compare_plot(self, filename="band_comparison.eps",img_format="eps"):
        plotter = PhononBSPlotter(bs=self.bs1)
        plotter2 = PhononBSPlotter(bs=self.bs2)
        new_plotter=plotter.plot_compare(plotter2)
        new_plotter.savefig(filename, format=img_format)
        new_plotter.close()

    def rms_overall(self):
        diff = self.bands1 - self.bands2
        return np.sqrt(np.mean(diff ** 2))
     
    def rms_overall_second_definition(self):
        #makes sure the frequencies are sorted by energy
        band1= np.sort(self.bands1,axis=0)
        band2= np.sort(self.bands2,axis=0)
        diff=band1-band2
        return np.sqrt(np.mean(diff ** 2)) 
        
    def mean_absolute_error(self):
        band1= np.sort(self.bands1,axis=0)
        band2= np.sort(self.bands2,axis=0)
        diff_perc=(np.abs(band1-band2))/band1
        return np.mean(np.abs(diff_perc))
          
   
       
      
       
     



class PhonopyFiniteDisplacements:
    """
    Class with some standard methods that can be reused in PhononsQuippyFiniteDisplacements
    """

    def __init__(self, filename="phonopy.yaml", npoints_band=51, kpoint_density=12000):
        self.filename = filename
        self.npoints_band = npoints_band
        self.kpoint_density = kpoint_density

    def run_all(self):
        self.phonon = self.read_yaml(self.filename)

        self.optimized_structure = get_pmg_structure(self.phonon.get_unitcell())
        self.supercell_matrix = self.phonon.get_supercell_matrix()
        self.phonon_band_structure_pymatgen = self._get_bandstructure_calc(structure=self.optimized_structure,
                                                                           phonon=self.phonon,
                                                                           npoints_band=self.npoints_band)
        self.phonon_dos_pymatgen = self._get_dos_calc(structure=self.optimized_structure, phonon=self.phonon,
                                                      kpoint_density=self.kpoint_density)

    def read_yaml(self, filename):
        return load(filename)

    def _get_primitive_cell(self, structure):
        # returns a primitive cell
        return get_prim(structure)

    def _get_ase_from_pmg(self, structure):
        return get_ase_atoms(structure)

    def _get_pmg_from_ase(self, atoms):
        return get_structure(atoms)

    def _get_bandstructure_calc(self, structure, phonon, npoints_band=51):
        tempfilename = tempfile.gettempprefix() + '.yaml'
        kpath_dict, kpath_concrete = get_kpath(structure)
        qpoints, connections = get_band_qpoints_and_path_connections(kpath_concrete, npoints=npoints_band)
        phonon.run_band_structure(qpoints, path_connections=connections)
        phonon.write_yaml_band_structure(
            filename=tempfilename)
        bs_symm_line = get_ph_bs_symm_line(tempfilename, labels_dict=kpath_dict["kpoints"])
        os.remove(tempfilename)
        return bs_symm_line

    def _get_dos_calc(self, phonon, structure, kpoint_density):
        tempfilename = tempfile.gettempprefix() + '.yaml'
        kpoint = Kpoints.automatic_density(structure=structure, kppa=kpoint_density, force_gamma=True)
        phonon.run_mesh(kpoint.kpts[0])
        phonon.run_total_dos()
        phonon.write_total_dos(filename=tempfilename)
        dos = get_ph_dos(tempfilename)
        os.remove(tempfilename)
        return dos

    def save_plot_band_structure(self, filename, img_format="eps", units="thz", ylim=None):
        plotter = PhononBSPlotter(bs=self.phonon_band_structure_pymatgen)
        plotter.save_plot(filename, img_format=img_format, units=units, ylim=ylim)

    def save_plot_dos(self, filename, img_format="eps", units="thz", label="total dos"):
        plotter2 = PhononDosPlotter(stack=False, sigma=None)
        plotter2.add_dos(label=label, dos=self.phonon_dos_pymatgen)
        plotter2.save_plot(filename=filename, img_format=img_format, units=units)

    def has_imag_modes(self, tol=1e-5):
        return self.phonon_band_structure_pymatgen.has_imaginary_freq(tol=tol)

    def get_zero_point_energy(self):
        return self.phonon_dos_pymatgen.zero_point_energy(self.optimized_structure)

    def get_free_energy(self, temperature):
        return self.phonon_dos_pymatgen.helmholtz_free_energy(t=temperature, structure=self.optimized_structure)


#
# class PhonopyVaspEvaluator(PhonopyFiniteDisplacements):
#     """
#     Class to evaluate a phonon calculation with phonopy and vasp
#     1. needs POSCAR and FORCE_SETS files
#     2. will generate phonopy object
#     3. will calculate bands/dos based on this
#     """
#
#     def __init__(self, poscarpath="POSCAR", forcespath="FORCE_SETS", smat=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#                  path_parameters="phonopy.yaml", save_parameters=True,
#                  npoints_band=51, kpoint_density=12000):
#         self.optimized_structure = Structure.from_file(poscarpath)
#         self.smat = smat
#         self.path_parameters = path_parameters
#         self.save_parameters = save_parameters
#         self.npoints_band = npoints_band
#         self.kpoint_density = kpoint_density
#         self.poscarpath = poscarpath
#         self.forcespath = forcespath
#
#     def run_all(self):
#         self.phonon = self._get_phononobject_phonopy(poscarpath=self.poscarpath, forcespath=self.forcespath,
#                                                      smat=self.smat, save_parameters=self.save_parameters,
#                                                      path=self.path_parameters)
#
#         self.phonon_band_structure_pymatgen = self._get_bandstructure_calc(structure=self.optimized_structure,
#                                                                            phonon=self.phonon,
#                                                                            npoints_band=self.npoints_band)
#         self.phonon_dos_pymatgen = self._get_dos_calc(structure=self.optimized_structure, phonon=self.phonon,
#                                                       kpoint_density=self.kpoint_density)
#
#     def _get_phononobject_phonopy(self, poscarpath, forcespath, smat, save_parameters=False, path="phonopy.yaml"):
#         phonon = load(supercell_matrix=smat, unitcell_filename=poscarpath, force_sets_filename=forcespath,
#                       symprec=10e-5)
#         if save_parameters:
#             phonon.save(path)
#         return phonon


class PhononsQuippyFiniteDisplacements(PhonopyFiniteDisplacements):
    """
    Class to optimize structure and calculate phonon properties
    1. will take a structure object from pymatgen
    2. will calculate primitive cell
    3. will optimize primitive cell with potential and ase -> will be able to export this structure!
    4. will get displacements with phonopy and calculate forces -> will be able to export it as FORCES_SETS
    5. will get kpath for the calculation of the band structure
    6. will calculate dos
    7. will export dos and bandstructure to file
    8. will have the option to read in phonopy stuff to pymatgen -> will be able to export bands
    9. will be able to produce nice graphs and
    10. will calculate everything related to thermal conductivity if I say so
    """

    def __init__(self, structure: Structure, potential, smat=None, path_parameters="phonopy.yaml",
                 npoints_band: int =51, kpoint_density: float=12000, kpoint_density_phono3py=1000, displacementdistance=0.01,
                 set_phonons=True,
                 set_thermal_conductivity: bool = False,
                 displacementdistancephono3py: float = 0.03, work_with_primitive: bool = True,
                 max_distance_third_order=None,
                 temperature_range_kappa=range(50, 1001, 5)):
        """
        Class to optimize structure and calculate phonon properties with (GAP) potentials
        Args:
            structure:
            potential:
            smat:
            path_parameters:
            npoints_band:
            kpoint_density:
            kpoint_density_phono3py:
            displacementdistance:
            set_phonons:
            set_thermal_conductivity:
            displacementdistancephono3py:
            work_with_primitive:
            max_distance_third_order:
            temperature_range_kappa:
        """
        self.initial_structure = structure
        self.potential = potential
        self.smat = smat
        self.path_parameters = path_parameters
        self.npoints_band = npoints_band
        self.kpoint_density = kpoint_density
        self.kpoint_density_phono3py = kpoint_density_phono3py
        self.displacementdistance = displacementdistance
        self.displacementdistancephono3py = displacementdistancephono3py
        self.set_phonons = set_phonons
        self.set_thermal_conductivity = set_thermal_conductivity
        self.work_with_primitive = work_with_primitive
        self.temperature_range_kappa = temperature_range_kappa
        self.max_distance_third_order = max_distance_third_order

    def run_all(self):
        if self.work_with_primitive:
            structure = self._get_primitive_cell(self.initial_structure)
        else:
            structure = self.initial_structure
        self.intial_structure_primitive = structure
        # print(self.initial_structure.lattice)
        self.optimized_structure = self._get_optimized_cell(structure, self.potential)

        if self.set_phonons:
            if self.smat == None:
                self.smat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
                for ilength, length in enumerate(self.optimized_structure.lattice.abc):
                    self.smat[ilength][ilength] = math.ceil(15.0 / length)

            self.energy_optimized_structure = self._get_potential_energy(self.optimized_structure, self.potential)
            self.phonon = self._get_phononobject_phonopy(self.optimized_structure, self.potential, smat=self.smat,
                                                         save_parameters=True, path=self.path_parameters,
                                                         displacement_distance=self.displacementdistance)

            self.phonon_band_structure_pymatgen = self._get_bandstructure_calc(structure=self.optimized_structure,
                                                                               phonon=self.phonon,
                                                                               npoints_band=self.npoints_band)

            self.phonon_dos_pymatgen = self._get_dos_calc(structure=self.optimized_structure, phonon=self.phonon,
                                                          kpoint_density=self.kpoint_density)
        if self.set_thermal_conductivity:
            # will do everything to calculate thermal conductivity
            self.phono3py = self._get_phono3pyobject_phono3py(self.optimized_structure, self.potential,
                                                              kpoint_density=self.kpoint_density_phono3py,
                                                              displacementdistancephono3py=self.displacementdistancephono3py,
                                                              max_distance_third_order=self.max_distance_third_order)
            self.kappa = self._get_thermal_conductivity_matrix(temperatures=self.temperature_range_kappa)
            self.kappa_xx = []
            self.kappa_yy = []
            self.kappa_zz = []
            self.kappa_yz = []
            self.kappa_xz = []
            self.kappa_xy = []
            self.kappa_mean = []
            for value in self.kappa[0]:
                self.kappa_xx.append(value[0])
                self.kappa_yy.append(value[1])
                self.kappa_zz.append(value[2])
                self.kappa_yz.append(value[3])
                self.kappa_xz.append(value[4])
                self.kappa_xy.append(value[5])
                self.kappa_mean.append((value[0] + value[1] + value[2]) / 3.0)

    def kappa_at_temperature(self, temperature=100, whichvalue='xx'):
        # return a  kappa at a certain temperature
        # print(self.temperature_range_kappa)
        # print(self.kappa_xx)
        for itemp, temp in enumerate(self.temperature_range_kappa):
            if temp == temperature:
                if whichvalue == "xx":
                    return self.kappa_xx[itemp]
                elif whichvalue == "yy":
                    return self.kappa_yy[itemp]
                elif whichvalue == "zz":
                    return self.kappa_zz[itemp]
                elif whichvalue == "yz":
                    return self.kappa_yz[itemp]
                elif whichvalue == "xz":
                    return self.kappa_xz[itemp]
                elif whichvalue == "xy":
                    return self.kappa_xy[itemp]

    def calculate_rms(self, exp_data_x=None, exp_data_y=None):
        to_compare = []
        exp_data_y2 = []
        for ix, x in enumerate(exp_data_x):
            for itemp, temp in enumerate(self.temperature_range_kappa):
                if x == temp:
                    to_compare.append(self.kappa_mean[itemp])
                    exp_data_y2.append(exp_data_y[ix])
        diff = np.array(to_compare) - np.array(exp_data_y2)
        return np.sqrt(np.mean(diff ** 2))

    def calculate_rms_xyz(self, exp_data_T=None, exp_data_xx=None, exp_data_yy=None, exp_data_zz=None):
        to_compare = []
        exp_data_y2 = []
        for number in range(0, 3):
            for ix, x in enumerate(exp_data_T):
                for itemp, temp in enumerate(self.temperature_range_kappa):
                    if x == temp:
                        if number == 0:
                            to_compare.append(self.kappa_xx[itemp])
                            exp_data_y2.append(exp_data_xx[ix])
                        elif number == 1:
                            to_compare.append(self.kappa_yy[itemp])
                            exp_data_y2.append(exp_data_yy[ix])
                        elif number == 2:
                            to_compare.append(self.kappa_zz[itemp])
                            exp_data_y2.append(exp_data_zz[ix])

        diff = np.array(to_compare) - np.array(exp_data_y2)
        return np.sqrt(np.mean(diff ** 2))

    def save_kappa_plot(self, mean_exp_data_T=None, mean_exp_data=None, exp_data_xx=None, exp_data_yy=None,
                        exp_data_zz=None, mean_xx_yy_zz=True, xx=True, yy=True, zz=True, yz=True,
                        xz=True, xy=True, filename="kappa.eps", format='eps'):

        import matplotlib.pyplot as plt
        plt.close("all")

        if xx:
            plt.plot(list(self.temperature_range_kappa), self.kappa_xx, label="xx")
        if yy:
            plt.plot(list(self.temperature_range_kappa), self.kappa_yy, label="yy")
        if zz:
            plt.plot(list(self.temperature_range_kappa), self.kappa_zz, label="zz")
        if yz:
            plt.plot(list(self.temperature_range_kappa), self.kappa_yz, label="yz")
        if xz:
            plt.plot(list(self.temperature_range_kappa), self.kappa_xz, label="xz")
        if xy:
            plt.plot(list(self.temperature_range_kappa), self.kappa_xy, label="xy")
        if mean_xx_yy_zz:
            plt.plot(list(self.temperature_range_kappa), self.kappa_mean, label="mean_xx_yy_zz")
        if mean_exp_data_T is not None and mean_exp_data is not None:
            plt.plot(mean_exp_data_T, mean_exp_data, label="benchmark")
        if mean_exp_data_T is not None and exp_data_xx is not None:
            plt.plot(mean_exp_data_T, exp_data_xx, label="benchmark_xx")
        if mean_exp_data_T is not None and exp_data_yy is not None:
            plt.plot(mean_exp_data_T, exp_data_yy, label="benchmark_yy")
        if mean_exp_data_T is not None and exp_data_zz is not None:
            plt.plot(mean_exp_data_T, exp_data_zz, label="benchmark_zz")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Kappa (W/(mK))")
        plt.legend()
        plt.savefig(filename, format=format)

    def _get_thermal_conductivity_matrix(self, temperatures, write_kappa=True):
        self.phono3py.run_thermal_conductivity(temperatures=temperatures, boundary_mfp=1e6, write_kappa=write_kappa)
        return self.phono3py.get_thermal_conductivity().get_kappa()

    def _get_phono3pyobject_phono3py(self, structure, potential, kpoint_density, displacementdistancephono3py,
                                     max_distance_third_order):
        cell = get_phonopy_structure(structure)

        # TODO: get a mesh with a high-enough accuracy
        kpoint = Kpoints.automatic_density(structure=structure, kppa=kpoint_density, force_gamma=True)
        mesh = kpoint.kpts[0]
        # print(mesh)
        phono3py = Phono3py(cell,
                            self.smat,
                            primitive_matrix=[[1, 0., 0.],
                                              [0., 1, 0.],
                                              [0., 0., 1]],
                            mesh=mesh,
                            log_level=1)

        phono3py.generate_displacements(distance=displacementdistancephono3py,
                                        cutoff_pair_distance=max_distance_third_order)
        scells_with_disps = phono3py.get_supercells_with_displacements()
        #for scell in scells_with_disps:
        #    print(scell)

        disp_dataset = phono3py.get_displacement_dataset()
        # print(disp_dataset)
        numatoms = len(scells_with_disps[0].get_scaled_positions())
        dummy_force = np.zeros((numatoms, 3))
        # print(dummy_force)
        set_of_forces = []
        for scell in scells_with_disps:
            if scell is not None:
                cell = Atoms(symbols=scell.get_chemical_symbols(),
                             scaled_positions=scell.get_scaled_positions(),
                             cell=scell.get_cell(),
                             pbc=True)
                cell.set_calculator(potential)
                forces = cell.get_forces()
                # print(forces)

                drift_force = forces.sum(axis=0)
                print(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force))
                # Simple translational invariance
                # TODO check if this is also correct for phono3py
                for force in forces:
                    force -= drift_force / forces.shape[0]
                set_of_forces.append(forces)
            else:
                set_of_forces.append(dummy_force)
        phono3py.produce_fc3(set_of_forces,
                             displacement_dataset=disp_dataset,
                             symmetrize_fc3r=True)

        fc3 = phono3py.get_fc3()
        # fc2 = phono3py.get_fc2()

        show_drift_fc3(fc3)
        # show_drift_force_constants(fc2, name='fc2')
        return phono3py

    def _get_optimized_cell(self, structure, potential):
        # TODO: optimize this cell!
        # will optimize the cell with ase and any potential
        # will also optimize cell
        atoms = self._get_ase_from_pmg(structure)
        atoms.set_calculator(potential)

        # structure optimization
        try:

            sf = UnitCellFilter(atoms)
            dyn = BFGS(sf)
            dyn.run(fmax=1e-5)
            sf = StrainFilter(atoms)
            dyn = BFGS(sf)
            dyn.run(fmax=1e-5)
            dyn = BFGS(atoms)
            dyn.run(fmax=1e-5)
            sf = UnitCellFilter(atoms)
            dyn = BFGS(sf)
            dyn.run(fmax=1e-5)

        except AttributeError:
            warnings.warn("No optimization is performed.")
        return self._get_pmg_from_ase(atoms)

    def _get_potential_energy(self, structure, potential):
        atoms = self._get_ase_from_pmg(structure)
        atoms.set_calculator(potential)
        return atoms.get_potential_energy()

    def _get_phononobject_phonopy(self, structure, potential, smat, save_parameters, path, displacement_distance=0.01):
        cell = get_phonopy_structure(structure)
        phonon = Phonopy(cell, smat, primitive_matrix=[[1.0, 0.0, 0.0],
                                                       [0.0, 1, 0.0],
                                                       [0., 0., 1.0]], factor=VaspToTHz)

        # displacements
        phonon.generate_displacements(distance=displacement_distance)
        print("[Phonopy] Atomic displacements:")
        disps = phonon.get_displacements()
        for d in disps:
            print("[Phonopy] %d %s" % (d[0], d[1:]))

        supercells = phonon.get_supercells_with_displacements()
        # Force calculations by calculator
        set_of_forces = []
        for scell in supercells:
            cell = Atoms(symbols=scell.get_chemical_symbols(),
                         scaled_positions=scell.get_scaled_positions(),
                         cell=scell.get_cell(),
                         pbc=True)
            cell.set_calculator(potential)
            forces = cell.get_forces()
            drift_force = forces.sum(axis=0)
            print(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force))
            # Simple translational invariance
            # TODO: understand this part!
            for force in forces:
                force -= drift_force / forces.shape[0]
            set_of_forces.append(forces)

        phonon.produce_force_constants(forces=set_of_forces)
        if save_parameters:
            phonon.save(path)
        return phonon
