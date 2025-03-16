"""
Class to analyze categorized ints
"""

import copy
import math

import numpy as np
from skimage import morphology

from molecule_categorization import AFMMolecule
import analysis_functions as analysis


class Inta(AFMMolecule):

    def __init__(self, mol, dna_bp, pixel_size, background_1, background_2, min_area, nuc_min_height, nuc_min_area, mol_pars):
        super().__init__(mol, dna_bp, pixel_size, background_1, background_2, min_area, nuc_min_height, nuc_min_area)
        # Copy the variable, otherwise they are also changed in the AFMMolecule instances
        self.mol_pars = copy.deepcopy(mol_pars)
        self.results = {}
        self.results.update({'position_row': self.mol_bbox[0],
                             'position_col': self.mol_bbox[1],
                             'failed': False})

        if self.results['failed'] is False:
            self.calculate_rog()

        if math.isnan(self.results['com_core_r']):
            self.results.update({'com_core_r': 0,
                                 'com_core_c': 0,
                                 'center_of_mass_core': (0, 0)})

        if self.results['radius_of_gyration_core'] != 0:
            self.ellipsoid_fit()
        else:
            self.results.update({'ellipsoid_coeff': 0,
                                 'ellipsoid_width_a': 0,
                                 'ellipsoid_width_b': 0,
                                 'ellipsoid_height': 0,
                                 'ellipsoid_angle': 0,
                                 'ellipsoid_var_matrix': 0,
                                 'mol_nuc_ellipsoid_cut': 0,
                                 'ellipsoid_pixels': 0})

        if self.results['ellipsoid_width_a'] != 0:
            self.volume()
        else:
            self.results.update({'total_volume': 0,
                                 'int_volume': 0,
                                 'int_max_height': 0,
                                 'int_max_height_avg': 0})

    def calculate_rog(self):

        # Calculate radius of gyration and center of mass of the whole molecule
        rog, com = analysis.radius_of_gyration(self.mol_filtered)

        # Calculate the same for the nucleosome core particle
        mol_labelled = copy.deepcopy(self.mol_filtered)
        mol_labelled[mol_labelled < self.nuc_min_height] = 0
        mol_labelled[mol_labelled >= self.nuc_min_height] = 1
        mol_labelled = morphology.label(mol_labelled, connectivity=2)
        if np.amax(mol_labelled) != 1:
            mol_labelled = morphology.remove_small_objects(mol_labelled, self.mol_pars['max_area_over_height'])
        mol_core_part = copy.deepcopy(self.mol_filtered)
        mol_core_part[mol_labelled == 0] = 0

        rog_core, com_core = analysis.radius_of_gyration(mol_core_part)

        self.results.update({'radius_of_gyration': rog * self.pixel_size,
                             'radius_of_gyration_core': rog_core * self.pixel_size,
                             'center_of_mass': com,
                             'center_of_mass_core': com_core,
                             'com_r': com[0],
                             'com_c': com[1],
                             'com_core_r': com_core[0],
                             'com_core_c': com_core[1]})

        return self

    def ellipsoid_fit(self):

        grid_size = 10
        if int(2 * self.results['radius_of_gyration_core'] / self.pixel_size) > 10:
            grid_size = int(2 * self.results['radius_of_gyration_core'] / self.pixel_size)
        start = [self.results['radius_of_gyration_core'] / self.pixel_size,
                 self.results['radius_of_gyration_core'] / self.pixel_size,
                 4.]
        coeff, var_matrix, mol_nuc_ellipsoid_cut, ellipsoid_pixels = analysis.\
            nuc_core_ellipsoid_fit(self.mol_filtered,
                                   self.results['center_of_mass_core'],
                                   grid_size=grid_size,
                                   start=start)
        try:
            self.results.update({'ellipsoid_coeff': coeff,
                                 'ellipsoid_width_a': coeff[2] * self.pixel_size,
                                 'ellipsoid_width_b': coeff[3] * self.pixel_size,
                                 'ellipsoid_height': coeff[4],
                                 'ellipsoid_angle': coeff[5],
                                 'ellipsoid_var_matrix': var_matrix,
                                 'mol_nuc_ellipsoid_cut': mol_nuc_ellipsoid_cut,
                                 'ellipsoid_pixels': ellipsoid_pixels})
        except:
            self.results.update({'failed': True,
                                 'failed_reason': 'Ellipsoid fit'})

        return self

    def volume(self):
        """ Calculate the integrase volume and height based on the pixels within the fitted ground ellipse """

        ell_coeff = self.results['ellipsoid_coeff']
        mol_pixel_locs = np.where(self.mol_filtered != 0)

        # Pre-filter potential nucleosome pixels by selecting the ones within the range of the largest ellipse axis
        mol_pixels = [np.array([r, c]) for r, c in zip(mol_pixel_locs[0], mol_pixel_locs[1])
                      if np.linalg.norm((np.array([r, c]) - ell_coeff[0:2])) <= np.amax(ell_coeff[2:4])]

        # Find all pixels whose centers lie within the ground ellipse of the fitted ellipsoid
        inner_pixels = [pixel for pixel in mol_pixels
                        if np.linalg.norm((analysis.ellipse_arm_pixel([pixel], ell_coeff, z_h=0) - ell_coeff[0:2]))
                        >= np.linalg.norm((pixel - ell_coeff[0:2]))]

        pixel_heights = [self.mol_filtered[r, c] for r, c in inner_pixels]
        self.results.update({'total_volume': sum(self.mol_filtered[self.mol_filtered != 0]) * self.pixel_size**2,
                             'int_volume': sum(pixel_heights) * self.pixel_size**2,
                             'int_max_height': max(pixel_heights),
                             'int_max_height_avg': np.mean(sorted(pixel_heights)[-5:])})

        return self
