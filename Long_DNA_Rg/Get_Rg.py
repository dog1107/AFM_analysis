"""
Script to run image analysis on AFM images of bare_DNA and mono-nucleosomes
"""

from tqdm import tqdm
import numpy as np
from skimage.transform import rescale


import sys
import os 

sys.path.append(os.path.join(os.path.expanduser("~"), "Desktop", "classes"))
sys.path.append(os.path.join(os.path.expanduser("~"), "Desktop", "lib"))


import import_custom
import export_custom
import molecule_categorization as cat
import analysis_bare_DNA
import analysis_nucleosome
import analysis_nucleosome_eb
import analysis_ints
import plot_functions as plotting

# Analysis decision making
# Note: Endbound nucleosomes will only be analyzed if both manual filtering and analyze_nucleosomes_eb are True
manual_filtering = False                    # Manually inspect all uncategorized molecules
kwargs = {'analyze_bare_DNA': False,
          'analyze_nucleosomes': False,
          'analyze_nucleosomes_eb': False,  # Analyze endbound nucleosomes - True or False
          'analyze_ints': True}

# Data export decision making
save_close_ups = True                  # Save close-up shots of each molecule in a png - True or False
export_data = True                     # Store the analyzed data in an Excel sheet - True or False
export_select_manually = True          # Select the data to be exported as valid manually
output_file = 'int_export.xlsx'         # Name of the Excel Workbook to store the data in
add_file_name_pars = False              # Adds the image specifications to the excel sheet - proper image name necessary !!!
afm_type = 'jpk'                        # 'jpk' or 'nanoscope'

par_back_1 = 0.1           # Height value for the first background threshold - should be chosen once and kept the same
par_back_2 = 0.2           # Height value for the second background threshold - should be chosen once and kept the same
par_min_area = 500          # Minimum amount of pixels for a molecule to be considered for analysis
par_dna_bp = 2500            # Length of the bare DNA in bp (used to help categorizing molecules as too short and too long)
par_nuc_min_area = 500      # Minimum amount of pixels above nuc_min_height for a molecule to be considered as nucleosome
nuc_min_height = 500       # Height threshold for pixels to be categorized as nucleosome pixels

# Import the desired .ascii file (manual selection)
img_original, file_name, file_path, x_pixels, y_pixels, x_length = import_custom.import_ascii()

# Find the molecules in the image, each molecule is stored as an entry in the molecules list
img_filtered, molecules = import_custom.find_molecules(img_original, x_pixels, y_pixels,
                                                       background_1=par_back_1,
                                                       background_2=par_back_2,
                                                       min_area=par_min_area,
                                                       afm_type=afm_type)

# Create an AFM molecule instance for each individual molecule
par_pixel_size = x_length/y_pixels
afm_molecules = [cat.AFMMolecule(mol, par_dna_bp, par_pixel_size, background_1=par_back_1,
                                 background_2=par_back_2, min_area=par_min_area,
                                 nuc_min_height=nuc_min_height, nuc_min_area=par_nuc_min_area,
                                 categorize=True) for mol in molecules]

# Split the AFM molecules into lists depending on their type
mol_bare_DNA = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Bare DNA']
mol_nucleosome = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Nucleosome']
mol_int = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Int']
mol_trash = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Trash']
print('\nMolecules found:')
print('{} bare DNA strands'.format(len(mol_bare_DNA)))
print('{} nucleosomes'.format(len(mol_nucleosome)))
print('{} ints'.format(len(mol_int)))
print('{} undefined molecules'.format(len(mol_trash)))

if manual_filtering is True:
    # Give possibility to manually separate trashed molecules
    afm_molecules = cat.manual_trash_analysis(afm_molecules)

# Again split the AFM molecules into lists depending on their type after manually helping categorization
mol_bare_DNA = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Bare DNA']
mol_nucleosome = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Nucleosome']
mol_nucleosome_eb = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Nucleosome endbound']
mol_int = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Int']
mol_trash = [mol for mol in afm_molecules if mol.mol_pars['type'] == 'Trash']
if manual_filtering is True:
    print('\nMolecules found after manual separation:')
    print('{} bare DNA strands'.format(len(mol_bare_DNA)))
    print('{} nucleosomes'.format(len(mol_nucleosome)))
    print('{} nucleosomes endbound'.format(len(mol_nucleosome_eb)))
    print('{} ints'.format(len(mol_int)))
    print('{} undefined molecules'.format(len(mol_trash)))

results_final = {}
# Analysis of Bare DNA
if kwargs['analyze_bare_DNA'] is True:
    print('\nAnalyzing bare DNA:')
    analyzed_bare_DNA = [analysis_bare_DNA.BareDNA([mol.mol_original, mol.mol_filtered, mol.mol_bbox],
                                                   par_dna_bp, par_pixel_size, par_back_1, par_back_2,
                                                   par_min_area, nuc_min_height, mol.mol_pars)
                         for mol in tqdm(mol_bare_DNA)]
    results_final.update({'analyzed_bare_DNA': analyzed_bare_DNA})


    # Write the avg. length values
    bare_lengths = [mol.results['length_avg'] for mol in analyzed_bare_DNA if mol.results['length_avg'] is not False]
    print('\n# BARE DNA RESULTS #')
    print('Length was measured for {} molecules'.format(len(bare_lengths)))
    print('Average length: {} nm. (Sigma = {} nm)'.format(
        np.round(np.average(bare_lengths), decimals=1),
        np.round(np.std(bare_lengths), decimals=1)))
    print('Topmost 50%: {} nm. (Sigma = {} nm)'.format(
        np.round(np.average(bare_lengths[0:int(len(bare_lengths)/2)])),
        np.round(np.std(bare_lengths[0:int(len(bare_lengths)/2)]))))
    print('Botmost 50%: {} nm. (Sigma = {} nm)'.format(
        np.round(np.average(bare_lengths[int(len(bare_lengths)/2):len(bare_lengths)]), decimals=1),
        np.round(np.std(bare_lengths[int(len(bare_lengths)/2):len(bare_lengths)])), decimals=1))

# Analysis of normal mono-nucleosomes
if kwargs['analyze_nucleosomes'] is True:
    print('\nAnalyzing mono-nucleosomes:')
    analyzed_nucleosomes = [analysis_nucleosome.Nucleosome([mol.mol_original, mol.mol_filtered, mol.mol_bbox],
                                                           par_dna_bp, par_pixel_size, par_back_1, par_back_2,
                                                           par_min_area, nuc_min_height, par_nuc_min_area, mol.mol_pars)
                            for mol in tqdm(mol_nucleosome)]
    results_final.update({'analyzed_nucleosomes': analyzed_nucleosomes})

# Analysis of endbound mono-nucleosomes
if kwargs['analyze_nucleosomes_eb'] is True:
    print('\nAnalyzing endbound mono-nucleosomes:')
    analyzed_nucleosomes_eb = [analysis_nucleosome_eb.NucleosomeEB([mol.mol_original, mol.mol_filtered, mol.mol_bbox],
                                                                   par_dna_bp, par_pixel_size, par_back_1, par_back_2,
                                                                   par_min_area, nuc_min_height, par_nuc_min_area,
                                                                   mol.mol_pars)
                               for mol in tqdm(mol_nucleosome_eb)]
    results_final.update({'analyzed_nucleosomes_eb': analyzed_nucleosomes_eb})

# Analysis of Ints
if kwargs['analyze_ints'] is True:
    print('\nAnalyzing Ints:')
    analyzed_ints = [analysis_ints.Inta([mol.mol_original, mol.mol_filtered, mol.mol_bbox],
                                        par_dna_bp, par_pixel_size, par_back_1, par_back_2,
                                        par_min_area, nuc_min_height, par_nuc_min_area,
                                        mol.mol_pars)
                     for mol in tqdm(mol_int)]
    results_final.update({'analyzed_ints': analyzed_ints})

# Export the data
if export_data is True:
    if export_select_manually is True:
        results_final = export_custom.manual_export_selection(results_final)

    export_custom.export_to_excel(output_file,
                                  results_final,
                                  file_path,
                                  add_file_name_pars=add_file_name_pars,
                                  **kwargs)

# Plot the results
plotting.plot_overview_image(img_filtered,
                             file_name,
                             results_final,
                             **kwargs)

if save_close_ups is True:
    results_final.update({'mol_trash': mol_trash})
    plotting.plot_save_close_ups(results_final,
                                 file_name,
                                 **kwargs,
                                 plot_trash=True)
