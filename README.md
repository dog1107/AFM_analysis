# AFM analysis
Code used for analyzing AFM images of DNA. 

This code is not mine.


The Short_DNA folder is only suitable for DNA with smaller base pair lengths up to 1000. It is not capable of tracing the contour length of longer DNA. It can be used to measure both end to end distance and radius of gyration. It is also necessary to download and use Environment_AFM_analysis, which uses a downgraded version of spyder because the code is written with outdated commands. The environment can be uploaded to anaconda and spyder opened from there. It is also important to adjust the DNA/nucleosome length and be sure that a suitable threshold is used-- and that this threshold is kept constant across all measurements. Measurements are saved to an excel file. 

The Long_DNA_Rg folder can be used for longer length DNA, but only computes the radius of gyration. It does not require the environment used for the other code. For instructions for use, see the readme document in the folder containing the code.


All images must be in ascii format when uploaded for analysis
