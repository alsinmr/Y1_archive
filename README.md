# Y1_archive
Archive of scripts used for Y1 MD simulation analysis


This archive provides the python scripts used to create figures in the paper:

"Fine-Tuning of the Activity of the Neuropeptide Y1 G Protein-Coupled Receptor by the Tryptophan6.48 ”Toggle Switch”"

by

Matthias Voitel, Maik Pankonin, Alexander Vogel, Anette Kaiser, Daniel Huster, Peter Hildebrand, and Albert A. Smith

This is only the Python code, in order to show how these analyses are performed. In order to run the code, it is necessary to request the MD trajectories. 

Note that Python files are names as FigXX or SI_FigXX, corresponding to which plots they produce. Other files either provide objects that need to be imported to create the figures, or analyses that need to be run before figure creation.

A copy of pyDR has been provided in this archive, which includes our code for performing detector analysis and principal component analysis.

There is NO INSTALLATION required for the provided code. Just place everything in a folder, navigate there, and run with python3. However, python3 and the following modules must be installed from other sources to run pyDR (these are the tested versions, although other versions may work).

Python v. 3.7.3, numpy v. 1.19.2, scipy v. 1.5.2, MDAnalysis v. 1.0.0, matplotlib v. 3.4.2

Additionally, ChimeraX must be installed and its location provided to pyDR. 

All files are copyrighted under the GNU General Public License. A copy of the license has been provided in the file LICENSE

Copyright 2025 Albert Smith-Penzel

This study was funded through the Deutsche Forschungsgemeinde (DFG) through CRC 1423, project number 421152132, subprojects A04 (DH) and C01 (PWH) and DFG grant 45014882 (AAS). 
