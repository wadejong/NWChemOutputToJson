# NWChemOutputToJson

Python files for reading NWChem output and converting to Json

Use case:

python NWChemJsonConversion.py [noOrbitals] <one or more output files>   

(noOrbitals is an optional keyword supressing storing orbitals in the JSON file. 
 Note, the orbital coefficients are read from the output, and by default only
 large coefficients are printed. Make sure to adjust the print level to also
 print smaller coefficients, to avoid getting funky looking orbitals. The
 alternative is to use the NWChem ready version writing JSON files, which uses
 all coefficients.)
 
 If you find any problems, please let me know by sending the output file and I'll
 fix the Python code.
