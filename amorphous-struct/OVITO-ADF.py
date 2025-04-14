from ovito.io import import_file, export_file
from ovito.modifiers import CreateBondsModifier, BondAnalysisModifier
import numpy as np
import matplotlib.pyplot as plt
from CalculateBondAnglesAtParticle import CalculateBondAnglesAtParticle


# Load LAMMPS data structure with style=bonds, choosing an appropriate cutoff for each bond type.
pipeline = import_file('/u/vld/magd5247/lammps_gap/ACE-iter5-4/no-high-E/BZS-crystal-BaS-bonds-3p8cutoff.data')

pipeline.modifiers.append(BondAnalysisModifier(bins = 500))

# Export bond angle distribution to an output text file.
export_file(pipeline, 'bond_angles-crystal-BaS-bonds-3p8cutoff.txt', 'txt/table', key='bond-angle-distr')

# Convert bond length histogram to a NumPy array and print it to the terminal.
data = pipeline.compute()
print(data.tables['bond-length-distr'].xy())
