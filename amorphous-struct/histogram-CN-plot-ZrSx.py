from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager

# Read in the structures of the amorphous and crystalline phases (LAMMPS data files)
amorphous_file = "/u/vld/magd5247/lammps_gap/ACE-iter5-4/no-high-E/BZS-2500K-80ps-10^14quench-1500K-15ps-300K-10^15quench.data"
crystalline_file = "/u/vld/magd5247/lammps_gap/ACE-iter5-4/no-high-E/BZS-cif-MP-10000atoms.data"  # Change this to your actual file

# Load structures into OVITO pipeline
amorphous_pipeline = import_file(amorphous_file)
crystalline_pipeline = import_file(crystalline_file)

# Define the cutoff for the bond (angstroms)
zr_s_cutoff = 3.1  # Adjust if necessary

# Compute data for both structures
amorphous_data = amorphous_pipeline.compute()
crystalline_data = crystalline_pipeline.compute()

# Enable periodic boundary conditions
amorphous_data.cell_.pbc = (True, True, True)
crystalline_data.cell_.pbc = (True, True, True)

# Identify particle types
particle_types_amorphous = amorphous_data.particles['Particle Type']
particle_types_crystalline = crystalline_data.particles['Particle Type']

# Extract type IDs
type_ids = {}
for t in amorphous_data.particles.particle_types.types:
    type_ids[t.name] = t.id

if "Zr" not in type_ids:
    raise ValueError("Zr type not found in the dataset.")

zr_type_id = type_ids["Zr"]

# Function to compute coordination percentages for Zr
def compute_zr_coordination_percentages(data, zr_type_id, cutoff):
    neighbor_finder = CutoffNeighborFinder(cutoff, data)
    zr_coordination = []

    for i in range(data.particles.count):
        if data.particles['Particle Type'][i] == zr_type_id:
            coord_number = sum(1 for _ in neighbor_finder.find(i))  # Count neighbors
            zr_coordination.append(coord_number)

    total_zr = len(zr_coordination)  # Total number of Zr atoms
    if total_zr == 0:
        return {5: 0, 6: 0, 7: 0}  # Avoid division by zero

    # Compute percentage of Zr atoms with 5, 6, and 7 coordination
    percentages = {c: (zr_coordination.count(c) / total_zr) * 100 for c in [5, 6, 7]}
    return percentages

# Get coordination statistics for both phases
amorphous_percentages = compute_zr_coordination_percentages(amorphous_data, zr_type_id, zr_s_cutoff)
crystalline_percentages = compute_zr_coordination_percentages(crystalline_data, zr_type_id, zr_s_cutoff)

# Prepare data for plotting
categories = [5, 6, 7]  # Coordination numbers of interest in the case of ZrSx polyhedra
amorphous_values = [amorphous_percentages.get(c, 0) for c in categories]
crystalline_values = [crystalline_percentages.get(c, 0) for c in categories]

#optional: adjust font size and style of figure

#rcParams['font.family'] = 'DeJavu Serif'
#rcParams['font.serif'] = ['Arial']
#plt.rcParams.update({
 #   'axes.labelsize': 8,
 #   'axes.facecolor': 'w',
 #   'xtick.labelsize': 6,
 #   'ytick.labelsize': 6,
 #   'legend.fontsize': 6,
 #   'legend.fancybox': False,
 #  'legend.edgecolor': 'k',
 #  'legend.borderaxespad': 1.5,
 #  'lines.linewidth': 1.0,
 #  'axes.linewidth': 0.65,
 #  'xtick.direction': 'out',
 #   'xtick.major.size': 5,
 #  'xtick.major.width': 0.65,
 #   'xtick.minor.size': 3,
 #   'xtick.minor.width': 0.65,
 #   'ytick.direction': 'out',
 #   'ytick.major.size': 5,
 #   'ytick.major.width': 0.65,
#  'ytick.minor.size': 3,
 #   'ytick.minor.width': 0.65,
#})
#plt.rcParams['mathtext.fontset'] = 'custom'
#plt.rcParams['mathtext.rm'] = 'Arial'
#plt.rcParams['mathtext.it'] = 'Arial:italic'
#plt.rcParams['mathtext.bf'] = 'Arial:bold'

# Plot histogram comparing amorphous vs crystalline Zr coordination (percentages)
bar_width = 0.4
x = np.arange(len(categories))
fig, ax= plt.subplots(figsize=(8, 6))
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2) 
plt.bar(x - bar_width/2, amorphous_values, width=bar_width, label="Amorphous", color="orange", alpha=0.7)
plt.bar(x + bar_width/2, crystalline_values, width=bar_width, label="Crystalline", color="red", alpha=0.5)

# Labels and formatting
plt.xticks(x, [f"CN={c}" for c in categories], fontsize=20)
plt.yticks(fontsize=18)
plt.ylabel("Percentage (%)", fontsize=20)
plt.legend(fontsize=20, bbox_to_anchor=(0.5, 1.07))  


plt.savefig("histogram-Ba-coordination.png")
plt.show()

