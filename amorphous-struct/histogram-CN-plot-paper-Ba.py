from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager as font_manager

# Load Arial font
font_path = '/u/vld/magd5247/software/miniconda3_test/envs/ovito_env/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path, size=8)


# File paths for amorphous and crystalline structures
amorphous_file = "/u/vld/magd5247/lammps_gap/ACE-iter5-4/BZS-10000atoms-10^13quench-correct.data"
crystalline_file = "/u/vld/magd5247/lammps_gap/ACE-iter5-4/no-high-E/BZS-cif-MP-10000atoms.data"  # Change this to your actual file

# Load structures
amorphous_pipeline = import_file(amorphous_file)
crystalline_pipeline = import_file(crystalline_file)

# Define Zr-S cutoff (angstroms)
zr_s_cutoff = 3.8  # Adjust if necessary

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

if "Ba" not in type_ids:
    raise ValueError("Ba type not found in the dataset.")

zr_type_id = type_ids["Ba"]

# Function to compute coordination percentages for Zr
def compute_zr_coordination_percentages(data, zr_type_id, cutoff):
    neighbor_finder = CutoffNeighborFinder(cutoff, data)
    zr_coordination = []

    particle_types = data.particles['Particle Type']

    # Identify S type
    s_type_id = None
    for t in data.particles.particle_types.types:
        if t.name == 'S':
            s_type_id = t.id
            break

    if s_type_id is None:
        raise ValueError("S type not found in the dataset.")

    for i in range(data.particles.count):
        if particle_types[i] == zr_type_id:
            coord_number = sum(1 for neigh in neighbor_finder.find(i) if particle_types[neigh.index] == s_type_id)
            zr_coordination.append(coord_number)

    total_zr = len(zr_coordination)
    if total_zr == 0:
        return {c: 0 for c in range(2, 13)}

    percentages = {c: round((zr_coordination.count(c) / total_zr) * 100, 1) for c in range(5, 10)}
    return percentages


# Get coordination statistics for both phases
amorphous_percentages = compute_zr_coordination_percentages(amorphous_data, zr_type_id, zr_s_cutoff)
crystalline_percentages = compute_zr_coordination_percentages(crystalline_data, zr_type_id, zr_s_cutoff)

print(amorphous_percentages)

# Prepare data for plotting
categories = [5, 6, 7, 8, 9]  # Coordination numbers of interest
amorphous_values = [amorphous_percentages.get(c, 0) for c in categories]
crystalline_values = [crystalline_percentages.get(c, 0) for c in categories]


#rcParams['font.family'] = 'DeJavu Serif'
#rcParams['font.serif'] = ['Arial']
#plt.rcParams.update({
 #   'axes.labelsize': 8,
#    'axes.facecolor': 'w',
#    'xtick.labelsize': 6,
#    'ytick.labelsize': 6,
#    'legend.fontsize': 6,
#    'legend.fancybox': False,
#    'legend.edgecolor': 'k',
#    'legend.borderaxespad': 1.5,
#    'lines.linewidth': 1.0,
#   'axes.linewidth': 0.65,
 #   'xtick.direction': 'out',
 #   'xtick.major.size': 5,
 #   'xtick.major.width': 0.65,
#    'xtick.minor.size': 3,
#   'xtick.minor.width': 0.65,
#    'ytick.direction': 'out',
#    'ytick.major.size': 5,
#    'ytick.major.width': 0.65,
#    'ytick.minor.size': 3,
#    'ytick.minor.width': 0.65,
#})
#plt.rcParams['mathtext.fontset'] = 'custom'
#plt.rcParams['mathtext.rm'] = 'Arial'
#plt.rcParams['mathtext.it'] = 'Arial:italic'
#plt.rcParams['mathtext.bf'] = 'Arial:bold'

# Plot histogram comparing amorphous vs crystalline Zr coordination (percentages)
bar_width = 0.4
x = np.arange(len(categories))
fig, ax = plt.subplots(figsize=(2.5, 1.875)) 
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(0.65) 
plt.bar(x - bar_width/2, amorphous_values, width=bar_width, label="Amorphous", color="orange", alpha=0.7)
plt.bar(x + bar_width/2, crystalline_values, width=bar_width, label="Crystalline", color="red", alpha=0.5)

# Labels and formatting
plt.xticks(x, [f"{c}" for c in categories], fontproperties=prop, fontsize=6)
plt.xlabel("Coordination Number", fontproperties=prop, fontsize=8,)
ax.tick_params(axis='both', direction='out', width=0.65)
plt.yticks(np.arange(0, 101, 20), fontproperties=prop, fontsize=6)  # Ticks at 0, 20, 40, 60, 80, 100%
plt.ylabel("Percentage",fontproperties=prop, fontsize=8)

plt.legend(
    prop=prop, frameon=False, handletextpad=0.3, borderpad=0.1, 
    loc='upper left', bbox_to_anchor=(0.02, 0.98)  # Fine-tune positioning
)


# Adjust subplot layout to push the plot downward
 # Leaves space at the top
plt.tight_layout()

plt.savefig("histogram-Ba-coordination-10^14quench.png", dpi=600)
plt.show()

