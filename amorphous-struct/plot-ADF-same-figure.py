import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Read the text files for the crystalline and amorphous structure generated using OVIT-ADF.py
crystal_file = '/u/vld/magd5247/Crystal-structures/python-files/bond_angles-amorphous-BaS-bonds.txt'
amorphous_file = '/u/vld/magd5247/Crystal-structures/python-files/bond_angles-crystal-BaS-bonds-3p8cutoff.txt'



# Function to read and normalize angle data
def read_norm_data(file_name):
    angles = []
    counts = []
    with open(file_name, 'r') as file:
        for line in file:
            # Skip comments
            if line.startswith('#'):
                continue
            # Split the line into angle and count
            parts = line.split()
            if len(parts) == 2:
                angles.append(float(parts[0]))
                counts.append(float(parts[1]))
    angles = np.array(angles)
    counts = np.array(counts)
    # Normalize counts
    counts /= counts.max() 
    return angles, counts

font = {'size'   : 8}

plt.rc('font', **font)

# Read and normalize data for crystal and amorphous phases
crystal_angles, crystal_counts = read_norm_data(crystal_file)
amorphous_angles, amorphous_counts = read_norm_data(amorphous_file)

# Start plotting ADF 

fig, ax = plt.subplots(figsize=(8, 6))
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2) 

# crystal phase
ax.plot(crystal_angles, crystal_counts, label='Crystalline', color='red', linewidth=1.5)
ax.fill_between(crystal_angles, crystal_counts, color='red', alpha=0.4)


# amorphous phase
ax.plot(amorphous_angles, amorphous_counts, label='Amorphous', color='red', linewidth=1.5, alpha=0.7)
ax.fill_between(amorphous_angles, amorphous_counts, color='red', alpha=0.4)

#set figure font settings for publication
rcParams['font.family'] = 'DeJavu Serif'
rcParams['font.serif'] = ['Arial']
plt.rcParams.update({
    'axes.labelsize': 8,
    'axes.facecolor': 'w',
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'legend.fancybox': False,
    'legend.edgecolor': 'k',
    'legend.borderaxespad': 1.5,
    'lines.linewidth': 1.0,
    'axes.linewidth': 0.65,
    'xtick.direction': 'out',
    'xtick.major.size': 5,
    'xtick.major.width': 0.65,
    'xtick.minor.size': 3,
    'xtick.minor.width': 0.65,
    'ytick.direction': 'out',
    'ytick.major.size': 5,
    'ytick.major.width': 0.65,
    'ytick.minor.size': 3,
    'ytick.minor.width': 0.65,
})
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.bf'] = 'Arial:bold'


plt.xlabel('Bond Angle (°)')
plt.xlim(0, 180)
plt.ylabel('Angle Distribution (a.u.)')
plt.yticks([], [])
plt.xticks(np.arange(0, 225, 45))

plt.legend(fontsize=22, loc='upper right', bbox_to_anchor=(1.05, 1.08))
plt.tight_layout()

plt.savefig('ADF-Crystal-vs-Amorphous-BaS-bonds.png')
plt.show()
