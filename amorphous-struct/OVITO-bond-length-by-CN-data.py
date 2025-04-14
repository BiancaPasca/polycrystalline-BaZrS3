from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
from ovito.data import CutoffNeighborFinder
import statistics
import numpy as np

# Load a LAMMPS .data file and apply the CoordinationAnalysisModifier to compute coordination numbers.
pipeline = import_file("/u/vld/magd5247/lammps_gap/ACE-iter5-4/no-high-E/BZS-amorphous-10000atoms-10^13quench.data")

modifier = CoordinationAnalysisModifier(cutoff=3.1, number_of_bins=200)
pipeline.modifiers.append(modifier)

# Helper function to compute distance considering PBC
def compute_distance_with_pbc(pos1, pos2, cell):
    delta = pos1 - pos2
    delta -= np.round(delta / cell.diagonal()) * cell.diagonal()  # Apply minimum image convention
    return np.linalg.norm(delta)

# Compute data for the static frame
data = pipeline.compute()
data.cell_.pbc = (True, True, True)
particle_types = data.particles['Particle Type']
positions = data.particles['Position']
cell = data.cell_.matrix[:3, :3] 

# Identify particle type IDs for Zr and S
zr_type_id = None
s_type_id = None
for t in data.particles.particle_types.types:
    if t.name == 'Zr':
        zr_type_id = t.id
    elif t.name == 'S':
        s_type_id = t.id

# Ensure both Zr and S types are found
if zr_type_id is None or s_type_id is None:
    raise ValueError("Zr or S type not found in the dataset.")

# Initialize variables for categorizing bond lengths
bond_lengths_by_coord = { "4": [], "5": [], "6": [], "7": [] }
coordination_counts = { "4": 0, "5": 0, "6": 0, "7" : 0}
total_zr_atoms = 0

# Create a neighbor finder object for the current frame
neighbor_finder = CutoffNeighborFinder(3.1, data)

# Loop through all particles to calculate coordination numbers and bond lengths
for i in range(data.particles.count):
    if particle_types[i] == zr_type_id:
        total_zr_atoms += 1
        zr_coord = 0
        zr_bond_lengths = []
        for neigh in neighbor_finder.find(i):
            if particle_types[neigh.index] == s_type_id:
                zr_coord += 1
                zr_s_distance = compute_distance_with_pbc(positions[i], positions[neigh.index], cell)
                zr_bond_lengths.append(zr_s_distance)
        
        # Categorize the bond lengths and coordination numbers
        if zr_coord == 4:
            bond_lengths_by_coord["4"].extend(zr_bond_lengths)
            coordination_counts["4"] += 1
        elif zr_coord == 5:
            bond_lengths_by_coord["5"].extend(zr_bond_lengths)
            coordination_counts["5"] += 1
        if zr_coord == 6:
            bond_lengths_by_coord["6"].extend(zr_bond_lengths)
            coordination_counts["6"] += 1
        elif zr_coord == 7:
            bond_lengths_by_coord["7"].extend(zr_bond_lengths)
            coordination_counts["7"] += 1

# Calculate the average bond lengths for each category
avg_bond_lengths_by_coord = {
    key: (statistics.fmean(lengths) if lengths else 0)
    for key, lengths in bond_lengths_by_coord.items()
}

# Calculate the percentages of coordination categories
coordination_percentages = {
    key: (count / total_zr_atoms) * 100 for key, count in coordination_counts.items()
}

# Output the results
print("Average Zr-S bond lengths by coordination number:")
for coord_category, avg_length in avg_bond_lengths_by_coord.items():
    print(f"Coordination {coord_category}: {avg_length:.3f}")

print("\nPercentages of Zr atoms by coordination number:")
for coord_category, percentage in coordination_percentages.items():
    print(f"Coordination {coord_category}: {percentage:.2f}%")
