from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
from ovito.data import CutoffNeighborFinder
import statistics
import numpy as np
import matplotlib.pyplot as plt

# Load a particle dataset, apply the modifier, and evaluate pipeline.
pipeline = import_file("/u/vld/magd5247/lammps_gap/ACE-iter5-4/traj-BZS-relaxed-hard-sphere-1500K-50ps-300K-10-13quench.dump")
modifier = CoordinationAnalysisModifier(cutoff=3.1, number_of_bins=200)
pipeline.modifiers.append(modifier)

coord_Zr = []
zr_s_distances = []

def compute_distance_with_pbc(pos1, pos2, cell):
    delta = pos1 - pos2
    delta -= np.round(delta / cell.diagonal()) * cell.diagonal()  # Apply minimum image convention
    return np.linalg.norm(delta)

# Modify loop to process as many timesteps as desired

timesteps = list(range(0, pipeline.source.num_frames, 1))

for frame in timesteps:
    data = pipeline.compute(frame)
    data.cell_.pbc = (True, True, True)
    particle_types = data.particles['Particle Type']
    positions = data.particles['Position']
    cell = data.cell_.matrix[:3, :3]
    
    # Get the particle type IDs for Zr and S
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
    
    neighbor_finder = CutoffNeighborFinder(3.1, data)
    
    zr_coord_numbers = []
    frame_distances = []

    # Loop over all Zr particles
    for i in range(data.particles.count):
        if particle_types[i] == zr_type_id:
            zr_coord = 0
            zr_distances = []
            
            for neigh in neighbor_finder.find(i):
                if particle_types[neigh.index] == s_type_id:
                    zr_coord += 1
                    zr_s_distance = compute_distance_with_pbc(positions[i], positions[neigh.index], cell)
                    zr_distances.append(zr_s_distance)

            # Append to the lists only if neighbors are found
            zr_coord_numbers.append(zr_coord)
            frame_distances.extend(zr_distances)

    # Store coordination and bond length data
    coord_Zr.append(zr_coord_numbers)
    zr_s_distances.append(frame_distances)

# Calculate the average coordination and bond length per frame
avg_coord = [statistics.fmean(coords) if coords else 0 for coords in coord_Zr]
avg_bond_lengths = [statistics.fmean(distances) if distances else 0 for distances in zr_s_distances]

# Print results

print("Average Zr-S coordination number per frame:", avg_coord)
print("Average Zr-S bond lengths per frame:", avg_bond_lengths)

