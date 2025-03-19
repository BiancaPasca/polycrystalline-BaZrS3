import numpy as np
from ase.io import read
from quippy.descriptors import Descriptor
import matplotlib.pyplot as plt
from umap import UMAP

# Load structures from extxyz file
structures = read("/u/vld/magd5247/lammps_gap/ACE-iter5-4/no-high-E/plus-GB/ACE-iter5-4-GB/train_with_sigma-E-under-1.extxyz", index=":")

# Define SOAP descriptor parameters
soap_descriptor = Descriptor("soap l_max=8 n_max=8 atom_sigma=0.5 cutoff=5 n_Z=3 Z={56 40 16} n_species=3 species_Z={56 40 16} average=T")

# Map configuration types (config_type in the exyxyz) to iteration labels
config_type_map = {
    "De novo exploration": ["bulk", "dimer"], 
    "High-T MD": ["MD"],
    "Crystalline structures": ["crys"],
    "NVT melting": ["4000-1p10vol", "4500-1p05vol", "4500-1p10vol", "4000-1p05vol", "3500-1p20vol", "4500-1p20vol", "high-T-1p10vol", "4000-1p20vol", "high-T-1p05vol", "high-T"],
    "Melt-quench NPT": ["npt", "10^12-quench", "npt-10^15", "iter4-3", "npt-2p5x10^14", "npt-5x10^14", "npt-10^14"],
    "Crystalline-amorphous interfaces": ["Crys-amorph", "crys-amorph", "crys-amorph-iter5-2", "crys-amorph-iter5-3", "crys-amorph-iter5-4", "iter5-5", "iter5-5-extra"]
}

# Function to get color label based on config_type
def get_iter_label(config_type):
    for iter_label, types in config_type_map.items():
        if config_type in types:
            return iter_label
    return None  

# Compute SOAP descriptors for each structure and map to color labels
soap_vectors = []
labels = []
for structure in structures:
    soap_vector = soap_descriptor.calc(structure)['data']
    soap_vectors.append(soap_vector)
    config_type = structure.info.get("config_type")
    labels.append(get_iter_label(config_type))

# Filter out structures with an unknown label, if they exist, to prevent errors (make sure you include all labels in config_type_map)
soap_vectors = np.vstack(soap_vectors)
valid_indices = []
for i, label in enumerate(labels):
    if label is None:
        print(f"Unknown label found: {structures[i].info.get('config_type')}")
    else:
        valid_indices.append(i)
soap_vectors = soap_vectors[valid_indices]
labels = [labels[i] for i in valid_indices]

# Compute the similarity kernel as a dot product of the SOAP vectors
similarity_kernel = np.dot(soap_vectors, soap_vectors.T)

# Reduce to 2D using UMAP (compute once for all data points) - change the parameters as desired to see how this affects the plot
umap_reducer = UMAP(metric='precomputed', min_dist=0.5, n_neighbors=20, n_components=2, random_state=42)
reduced_data_umap = umap_reducer.fit_transform(similarity_kernel)

# Map each label to a color
label_colors = {
    "De novo exploration": "lightskyblue",
    "High-T MD": "steelblue",
    "Crystalline structures": "thistle",
    "NVT melting": "orchid",
    "Melt-quench NPT": "palevioletred",
    "Crystalline-amorphous interfaces": "mediumpurple"
}
legend_labels = {v: k for k, v in label_colors.items()}

# Create incremental plots
selected_labels = ["De novo exploration"]  # Start with the first label
for iter_label in config_type_map.keys():
    selected_labels.append(iter_label)  # Add new label for this iteration

    # Filter points based on selected labels
    mask = [label in selected_labels for label in labels]
    filtered_data = reduced_data_umap[mask]
    filtered_colors = [label_colors[label] for label in labels if label in selected_labels]

    # Plot
    plt.figure(figsize=(10, 8))
    plt.scatter(filtered_data[:, 0], filtered_data[:, 1], c=filtered_colors, rasterized=True)
    plt.xlim(-4, 8)
    plt.ylim(-2, 10)
    plt.xticks([], [])
    plt.yticks([], [])
    
    #uncomment these lines if you want to add labels and title
    #plt.xlabel("UMAP Component 1")
    #plt.ylabel("UMAP Component 2")
    #plt.title(f"UMAP with {', '.join(selected_labels)}")
    
    # Create legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=legend_labels[color],
                          markerfacecolor=color, markersize=10) for color in label_colors.values() if legend_labels[color] in selected_labels]
    #include legend:
    #plt.legend(handles=handles)
    
    plt.savefig(f"UMAP_{'_'.join(selected_labels)}-check.pdf")
    plt.close()
