# Machine-learning-driven modelling of amorphous and polycrystalline BaZrS₃

Chalcogenide perovskites, and in particular the sulfide-based BaZrS₃ composition, has gained significant attention for its potential applications in thin-film photovoltaics and catalysis. 
The present study shows how ML-driven modelling can be leveraged to study both the material's amorphous precursor phase and polycrystalline structures with complex grain boundaries.
Using a tailored machine-learned interatomic potential (MLIP) for BaZrS₃ based on the efficient Atomic Cluster Expansion (ACE) framework, we investigate the atomic structure of the amorphous phase, evaluate grain-boundary formation energies, and generate realistic polycrystalline models for comparison with experimental data. Beyond the explicit BaZrS₃ composition, the potential fitting and analysis methods proposed in this work can be more widely applied to model polycrystalline functional materials, and make a step towards advancing large-scale simulations of emerging optoelectronic and photovoltaic materials.
This repository provides access to the trained ACE models, the dataset used for training, and the analysis code supporting the results presented in the paper.

## Dataset construction and fitting protocol 

The training dataset is constructed using an established protocol based on de novo structures, using an initial Random Structure Search (RSS) training dataset that is then iteratively expanded using a potential fitted on the growing dataset (in this case, GAP-RSS). "Domain-specific" structures, such as crystalline–amorphous interfaces or high-temperature structures are additionally included  to ensure that the training dataset spans the relevant chemical space. To visualise the composition of the dataset interactively, please see the Jupyter Notebook in the [dataset directory](https://github.com/BiancaPasca/polycrystalline-BaZrS3/blob/main/dataset/UMAP-visualise.ipynb).


## Bespoke potential for polycrystalline BaZrS₃

The MLIP enables computational exploration of the different phases of BaZrS₃ in a sequence that aligns with experimental synthesis processes.
It can first be applied to simulate the amorphous phase, representing precursor states formed during deposition.
The potential also facilitates a quantitative analysis of grain boundaries, ensuring their accurate description for modeling polycrystalline structures.
With both the amorphous and polycrystalline phases accessible, the model further allows for running simulations of structures with varying grain sizes, offering a direct link to experimental scattering data, as demonstrated in the publication. 

Publication currently in preparation.
