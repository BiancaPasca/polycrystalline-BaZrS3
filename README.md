# Machine-learning-driven modelling of amorphous and polycrystalline BaZrS3

The chalcogenide perovskite material BaZrS3 is of growing interest for emerging thin-film photovoltaics, as well as in the field of catalysis. 
In this work we show how ML-driven modelling can be used to describe the material’s amorphous precursor, together with polycrystalline structures containing complex grain boundaries. 
Using a bespoke, multi-purpose machine-learned interatomic potential (MLIP) model for BaZrS3,we study the atomic structure of the amorphous phase, quantify grain-boundary
formation energies, and create realistic-scale polycrystalline structural models which can be compared to experimental data. Beyond BaZrS3, our work marks a step 
towards realistic device-scale simulations of emerging optoelectronic and photovoltaic materials. 
In this repository, the trained ACE models, the constructed dataset and the code used to analyse the data presented in the paper are made publicly available. 

## Dataset construction and fitting protocol 

<img width="1068" alt="overview_plot" src="https://github.com/user-attachments/assets/444d4e67-3d67-4f12-8b6f-d691937cc1c5" />


In this work, we use an established protocol for training the potential starting from de novo structures, using an initial Random Structure Search (RSS)
training dataset that is then iteratively expanded using a potential fitted on the growing dataset (in this case, GAP-RSS). "Domain-specific" structures, such as
high-temperature snapshots obtained from ML-driven MD, as well as crystalline–amorphous interfaces were further added to ensure that the training dataset spans the
relevant chemical space.
To [visualise the composition of the dataset interactively]([url](https://jupyter.org/try-jupyter/notebooks/?path=UMAP-visualise.ipynb)), please see the Jupyter Notebook in the [dataset directory]([url](https://github.com/BiancaPasca/polycrystalline-BaZrS3/blob/main/dataset/UMAP-visualise.ipynb)). 


## Bespoke potential for polycrystalline BaZrS3

The potential can be used to computationally study BaZrS3 in the same sequence as would be relevant in experiment. 
First, the MLIP can be used simulate the amorphous phase, corresponding to precursor phases that have been deposited in experimental synthesis routes. 
The potential can also be used in the quantitative study of grain boundaries, which need to be accurately described so that the model can be applied to polycrystalline samples. 
Given that both of these are available, the model proves to additionally be useful for simulating structures with different grain sizes, providing a direct connection to experimental scattering data,
as shown in the publication.
