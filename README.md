# ODT with Lippmann-Schwinger model

## Description

This repository provides codes for 2D and 3D optical diffraction tomography (ODT) with iterative Lippmann-Schwinger model. More details can be found in the following paper:

[1] <a href="https://www.osapublishing.org/oe/abstract.cfm?uri=oe-25-18-21786" target="_blank">Efficient inversion of multiple-scattering model for optical diffraction tomography.</a>, <br />
Optics Express, 25, 21786–21800 (2017), <br />
E. Soubies, T-a. Pham, and M. Unser.

[2] <a href="https://ieeexplore.ieee.org/document/8970570" target="_blank">Three-Dimensional Optical Diffraction Tomography with Lippmann-Schwinger Model.</a>, <br />
IEEE Transactions on Computational Imaging (2020), <br />
T-a. Pham, E. Soubies, A. Ayoub, J. Lim, D. Psaltis, and M. Unser.

[3] <a href="https://www.sciencedirect.com/science/article/pii/S003040181930776X" target="_blank">A method for assessing the fidelity of optical diffraction tomography reconstruction methods using structured illumination.</a>, <br />
Optics Communications, 454, 124486 (20209, <br />
A. B. Ayoub, T-a. Pham, J. Lim, M. Unser, and D. Psaltis.

## Requirements

The code requires the GlobalBioIm library v1.1.2 (or more recent releases) <br />
https://biomedical-imaging-group.github.io/GlobalBioIm/

## Repository content

The repository is organized as follows.
* Folder 2D: two-dimensional ODT [1]
  * In folder code, the script **main_fig7.m** reproduces Figure 7 of [1].
  * In folder code, the script **ExampleFwdModel.m** provides an example on how to use the iterative Lippmann-Schwinger forward model, 
  * The folder mat_files contains data. 
* Folder 3D: three-dimensional ODT [2]
  * The script **DiscretizationGreen.m** reproduces Figure 3 of [2],
  * The script **ExampleFwdModel.m** provides an example on how to use the iterative Lippmann-Schwinger forward model, 
  * The script **Compute_uin.m** reproduces Figure 5 of [2],
  * The script **ODT_Reconstruction.m** contains the whole reconstruction pipeline.
  * The script **Assessment_method.m** implements the method of assessment of [3]
  **Note** The data for ODT_Reconstruction.m can be downloaded here: https://www.irit.fr/~Emmanuel.Soubies/Data/Data_LS3D.zip
