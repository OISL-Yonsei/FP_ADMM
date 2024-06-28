# Fourier Ptychogrphic Phase Retrieval with Convex Optimization

Repo for the project - **Anisotropic Regularization for Sparsely Sampled and Noise-Robust Fourier Ptychography**
#### [Project Page](https://kyungchullee.com/projects/6_project/) | [Paper](https://opg.optica.org/oe/fulltext.cfm?uri=oe-32-14-25343&id=552914)


[Kyung Chul Lee](https://kyungchullee.com/)\ $^{1}$, Hyesuk Chae\ $^{1}$, Shiqi Xu\ $^{2}$, Kyungwon Lee\ $^{1}$, Roarke Horstmeyer\ $^{2}$, Byung-Woo Hong\ $^{3}$, [Seung Ah Lee](https://biomedia.yonsei.ac.kr/) $^{1}$. <br>
$^1$ Yonsei University, $^2$ Duke University, $^3$ Chung-Ang University.

:pushpin: Related paper accepted to [SIGGRAPH ASIA 2023](https://asia.siggraph.org/2023/submissions/technical-papers/).




--------------
## What We Contribute?

- derivation of the algorithm for TV-regularized FP reconstruction algorithm via ADMM.
- validation of Tikhonov regularization on the pupil function that ensures a smooth pupil phase profile.
- demonstration of the performance of our method with extremely low measurement SNR and Fourier-space subsampling for reduced number of measurements. 
- application of our algorithm using real experimental data, which shows the effectiveness of TV-regularizer for both super-resolution and phase imaging of biological specimens under challenging imaging conditions.



## How to Use this Repo?
See details in [code.md]().



## Citation

If you find our code or any of our materials useful, please cite our paper:
```bibtex
@article{Lee:24,
author = {Kyung Chul Lee and Hyesuk Chae and Shiqi Xu and Kyungwon Lee and Roarke Horstmeyer and Seung Ah Lee and Byung-Woo Hong},
journal = {Opt. Express},
keywords = {Computational imaging; Imaging systems; Imaging techniques; Phase imaging; Spatial resolution; Superresolution},
number = {14},
pages = {25343--25361},
publisher = {Optica Publishing Group},
title = {Anisotropic regularization for sparsely sampled and noise-robust Fourier ptychography},
volume = {32},
month = {Jul},
year = {2024},
url = {https://opg.optica.org/oe/abstract.cfm?URI=oe-32-14-25343},
doi = {10.1364/OE.529023},
}
```
