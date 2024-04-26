# Glare_Suppression_AlterProjTM

This repository stores the simulation codes for our recently published work "[Alternating projection-based phase optimization for arbitrary glare suppression through multimode fiber](https://www.sciencedirect.com/science/article/pii/S0143816622004195?via%3Dihub)"


**Abstract:** Wavefront shaping for glare suppression helps to reduce the speckle background through scattering media and enables imaging, sensing, and some speckle-related advanced applications through customizing the speckle light field. For glare suppression in a target region, current methods are either slow or not sufficiently generic. Here, an alternating projection method that fully exploits the transmission matrix of the scattering medium is proposed for fast and arbitrary glare suppression. Parallelly, multiple phase masks corresponding to different target regions can be computationally optimized without iterative hardware feedback. Numerical simulation shows that under phase-only modulation, a suppression factor of ~10^(-3) can be achieved for a large suppression region. With the use of a graphics processing unit and a digital micromirror device, experiments demonstrate fast and effective glare suppression for target regions of various shapes and sizes at the distal end of a multimode fiber. This technique could be promising for applications like speckle optical tweezer or endoscopy with speckle illumination inside deep tissues or other complex environments.   

## Schematic of the TM-based alternating projection method
<img src="https://github.com/Ford666/Glare-Suppression-AlterProjTM/blob/main/images/principle_of_AlterProjTM.png" width="800px">

## Simulated results of glare suppression with ER and HIO constraints

<img src="https://github.com/Ford666/Glare_Suppression_AlterProjTM/blob/main/images/glarSupn_ER_HIO.png" width="800px">

## Citation
If you think the codes are helpful for your research, please consider to cite our work:
```
@article{CHENG2023107368,
title = {Alternating projection-based phase optimization for arbitrary glare suppression through multimode fiber},
journal = {Optics and Lasers in Engineering},
volume = {161},
pages = {107368},
year = {2023},
issn = {0143-8166},
doi = {https://doi.org/10.1016/j.optlaseng.2022.107368},
url = {https://www.sciencedirect.com/science/article/pii/S0143816622004195},
author = {Shengfu Cheng and Tianting Zhong and Chi Man Woo and Puxiang Lai},
}
```
