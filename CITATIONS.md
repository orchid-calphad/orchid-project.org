# Software and data citations

These credits describe the ORCHID project and its examples. Cite the papers and software used in your analysis, together with the ORCHID release and original source of each observation. The list does not imply that every dataset was generated with every package. Preserve the versions and commit hashes recorded by each example; a repository link is not a reproducibility pin. Materials Project and OQMD energies retain separate references.

## pycalphad

Thermodynamic calculations and phase equilibria.

- R. Otis and Z.-K. Liu (2017), pycalphad: CALPHAD-based Computational Thermodynamics in Python. [Paper](https://doi.org/10.5334/jors.140).
- [Repository](https://github.com/pycalphad/pycalphad).

## ESPEI

Thermodynamic parameter fitting and uncertainty quantification.

- B. Bocklund et al. (2019), ESPEI for efficient thermodynamic database development, modification, and uncertainty quantification: application to Cu–Mg. [Paper](https://doi.org/10.1557/mrc.2019.59).
- [Repository](https://github.com/PhasesResearchLab/ESPEI).

## PhaseForge

MLIP-informed CALPHAD construction in the examples.

- S. Zhu, D. Sarıtürk, and R. Arróyave (2025), Machine learning potentials for alloys: a detailed workflow to predict phase diagrams and benchmark accuracy. [Paper](https://doi.org/10.1038/s41524-025-01814-z).
- [Repository](https://github.com/dogusariturk/PhaseForge).
- [Archived software](https://doi.org/10.5281/zenodo.15730912). Cite the version actually used.

## PhaseForge+

Construction and experimental tuning of CALPHAD models in the examples; repository name PhaseForgePlus.

- C. Kunselman, S. Zhu, D. Sarıtürk, and R. Arróyave (2025), Construction and Tuning of CALPHAD Models Using Machine-Learned Interatomic Potentials and Experimental Data: A Case Study of the Pt–W System. [Paper](https://doi.org/10.1007/s11669-025-01222-2).
- [Repository](https://github.com/dogusariturk/PhaseForgePlus).
- [Archived software](https://doi.org/10.5281/zenodo.17519087). Cite the version actually used.

## Materials Project

DFT structures and energies; repository link is the API client, not the database itself.

- A. Jain et al. (2013), Commentary: The Materials Project: A materials genome approach to accelerating materials innovation. [Paper](https://doi.org/10.1063/1.4812323).
- [Repository](https://github.com/materialsproject/api).
- [Data resource](https://materialsproject.org).

## OQMD

DFT structures and energies; qmpy is the supporting software.

- J. E. Saal et al. (2013), Materials Design and Discovery with High-Throughput Density Functional Theory: The Open Quantum Materials Database (OQMD). [Paper](https://doi.org/10.1007/s11837-013-0755-4).
- S. Kirklin et al. (2015), The Open Quantum Materials Database (OQMD): assessing the accuracy of DFT formation energies. [Paper](https://doi.org/10.1038/npjcompumats.2015.10).
- [Repository](https://github.com/wolverton-research-group/qmpy).
- [Data resource](https://oqmd.org).

The Cu–Zr gradient workflows additionally credit Courtney Kunselman’s ESPEI and littlemcmc branches, as pinned in the example requirements and provenance manifests: [analytic gradients](https://doi.org/10.1016/j.actamat.2025.121453), [gradient-informed HMC](https://doi.org/10.1016/j.actamat.2026.122412), [ESPEI fork](https://github.com/cjkunselman18/ESPEI), [littlemcmc fork](https://github.com/cjkunselman18/littlemcmc).
