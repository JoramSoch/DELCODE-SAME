# DELCODE_SAME

**FADE and SAME scores for multi-center DELCODE study**

This code belongs to the paper "Single-value fMRI scores reflect both stage and risk across the Alzheimer's disease continuum" by Soch et al. (2021), currently under review at *Brain* and publicly available from *medRxiv*. For instructions how to process these data, see below.

- Preprint: TBA
- data and code: https://github.com/JoramSoch/DELCODE_SAME
- original study: https://github.com/JoramSoch/FADE_SAME


### Requirements

This code was developed and run using the following software:
- Windows 10 Professional 64-bit
- [MATLAB R2021a](https://de.mathworks.com/help/matlab/release-notes.html) (Version 9.10)
- [MATLAB Statistics and ML Toolbox](https://de.mathworks.com/products/statistics.html) (Version 12.1)
- [ME_GLM](https://github.com/JoramSoch/MACS/blob/master/ME_GLM.m) function from the [MACS](https://github.com/JoramSoch/MACS) toolbox


### Overview

This repository contains the following files:

* `data/fMRI_scores.xls`: 4 single-value fMRI scores computed from 468 DELCODE subjects (used by virtually all scripts in this directory);
* `data/fMRI_scores_yFADE_ref.xls`: fMRI scores of the same subjects when calculated using a different reference sample (only used by `Figure_FADE_SAME_stability.m`);
* `data/fMRI_scores_AiA_subj.xls`: 4 fMRI scores for 259 subjects from the [original study](https://onlinelibrary.wiley.com/doi/10.1002/hbm.25559) (obtained using [FADE_SAME](https://github.com/JoramSoch/FADE_SAME), only used by `Figure_FADE_SAME.m` and `Figure_FADE_SAME_by_age.m`);
* `subjects/bhvr_data.xls`: behavioral response measures within the fMRI task for all 468 subjects (only used by `Figure_FADE_SAME_corrs.m`);
* `subjects/subj_covs.mat`: a number of covariates obtained for those subjects (used by virtually each script in this directory);
* `tools/*.m`: a number of in-house MATLAB helper functions assisting statistical analysis of the data
* `Figure_*.m`: scripts reproducing Figures from the main manuscript or supplementary material;
* `Table_*.m`: scripts reproducing Tables from the main manuscript or supplementary material;
* `Table_analyses_methods.pdf`: an overview of all Figures and Tables in manuscript and supplement (and their correspondence to the [original study](https://onlinelibrary.wiley.com/doi/10.1002/hbm.25559));
* `Figure_FADE_SAME.png`: a screenshot of the main Figure from the Paper (used as a [Graphical Abstract](#graphical-abstract) in this repository).


### Instructions

For re-running analyses reported in the paper, you need to perform the following steps:
1. Clone this repository to some folder on your local computer.
2. Open MATLAB and change your current directory to this folder.
3. Run scripts `Figure_*.m` or `Table_*.m` to reproduce Figures and Tables from the paper.

* For creating Figure S3 (`Figure_FADE_SAME_stability.m`), you need to either have the [MACS toolbox](https://github.com/JoramSoch/MACS) on your MATLAB path or download the [ME_GLM function](https://github.com/JoramSoch/MACS/blob/master/ME_GLM.m) into the repository folder.


### Graphical Abstract

<img src="https://github.com/JoramSoch/DELCODE_SAME/raw/master/Figure_FADE_SAME.png" alt="Graphical Abstract" width=1000>