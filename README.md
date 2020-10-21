# EEGmusic2020
Octave Etard (octave.etard11@imperial.ac.uk)

Analysis code to investigate neural responses to the temporal fine structure of continuous musical pieces in a bipolar EEG dataset. The corresponding dataset is available on [Zenodo](https://zenodo.org/).

Requires Matlab R2019b or newer (tested on R2019b & R2020a/b).

## Installation
This code is based on the [LM package](https://github.com/octaveEtard/LMpackage), and requires it to be available in your Matlab path.

Add the `functions` folder to your path. The code is structured as a main [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) contained in `functions/+JoNmusic2020` implementing the analysis of the EEG data. A second auxilliary package (`functions/+pltools`) provides functions to plot the results. The package organisation reduces the risk of shadowing any of the user's own functions.

By default, the following relative location of the analysis code and data is expected:

`someFolder`\
` ├─── JoNmusic2020.git`  (arbitrary name: folder containing this repository)\
` │     ├─── functions`\
` │     ├─── forwardModels`\
` │     ├─── ...`\
` └─── data`         (architecture as implemented in the [Zenodo dataset](https://zenodo.org/))\
`       ├─── EEG`\
`       ├─── stimuli`\
`       ├─── ...`\

This behaviour is implemented in `functions/+JoNmusic2020/getPath.m`, and a different data folder location can be simply specified by editing this function. The expected organisation of the data folder (EEG and feature file locations, etc.) is further specified by the functions `makePathEEGFolder`, `makePathFeatureFiles` and `makePathSaveResults`. The file naming convention is implemented by `makeNameEEGDataFile` and `makeNameEEGDataFile`.

## Description
This repository is structured as follow:

* `functions/+JoNmusic2020`: contains analysis functions. Essentially wrapping functions to load the relevant EEG and stimulus data and pass them to the `LMpackage` functions, as well as cross-validation procedures

* `functions/+pltools`: contains plotting functions

* `forwardModels`: contains the main scripts running the following analyses:

  * `music_SI_forward.m`: linear forward model deriving the neural response to the stimulus waveforms from the EEG data in each of the two SI conditions
  * `music_SI_forward_pooled.m`: linear forward model deriving the neural response to the stimulus waveforms from the EEG data pooled over the two SI conditions
  * `music_CI_forward.m `:linear forward model deriving the attended and ignored neural responses to the stimulus waveforms from the EEG data pooled over the two CI conditions
  
* `backwardModels`: contains the main scripts running the following analyses:

  * `music_SI_backward.m`: linear backward model relating the stimulus waveforms to the EEG data in each of the two SI conditions
  * `music_CI_backward.m `: linear backward model relating the attended and ignored stimulus waveforms to the EEG data in each of the two CI conditions
  
* `behav/compute_perf.m`: compute behavioural performance of the subjects

* `figures`: script plotting the results of the analyses detailed above
