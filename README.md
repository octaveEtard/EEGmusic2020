# JoNmusic2020

Analysis code for the Journal of Neuroscience submission setting out to investigate neural responses to the temporal fine structure of continuous musical pieces in a bipolar EEG dataset. The corresponding dataset is available on [Zenodo](https://zenodo.org/).

Requires Matlab R2019b or newer (tested on R2019b & R2020a/b).

## Installation
This code is based on the [LM package](https://github.com/octaveEtard/LMpackage), and requires it to be available in your Matlab path.

Add the `functions` folder to your path. The code is structured as a main [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) contained in `functions/+JoNmusic2020` to reduce the risk of shadowing any of your own functions. This package implements the analysis of the EEG data. A second auxilliary package (`functions/+pltools`) provides functions to plot the results.

## Description
This repository is structured as follow:

`functions/+JoNmusic2020`: contains analysis functions. Essentially wrapping functions to load the relevant EEG and stimulus data and pass them to the `LMpackage` functions, as well as cross-validation procedures.

`functions/+pltools`: contains plotting functions

`forwardModels`: contains the main scripts running the following analyses:

  * `music_SI_forward.m`: linear forward model relating the stimulus waveforms to the EEG data in each of the SI conditions
  * `music_SI_forward_pooled.m`: linear forward model relating the stimulus waveforms to the EEG data pooled over the two SI conditions
  * `music_CI_forward.m `: linear forward model relating the stimulus waveforms to the EEG data pooled over the two CI conditions
