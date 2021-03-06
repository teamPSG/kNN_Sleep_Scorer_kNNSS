# *k*-Nearest Neighbors Classifier-Based Sleep Staging (*k*NNSS)

## Welcome to *k*NNSS software

*k*NNSS is a Matlab package for automated sleep stage scoring using the *k*-nearest neighbors algorithm. Compared to other automated sleep scoring packages its main advantage is simplicity and the use of physiologically relevant, human-interpretable features. The code is documented in the following JOSS publication: [![DOI](https://joss.theoj.org/papers/10.21105/joss.02377/status.svg)](https://doi.org/10.21105/joss.02377)

## LICENSE

    /* This Source Code Form is subject to the terms of the Mozilla Public
     * License, v. 2.0. If a copy of the MPL was not distributed with this
     * file, You can obtain one at http://mozilla.org/MPL/2.0/.

## Using *k*NNSS
### Compatibility and dependencies
*k*NNSS was tested in Matlab R2014b and R2018a. The software uses the Signal Processing Toolbox and the Statistics and Machine Learning Toolbox from Matlab. Also, some functions from the [Chronux Toolbox](http://chronux.org)  are used to calculate spectral properties of data (see below for a list of files). European Data Format (EDF; https://www.edfplus.info/) files are read using Brett Shoelson's (brett.shoelson@mathworks.com) edfread() function. Hjorth parameters are calculated by hjorth() a function written by Alois Schloegl (a.schloegl@ieee.org) as part of the [BIOSIG-toolbox](http://biosig.sf.net). This GitHub project is self-containing, no further downloads are required.

#### Chronux Toolbox file list
change_row_to_column.m, dpsschk.m, getfgrid.m, getparams.m, mtfftc.m, mtspecgramc.m, mtspectrumc.m

### Folder content

- Example_Data: this directory contains three sub-directories. IntermRes is initially empty, intermediate files generated during the execution of the test scripts will be stored here. RawData contains polysomnographycal recordings from Thomas Kilduff's laboratory at SRI International in EDF format. Sub-folder EDF contain raw recordings, sub-folder FFT stores associated manual sleep scores used to train classifiers and test prediction efficacy.

- Function_Library: this directory stores all the necessary Matlab functions required by *k*NNSS. Functions inputs and outputs are described in detail in the function header. Use Matlab's help to read about details of each function. A list of high-level purpose of functions is found in the Readme in this folder.

- Software_Verification: a set of scripts to test the code on your system are in this folder. For a detailed description see the Readme in this folder.

- JOSS_Paper: contains the manuscript submitted to the Journal of Open Source Software. The journal paper will provide background to the software, a summary of our work, applications and references.

### Installation
All functions are stored in the Function_Library folder. Simply adding this folder to Matlab's function path will take care of installation. Alternatively, check out contents of the Software_Verification folder.

## Top-level inputs and outputs
*k*NNSS reads raw electrophysiological data from file. Currently the European Data Format and Matlab files are supported but see generate_statespace.m in Function_Library to see how additional formats can be added. Manual scores are read together with the data. These are typically text files that come in many format.

Once the algorithm finishes predicting labels the output is a Matlab structure type variable, containing a cell array of strings in its fields. Strings correspond to predicted labels.

## Getting help
For more information please do not hesitate to contact Tamás Kiss (kiss.t (at) wigner.hu).

Thank you!
