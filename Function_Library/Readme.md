# High-level function description for *k*NNSS

This readme describes what each function of the *k*NNSS package is used for. A detailed description of mandatory and optional inputs for functions, as well as their outputs can be found in the function headers.

## Top-level functions

generate_statespace.m
:  This is the first function to use in the processing chain. It loads raw data, cuts it into epochs and calculates feature describing each epoch. It also loads manual scores associated with the raw data and labels each epoch using manual scores (optional). Feature data (since its calculation can be time consuming) can be saved in Matlab files. This function is written with modularity in mind: if further features are to be included their calculation can be implemented in a separate function and values passed to generate_statespace that merges them and saves a unified feature file. An example use is shown in the FirstStep.m script in SoftwareVerification folder.

train_one_model.m
: This function uses features calculated by generate_statespace.m and trains a *k*-Nearest Neighbors (*k*NN) classifier. This function gathers feature information from multiple experiments (could be from the same or from different subjects) and uses the merged data to train a single classifier for all experiments. This way a general classifier is created that can be used across multiple experiments. Its use is presented in the SecondStep.m script in SoftwareVerification folder.

train_many_models.m
: Similarly to train_one_model.m, this function uses pre-calculated feature values to train a *k*NN classifier. It differs from train_one_model.m in that it trains a classifier for each experiment it gets. Individual classifiers are more accurate for any given subject, however, they do not generalize as well as models trained on multiple subjects or conditions. Its use is presented in the SecondStep.m script in SoftwareVerification folder.

evaluate_model_goodness.m
: This function calculates confusion matrices from manually and automatically determined label sets and displays a set of measures to describe how well the classifiers perform.

## Helper functions
### Functions from the [Chronux toolbox](http://chronux.org)  used for spectrum  calculations
change_row_to_column.m
: Helper routine to transform 1d arrays into column vectors that are needed by other routines in Chronux

dpsschk.m
: Helper function to calculate tapers

getfgrid.m
:  Helper function that gets the frequency grid associated with a given FFT based computation

getparams.m
: Helper function to convert structure params to variables used by the various routines in Chronux

mtfftc.m
: Multi-taper fourier transform - continuous data

mtspecgramc.m
: Multi-taper time-frequency spectrum - continuous process

mtspectrumc.m
: Multi-taper spectrum - continuous process

### Other helper functions

edfread.m
: Read European Data Format file into MATLAB. This function was written by Brett Shoelson, PhD (brett.shoelson@mathworks.com).

canonize_fieldname.m
: This function attempts to create properly formatted field names.

confusion2PerformanceMetrics.m
: This function converts the confusion matrix into different preformance metrics including SEN, FPR, SPC, ACC, PPV, NPR. For details see the [Wikipedia page for confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix).

estimate_feature_goodness.m
: This function uses [filter methods](https://en.wikipedia.org/wiki/Feature_selection#Filter_method) to select features used for classifier training that are most dissimilar.

hjorth.m
: Calculates the [Hjorth parameters](https://en.wikipedia.org/wiki/Hjorth_parameters) activity, mobility, and complexity. This function is part of the [BIOSIG-toolbox](http://biosig.sf.net/) and was written by Alois Schloegl (a.schloegl@ieee.org).

select_features.m
: The pair of estimate_feature_goodness.m, this function uses the [wrapper method](https://en.wikipedia.org/wiki/Feature_selection#Wrapper_method) (sequential feature selection here) to find best features to be used for training the classifier.

flag_outliers.m
: This function finds outliers in a feature table using median and mad.

metric_band_power.m
:  This function calculates power in given frequency bands using multitaper power estimation.

preprocess_features.m
: Feature pre-processing is performed by this function. It is called before fitting a classifier or before predicting new data. Pre-processing steps include removal of data with artifacts, missing data segments, removal of manually selected epochs, normalization, and log transformation.

train_classifier.m
: This function trains the KNN classifier.