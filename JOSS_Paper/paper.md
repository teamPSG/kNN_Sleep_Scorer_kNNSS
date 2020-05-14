---
title: 'Automated Sleep Stage Scoring Using *k*-Nearest Neighbors Classifier'

tags:
  - Supervised clustering
  - Polysomnography
  - Power spectrum
  - MATLAB

authors:
  - name: Tamás Kiss
    orcid: 0000-0001-6360-0714
    affiliation: 1,2

  - name: Stephen Morairty
    orcid: 0000-0002-0781-1645
    affiliation: 3

  - name: Michael Schwartz
    orcid: 0000-0002-5464-638X
    affiliation: 3

  - name: Thomas S.\ Kilduff
    orcid: 0000-0002-6823-0094
    affiliation: 3

  - name: Derek L.\ Buhl
    orcid: 0000-0003-4433-7150
    affiliation: 1,4

  - name: Dmitri Volfson
    orcid: 0000-0002-5167-7834
    affiliation: 1,4

affiliations:
 - name: Global Research and Development, Pfizer Inc, Groton, CT, USA
   index: 1

 - name: Department of Computational Sciences, Wigner Research Centre for Physics, Budapest, Hungary
   index: 2

 - name: Center for Neuroscience, SRI International, Menlo Park, CA, USA
   index: 3

 - name: Takeda Pharmaceuticals, Inc., Cambridge, MA, USA
   index: 4

date: 1 May 2020
bibliography: paper.bib
---

# Polysomnographic Sleep Stage Scoring

Many features of sleep, such as the existence of rapid eye movement
(REM) sleep or non-REM sleep stages, as well as some of the underlying
physiological mechanisms controlling sleep, are conserved across
different mammalian species. Sleep research is important to
understanding the impact of disease on circadian biology and optimal
waking performance, and to advance treatments for sleep disorders,
such as narcolepsy, shift work disorder, non-24 sleep-wake disorder,
and neurodegenerative disease.  Given the evolutionary relatedness of
mammalian species, sleep architecture and changes therein may provide
reliable translational biomarkers for pharmacological engagement in
proof-of-mechanism clinical studies.

Key physiological indicators in sleep include electroencephalography
(EEG) or electrocorticography, electrooculography (EOG), and
electromyography (EMG). Polysomnography (PSG) is the simultaneous
collection of some or all of these measurements and is typically
performed in a specialized sleep laboratory. Determination of the wake
or sleep stage someone is in (i.e., wake, REM sleep, or non-REM sleep,
which is broken down into stages 1, 2, or 3), relies on the judgment
of a trained professional who scores the data based on the
standardized criteria for the recording and staging of human PSG set
forth by @AASMM. Disagreement between individual recordings might
arise due to differences in instrumentation or to the subjective
opinion of the individual scoring the stages. Animal sleep studies
show even greater variability [@ROBERT1999111], as each laboratory
uses methods that best suit their individual needs (e.g.,
electrode/reference positions, muscle choice for EMG implantation, use
of EOG, etc.). While these technical differences make it difficult to
compare studies, the variability in scoring of sleep stages makes it
even more challenging. Although numerous scoring algorithms exist,
most are unreliable, especially following drug treatment. After nearly
half a century of PSG studies, the gold standard of scoring sleep
architecture remains a complete and thorough examination of the PSG
signals, which are scored in 4-, 10-, or 12-second epochs in animal
studies and 30-second epochs in human studies, making it very
difficult to screen through drugs in animal studies and cumbersome to
implement large clinical trials.

# Applications and Advantage

To expedite the tedious process of visually analyzing PSG signals and
to further objectivity in the scoring procedure, a number of sleep
staging algorithms have been developed both for animals
[@STEPHENSON2009263, @BASTIANINI2014277,@vladimir2020] and human
subjects [@PENZEL2000131,@gunn2020,@zhang2020] as reviewed most
recently by @fiorillo2019 and @faust2019. However, computer-based
methods are typically tested on data obtained from healthy subjects or
control animals, and performance is assessed only in a few cases in
subjects with sleep disorders or following drug treatment
[@BOOSTANI201777, @allocca2019]. Furthermore, scoring sleep for
hundreds of animals in a typical preclinical drug discovery effort
often becomes a bottleneck and a potential source of
subjectivity affecting research outcomes.

In this paper, we present an automated approach intended to eliminate
these potential issues. Building on features classically extracted
from EEG and EMG data and machine learning-based classification of
PSG, this approach is capable of staging sleep in multiple species
under control and drug-treated conditions, facilitating the detection
of treatment-induced changes or other manipulations (e.g.,
genetic). Using human interpretable features calculated from EEG and
EMG will be important to understand drug mechanisms, for prediction of
treatment outcomes, and as biomarkers or even translational
biomarkers. However, the initial application of our approach is for
basic and discovery research in which experiments are conducted in
large cohorts of rodents, with the expectation that results can be
translated to higher-order mammals or humans.

Multiple software applications have been developed to address the
problem of automated sleep stage scoring. In their comparative review,
@BOOSTANI201777 found that the best results could be achieved when
entropy of wavelet coefficients along with a random forest classifier
were chosen as feature and classifier, respectively. Another recent
method [@miladinovic2019] used cutting-edge machine learning methods
combining a convolutional neural network-based architecture to produce
domain invariant predictions integrated with a hidden Markov model to
constrain state dynamics based upon known sleep physiology. While our
method also builds on machine learning techniques, it is based on
interpretable features [@wiki:XAI] and uses a simpler algorithm for
classification -- which should make it an ideal choice for the broader
community as well as for sleep experts who might not be too familiar
with complex machine learning approaches. Furthermore, we chose not to
constrain the number of identifiable sleep/wake states or the
probability of transition from one state to another, as we and others
have found that drug interventions and disease processes tend to
change not only the amount of different sleep stages but their
transition probabilities as well. Finally, our method is a supervised
method that requires a training set. While this might seem to be a
disadvantage over non-supervised methods, we have found that drug
treatment or pathological conditions can result in sleep stages not
observed in healthy controls. Thus, the algorithm must be trained to
these new stages.

# Brief Software Description

Our software package, implemented in Matlab, is available for download
on GitHub [@kNNSS]. Automatic sleep staging consists of the classical
consecutive steps of machine learning-based sleep scoring algorithms
\autoref{fig:method}. First, offline stored EEG and EMG data are
loaded into memory to allow for the uniform processing of time-series
data and segmented into consecutive 10-second, non-overlapping epochs
that correspond to manually scored epochs. Second, features are
extracted from the raw signal for all epochs. Features consist of the
power contained in physiologically-relevant frequency bands, as well
as Hjorth parameters for both EEG and EMG data. Third, features
undergo a pre-processing step including the following operations:
unusable epochs that contain too much noise or contain no signal are
removed.  Features are then transformed using the logarithm function
making feature distributions more Gaussian-like, thereby
facilitating subsequent machine classification.  Finally, each feature
is normalized to its median wake value within an animal to enable
usability of the algorithm across laboratories.  Wake periods can be
identified before running the algorithm using the manually-scored
training set or an experiment can be performed such that a given
period is expected to be comprised of an extended period of
wakefulness.  Following feature extraction, a combined filter and
wrapper method-based feature selection step is applied.  This step
ensures that features with the most predictive value are chosen and
also helps to prevent over-fitting. For classification, the
*k*-nearest neighbors classifier is used on data pre-processed
following the procedure described above.

![Summary of training and using the *k*-nearest neighbors algorithm
 for predicting sleep stage labels.\label{fig:method}](fig1.png)

The algorithm was used to predict sleep stages in mice, rats and
non-human primates. Prediction accuracy was found to depend on a
number of parameters of the input data, including consistency of
manual scores and physiological signals, as well as the amount of
artifacts. Furthermore, relative frequency of predicted labels can
influence efficacy, with rare labels being harder to predict. The code
on GitHub [@kNNSS] accompanying this paper contains the abridged
version of a dataset described in detail in @schwartz2018. Three
labels were predicted: wake (W), non-REM sleep (NR), and REM sleep
(R), and prediction efficacy was calculated (\autoref{fig:effic}). The
model was first used to train a single classifier merging training
data from all animals (\autoref{fig:effic}, A), then individual models
were trained, one for each animal (\autoref{fig:effic}, B). The GitHub
repository includes additional information on prediction accuracy,
including detailed values of true and false positive rates, as well as
a method to deal with imbalanced data.

![Estimation of prediction accuracy. For each state (wake -- W,
 non-REM -- NR, REM -- R) and animal (labeled asterisks) true and
 false positive rates are calculated. Red crosses denote mean and SEM.
 In A, training data was merged and one single classifier was trained
 to predict sleep stages of all animals. In B, an individual
 classifier was trained for each animal
 separately.\label{fig:effic}](fig2.png)

# Acknowledgments

TK, DV, and DLB were full time employees and shareholders of Pfizer
Inc. during development of this software package. This work was supported
by Pfizer Inc. and SRI International.

# References