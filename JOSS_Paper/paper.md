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
    orcid: 0000-0000-0000-0000
    affiliation: 3

  - name: Michael Schwartz
    orcid: 0000-0002-5464-638X
    affiliation: 3

  - name: Thomas Kilduff
    orcid: 0000-0002-6823-0094
    affiliation: 3

  - name: Derek L. Buhl
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
physiological mechanisms controlling sleep are conserved across
different mammalian species. Sleep research is important to
understanding the impact of disease on circadian biology and to
advance treatments for sleep disorders, such as narcolepsy, shift work
disorder, non-24 sleep-wake disorder, and neurodegenerative disease.
Given the bridge between species, sleep architecture and changes
therein may provide a reliable translational biomarker for
pharmacological engagement proof-of-mechanism clinical studies.

Key physiological indicators in sleep include electroencephalography
(EEG) or electrocorticography, electrooculography (EOG), and
electromyography (EMG). Polysomnography (PSG) is the simultaneous
collection of some or all of these measurements and is typically
performed in a specialized sleep laboratory. Determination of the wake
or sleep stage someone is in (i.e., wake, REM sleep, or non-REM sleep,
which is broken down into stage 1, 2, or 3), relies solely on a
trained technician who scores the data based on the standardized
criteria for the recording and staging of human PSG set forth by
@AASMM. Disagreement between individual recordings might arise due to
differences in instrumentation or the subjective opinion of the
individual scoring the stages. Animal sleep studies show an even
greater variability [@ROBERT1999111], as each laboratory uses methods
that best suit their individual needs (e.g., electrode/reference
positions, muscle choice for EMG implantation, use of EOG,
etc.). While these differences make it difficult to compare studies,
the scoring of sleep stages makes it even more challenging. Although
numerous scoring algorithms exist, most are unreliable, especially
following drug treatment. After nearly half a century of PSG studies,
the gold standard of scoring sleep architecture remains a complete and
thorough examination of the PSG signals, typically in 30-second
epochs, which is next to impossible to perform to screen through drugs
in animal studies and cumbersome to implement for large clinical
trials.

# Applications and Advantage

To shorten the tedious process of visually analyzing PSG signals and
to further objectify the scoring procedure, a number of sleep staging
algorithms have been developed both for human subjects
[@PENZEL2000131] and animals [@STEPHENSON2009263,
@BASTIANINI2014277]. However, computer-based methods are typically
tested on data obtained from healthy subjects or control animals, and
performance is assessed only in a few cases in subjects with sleep
disorders [@BOOSTANI201777]. Furthermore, scoring sleep for hundreds
of animals in a typical preclinical drug discovery effort often
becomes the experimental bottleneck and a potential source of
subjectiveness affecting research outcomes.

In this paper we present an automated approach that eliminates these
potential issues. Using the same method, building on features
classically extracted from EEG and EMG data and machine learning-based
classification of PSG, this approach is capable of staging sleep in
multiple species under control and drug-treated conditions,
facilitating the detection of treatment-induced changes or
manipulations (e.g.\ genetic). It is worth to be noted that using
human interpretable features calculated from EEG and EMG is important
to understand drug mechanisms, prediction of treatment outcomes, and
can potentially be used as biomarkers or even translational
biomarkers. The primary use of our approach is in basic and discovery
research, in which most experiments are conducted in large cohorts of
rodents and interpretation of results are being translated to
higher-order mammals or humans.

Multiple software applications have been developed to solve the
problem of automated sleep stage scoring. In their comparative review
@BOOSTANI201777 found that the best results could be achieved when
entropy of wavelet coefficients along with random forest classifier
were chosen as feature and classifier, respectively. Another recent
method [@miladinovic2019] uses cutting-edge machine learning methods
combining a convolutional neural network-based architecture to produce
domain invariant predictions integrated with a hidden Markov model to
constrain state dynamics based upon known sleep physiology. While our
method also builds on machine learning techniques, it is based on
interpretable features [@wiki:XAI] and uses simpler algorithm for
classification making it the ideal choice for a broader community or
sleep experts not too familiar with complex machine learning
approaches. Furthermore, we chose not to constrain state label
sequences as we found that drug interventions as well as disease
processes tend to change not only the amount of different sleep stages
but their transition probability also. Finally, our method is a
supervised method that requires a training set. While this might seem
to be a disadvantage over non-supervised methods, we also found that
drug treatment or pathological conditions can result in sleep stages
not observed in healthy controls, thus these new stages need to be
trained to the algorithm.

# Brief Software Description

Our software package, implemented in Matlab, is available for download
on GitHub [@kNNSS]. Automatic sleep staging consists of the classical
consecutive steps of machine learning-based sleep scoring algorithms
\autoref{fig:method}. First, offline stored EEG and EMG data are
loaded into memory to allow for the uniform processing of time-series
data and segmented into 10 second long, non-overlapping, consecutive
epochs that fit manually scored epochs. Second, features are extracted
from the raw signal for all epochs. Features consist of the power
contained in given, physiologically relevant frequency bands, as well
as Hjorth parameters for both EEG and EMG data. Third, features
undergo a pre-processing step including the following operations:
unusable epochs that contain too much noise or contain no signal are
removed.  Then features are transformed using the logarithm function
making feature distributions be more Gaussian-like helping subsequent
machine classification.  Finally, each feature is normalized to its
median wake value within an animal to enable inter-center usability of
the algorithm.  Wake periods can initially be identified from the
manually scored training set or the experiment can be performed such
that a given period is prepared to contain wake periods only.
Following feature extraction a combined filter and wrapper
method-based feature selection step is applied.  This step ensures
that features with the most predictive value are chosen and also helps
to prevent over-fitting. For classification, the *k*-nearest neighbors
classifier is used on data pre-processed the way described above.

![Summary of training and using the *k*-nearest neighbors algorithm
 for predicting sleep stage labels.\label{fig:method}](fig1.png)

# Acknowledgments

TK, DV, and DLB were full time employees and shareholders of Pfizer
Inc. during development of this software package. This work was funded
by Pfizer Inc.

# References