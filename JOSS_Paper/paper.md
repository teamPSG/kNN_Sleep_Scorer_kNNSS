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

  - name: Clinton Bourbonais
    orcid: 0000-0002-7469-8902
    affiliation: 1

  - name: Dmitri Volfson
    orcid: 0000-0002-5167-7834
    affiliation: 1,4

  - name: Derek Buhl
    orcid: 0000-0003-4433-7150
    affiliation: 1,4

  - name: Thomas Kilduff
    orcid: 0000-0002-6823-0094
    affiliation: 3

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

# Polysomnographic Sleep Stage Scoring -- Summary

Sleep research comprises many different areas, including narcolepsy, circadian
rhythms, effects of sleep deprivation and shift work, and aging.  Many features
of sleep, such as the existence of rapid eye movement (REM) sleep or non-REM
sleep stages, as well as some of the underlying physiological mechanisms
controlling sleep are conserved across different mammalian species. Given this
bridge between species, sleep architecture may provide a reliable translational
biomarker for pharmacological engagement proof-of-mechanism clinical studies.

Key physiological indicators in sleep include electroencephalography (EEG) or
electrocorticography (ECoG), electrooculography (EOG), and electromyography
(EMG). Polysomnography (PSG) is the simultaneous collection of some or all of
these measurements and is typically performed in a specialized sleep laboratory
in humans. While standardized criteria for the recording and staging of human
PSG exists [@AASMM], differences in individual recordings might arise due to
different instrumentation or setup. Animal sleep studies show an even greater
variability [@ROBERT1999111], as each laboratory uses methods that best suit
their individual needs (e.g., electrode/reference positions, muscle choice for
EMG implantation, use of EOG, etc.). While these differences make it difficult
to compare studies, the scoring of sleep stages makes it even more
challenging. Although numerous scoring algorithms exist, most are unreliable,
especially following drug treatment. After nearly half a century of PSG studies,
the gold standard of scoring sleep architecture remains a complete and thorough
examination of the PSG signals, typically in 30-second epochs, which is next to
impossible to perform for even a small clinical trial.

To shorten the tedious process of visually analyzing PSG signals and to further
objectify the scoring procedure, a number of sleep staging algorithms have been
developed [@PENZEL2000131], however, computer-based methods are typically tested
on data obtained from healthy subjects and performance is assessed only in a few
cases in subjects with different sleep disorders [@BOOSTANI201777]. In this
paper we present an approach that employs a single method building on features
classically extracted from EEG and EMG data and machine learning-based
classification of PSG. This method is capable of staging sleep in multiple
species under control and drug-treated conditions.  The primary use of our
approach is in basic and discovery research, in which most experiments are
conducted in large cohorts of rodents and interpretation of results are being
translated to higher-order mammals or humans.

# Software Description

# Advantage and Applications

# Acknowledgements

# References