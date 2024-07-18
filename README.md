# B2G-22-005
#### Search for pair production of heavy particles decaying to a top quark and a gluon in the lepton+jets final state at $\sqrt{s}=13$ TeV

Code author: F. Labe

### Overview
This repository contains the analysis code of the B2G-22-005 analysis seraching for pair production of heavy particles t\* in the lepton + jets final state. Below, you can find all relevant links an an overview on the code, as well as how to run it.

* [CADI entry](https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=B2G-22-005&tp=an&id=2597&ancode=B2G-22-005)
* [PAS repository](https://gitlab.cern.ch/tdr/notes/B2G-22-005)
* [paper repository](https://gitlab.cern.ch/tdr/papers/B2G-22-005)
* [Analysis note](https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/068)

### Code overview
This analysis is designed to run within the [UHH2](https://github.com/UHH2/UHH2) framework, which is a C-based analysis framework. It runs on custom N-tuples based on MiniAOD, which can be produced using UHH2. More documentation on that can be found in the UHH2 repository Wiki.

To run this analysis, config files can be generated using the [xmlFileCreator.py](https://github.com/finnlabe/TstarTstar/blob/master/config/xmlFileCreator.py "xmlFileCreator.py") script, which creates config files designed to run the analysis.

Four steps (so-called AnalysisModules) are needed to produce the final results: Preselection, Selection, Analysis and DNN. Intermediate files are stored after each step. The DNN step requires a trained DNN, which is provided from outside of the UHH2 framework. Similarly, some required files are produced by the root macros found in `macros/rootmacros`. Some additional AnalysisModules are given as well, which were used for studies during development of the analysis.

Having run the final DNN step of this code, the combine setup documented in [this repository](https://gitlab.cern.ch/cms-analysis/b2g/b2g-22-005/combine_setup) follows to perform the statistical analysis.
