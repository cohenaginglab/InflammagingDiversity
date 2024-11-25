# Is inflammaging a universal aging process across human populations?

#### Maximilien Franck, Kamaryn Tanner, Robert L. Tennyson, Camille Daunizeau,   Luigi Ferrucci, Stefania Bandinell,  Benjamin C. Trumble, Hillard S. Kaplan, Jacob E. Aronoff, Thomas S. Kraft, Amanda J. Lea, Vivek V. Venkataraman, Ian J. Wallace, Yvonne A L Lim, Ng Kee Seong, Joe Poe Sheng Yeong, Roger Chun Man Ho,  Ameneh Mehrjerd, Eleftheria Charalambous, Allison E. Aiello, Graham Pawelec, Claudio Franceschi, Johannes Hertel, Tamàs Fülöp, Maël Lemoine, Michael Gurven, and Alan A. Cohen*

This repository contains code to replicate the main analyses and produce figures for this manucript. The datasets used in the manuscript are not able to be shared publicly so we have created some simulated data to illustrate the code. 
* To use this simulated data, copy sampleData1.RDS and sampleData2.RDS to your directory.
* The main code is contained in the file "Figures_2-3-4.R".  This code walks through the factor analysis, generation of biplots, calculation of factor scores, regression of factor scores on simulated health outcomes data and generation of the forest plots shown in the manuscript.
* "Figures_S1-S2-S4.R" includes code for generating the supplementary figures.
* The file "ImputeBLD.R" provides sample code for imputing biomarker values that are below the limit of detection (BLD).
