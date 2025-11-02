# M.-sodalis-LcWGS
Code for paper "Low-coverage whole-genome analysis of population structure, bottlenecks, and selection in Indiana bats before and after white-nose syndrome" 


This repository contains the bioinformatics and analysis code for the manuscript Low-coverage whole-genome analysis of population structure, bottlenecks, and selection in Indiana bats before and after white-nose syndrome.

This experiment included Indiana bats (Myotis sodalis) from Arkansas, Kentucky, Missouri, New Jersey, and New York. We collected samples after the advent of white-nose syndrome in each state and gathered samples prior to the emergence of WNS from an existing collection at Western Michigan University.

All samples were in the form of wing biopsies or previously extracted DNA.

We prepared sequencing libraries using a reduced volume iteration of the Illumina Nextera XT kit following the methods of Baym et al., 2015, Therkildsen & Palumbi 2017, and Gignoux-Wolfsohn et al., 2021.

Paired-end short read sequences were generated using an Illumina NovaSeq 6000.

Raw sequences are publicly available through NCBI (BioProject: PRJNA1354827)


Recommended order of code files: Bioinformatics, Population Structure, Bottleneck Analysis, Selection.SNPs, Selection.Fst, Selection.Tajima's D, Selection.CMH, Selection.Annotation


The Bioinformatics file contains code for pre- and post-variant filtering, alignment, and variant calling.

The population structure file includes filtering of PCA outlier samples along with the rest of the population structure analysis. The population structure analysis consists of a PCA, admixture analysis, and between state Fst.

The Bottleneck analysis consists of construction of site frequency spectra, calculation of genome-wide Tajima's D and linkage disequilibrium decay.


