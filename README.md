[![DOI](https://zenodo.org/badge/189260858.svg)](https://zenodo.org/badge/latestdoi/189260858)

# Data and code from the avian indel manuscript

Code and data from the paper on analyses of avian indels. 

## Data

### Indel binary encoded data

* [shortname-intronic-locus-indel-files.tar.gz](shortname-intronic-locus-indel-files.tar.gz)
* [uce-locus-indel-files-shortname.tar.gz](uce-locus-indel-files-shortname.tar.gz)


### Gene Trees

The folder [genetrees](genetrees) includes all the gene trees used in the paper. We also include shuffled gene trees used to compute RF among random pair of gene trees. 


### Species Trees

- The folder [species](species) includes all the species trees reported in this study. 
Each folder names corresponds to a tree, and a single newick file is given inside each folder for each species tree. 
This file structure is what DiscoVista takes as input for generating the Meta Table. 
    - See [newModel.txt](newModel.txt) for mapping between folder names and names that appear in the meta table. 
    - See [names-common-2.csv](names-common-2.csv) for mapping to common names 
    - See [clade-defs.txt](clade-defs.txt) and [annotation.txt](annotation.txt) for mapping to clade names used in the meta table
    - See [newTaxa.txt](newTaxa.txt) for mapping between names in clade definition file and names that appear in the meta table. 
- Folder [polytomytest](polytomytest) includes all the trees with the polytomy test applied. 


### Stat files

Stat files used by [draw-figures.R](draw-figures.R) and produced by [commands.sh](commands.sh) are also included here for completeness 
* RF-indel-indel-withlabels.stat, RF-indel-indel.stat, RF-nt-nt-withlabels.stat, RF-nt-nt.stat, and RF.stat (for nti-indel (same locus)): RF distances between pairs of gene trees
* RF.stats: RF distance between pairs of species trees
* qscore.stat: quartet scores of species trees
* nwstats.txt: stats on contracted gene trees
* bootstrap.stats: stats on bootstrap of gene trees



## Code and scripts

### [commands.sh](commands.sh)

This script includes the .nix instructions we used to manipulate data, including, 

* Contracting low support branches
* Combining gene trees
* Running ASTRAL
* Running polytomy test
* Drawing support values on existing trees
* Computing RF distances
* Computing Quartet scores 
* Computing stats on the contraction of gene trees and bootstrap support
* Running DiscoVista

The commands in this file use other tools, including,
* [newick utilities](http://cegg.unige.ch/newick_utils)
* [DiscoVista](https://github.com/esayyari/DiscoVista/)
* [ASTRAL](https://github.com/smirarab/ASTRAL); old version available [here](https://github.com/smirarab/astralhistory)
* scripts from https://github.com/smirarab/global


### [draw-figures.R](draw-figures.R)

Includes scripts used to draw most figures of the paper (using ggplot2 package of R)


### [Indel encoding](README-coding.md) 

Scripts for indel encoding are described in more detail [here](README-coding.md).

### [RAxML commands](README-raxml.md)

RAxML commands used for inferring gene trees from indel data.

### [PAUP commands](PAUP-readme.md)

Scripts for computing RI using  PAUP. 

