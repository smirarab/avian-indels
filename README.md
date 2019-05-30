# avian-indels

Code and data from the paper on analyses of avian indels. 

## Scripts

### `commands.sh`

This script includes all the unix instructions we used to manipulate data, including, 

* Contracting low support branches
* Combining gene trees
* Running ASTRAL
* Running polytomy test
* Drawing support values on existing trees
* Computing RF distances
* Computing Quartet scores 
* Computing stats on the contraction of gene trees and bootstrap support
* Running DiscoVista

The instruction use other tools, including,
* newick utilities
* DiscoVista
* ASTRAL
* scripts from https://github.com/smirarab/global


### `draw-figures.R`

Includes scripts used to draw most figures of the paper (using ggplot2 package of R)


## Stat files

Stat files used by `draw-figures.R` and produxed by `commands.sh` are also included here for completeness 



## Species Trees

The folder `species` includes all the species trees reported in this study. Each folder names corresponds to a tree, and a single newick file is given inside each folder for each species tree. 
- This file structure is what DiscoVista takes as input for generating the Meta Table. 
- See [newModel.txt](newModel.txt) for mapping between folder names and names that appear in the meta table. 
- See [names-common-2.csv](names-common-2.csv) for mapping to common names 
- See [clade-defs.txt](clade-defs.txt) and [annotation.txt](annotation.txt) for mapping to clade names used in the meta table

## Gene Trees

The folder [genetrees](genetrees) includes all the gene trees used in the paper. We also include shuffled gene trees used to compute RF among random pair of gene trees. 
