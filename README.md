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

- `clade-defs.txt`: shows the definition of clades used in the Metatable (input to DiscoVista) 
