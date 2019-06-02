# PAUP README


Each nexus file and tree file were located in their own directory to avoid overlapping files before running the script.

The PAUP command used to generate ensemble RI values is as follows:

```
begin paup; 
set autoclose=yes warntree=no warnreset=no; 
log start file=out.log replace; 
gettree file=tree.txt; 
Pscores / RI scoreFile=scores.txt; 
quit; 
end;
```

Where `tree.txt` is the generic gene tree name matching the locus nexus file it is paired with. (i.e. for locus 9999 the `tree.txt` file would represent the gene tree inferred for locus 9999)

Where `scores.txt` is the generic name given to the ensemble RI score later to be renamed to its assigned locus. (see above)


