#!/bin/bash

################################################### Note the files referenced below are scattered in different directories. Sorry for this. You may have to locate them and run the scripts from there or link them here. 

ASTRALDIR=~/workspace/ASTRAL

################################################### Contract low support branches

pushd genetrees
for c in 0 3 5 10 20 33; do 
	nw_ed intron-indel-genetrees.tre 'i & b <= '$c o > intron-indel-genetrees-$c-contract.tre ; echo $c;
done
for c in 0 3 5 10 20 33; do 
	nw_ed intron-nt-genetrees.tre  'i & b <= '$c o > intron-nt-genetrees-$c-contract.tre ; echo $c;
done
for c in 0 3 5 10 20 33; do 
	nw_ed uce-indel-genetrees.tre 'i & b <= '$c o > uce-indel-genetrees-$c-contract.tre ; echo $c;
done
popd genetrees


################################################### Combine: 

pushd genetrees
paste -d '\n' intron-indel-genetrees.tre intron-nt-genetrees.tre >intron-combined-genetrees.tre
for c in 0 3 5 10 20 33; do paste -d '\n' intron-indel-genetrees-$c-contract.tre intron-nt-genetrees-$c-contract.tre >intron-combined-genetrees-$c-contract.tre ;  done
paste -d '\n' uce-indel-genetrees.tre uce-nt-genetrees.tre > uce-combined-genetrees.tre
paste -d '\n' uce-indel-genetrees-5-contract.tre uce-nt-genetrees-5-contract.tre > uce-combined-genetrees-5-contract.tre
ls *5-contr*|grep -v comb |xargs cat > combined-combined-genetrees-5-contract.tre
popd genetrees


################################################### Run ASTRAL

java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-indel-genetrees.tre -o intron-nt-astraltree.tre 2> intron-indel-astraltree.log &
for c in 0 3 5 10 20 33; do 
	java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-indel-genetrees-$c-contract.tre -o intron-indel-astraltree-$c-contract.tre 2> intron-indel-astraltree-$c-contract.log &
done

java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-nt-genetrees.tre -o intron-nt-astraltree.tre 2> intron-nt-astraltree.log &
for c in 0 3 5 10 20 33; do 
	java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-nt-genetrees-$c-contract.tre -o intron-nt-astraltree-$c-contract.tre 2> intron-nt-astraltree-$c-contract.log &
done

java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i uce-indel-genetrees.tre -o uce-indel-astraltree.tre 2> uce-indel-astraltree.log &
for c in 0 3 5 10 20 33; do 
	java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i uce-indel-genetrees-$c-contract.tre -o uce-indel-astraltree-$c-contract.tre 2> uce-indel-astraltree-$c-contract.log &
done


#for c in 0 3 5 10 20 33; do
#	java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-combined-genetrees-$c-contract.tre -o intron-combined--astraltree-$c-contract.tre 2> intron-combined--astraltree-$c-contract.log & 
#done

java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-combined-genetrees-5-contract.tre -o intron-combined-astraltree-5-contract.tre -w2 2> intron-combined-astraltree-5-contract.log  ## Note; this was originally run by mistake with no -w2. Mistake was caught and fixed later. Be aware. 
java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i uce-combined-genetrees-5-contract.tre -o uce-combined-astraltree-5-contract.tre -w2 2> uce-combined-astraltree-5-contract.log
java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i combined-combined-genetrees-5-contract.tre -o combined-combined-asltraltree-5-contract.tre -w2 2> combined-combined-asltraltree-5-contract.log & 

################################################### score for polytomies

for d in intron-indel/intron-indel intron-nt/intron-nt uce-indel/uce-indel uce-nt/uce-nt; do 
	for c in 0 3 5 10 20 33; do 
		java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i $d-genetrees-$c-contract.tre -q $d-astraltree-$c-contract.tre -t10  -o $d-astraltree-$c-contract-polytometest.tre 2> $d-astraltree-$c-contract-polytometest.log &
	done;
	java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i $d-genetrees.tre -q $d-astraltree.tre -o $d-astraltree-polytometest.tre -t10 2> $d-astraltree.polytometest.log &
done
java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-indel-genetrees.tre -q intron-indel-astraltree.tre -o intron-indel-astraltree-polytometest.tre -t10 2> intron-indel-astraltree.polytometest.log &
java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i intron-nt-genetrees.tre -q intron-nt-astraltree.tre -o intron-nt-astraltree-polytometest.tre -t10 2> intron-nt-astraltree.polytometest.log &

################################################### score for quartets

for d in intron-indel/intron-indel intron-nt/intron-nt uce-indel/uce-indel uce-nt/uce-nt; do 
	for c in 0 3 5 10 20 33; do 
		java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i $d-genetrees-$c-contract.tre -q $d-astraltree-$c-contract.tre -t1  -o $d-astraltree-$c-contract-quartetscore.tre 2> $d-astraltree-$c-contract-quartetscore.log &
	done;
	java -jar $ASTRALDIR/Astral/astral.5.6.1.jar -i $d-genetrees.tre -q $d-astraltree.tre -o $d-astraltree-quartetscore.tre -t1 2> $d-astraltree.quartetscore.log &
done

################################################### Cleanup logs
sed -i '/Warning:/d' *.log    # astral.5.6.1 adds a spurious warning that should be removed; it makes the log files long otherwise



################################################### RF distance

# allastraltrees.txt has the path to all astral trees. 
# For example, if you ran commands above, run
#      ls *astraltree.tre > allastraltrees.txt
# or if you are using species trees we provide in DiscoVista format, run
#      ls species/*/*.tree > allastraltrees.txt
for x in `cat allastraltrees.txt`; do 
	for y in `cat allastraltrees.txt`; do 
		echo -n $x $y" "; compareTrees.missingBranch $x $y; echo; 
	done; 
done |tee RF.stats.txt
echo "V1 V2 V3 V4 V5" > RF.stats; cat RF.stats.txt |sed -e "s/-astraltree//g" -e "s/.tre//g" -e "s/\([^/ ]*\)\///g" -e "s/-$//g" -e "s/-\([0-9]\)-/-0\1-/g" -e "s/-contract/%/g" |sed '/^$/d'|tee -a RF.stats


################################################### RF among gene trees

# Create Shuffled gene trees
perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' < uce-indel-genetrees.tre > shuffled-uce-indel-genetrees.tre
perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' < uce-nt-genetrees.tre > shuffled-uce-nt-genetrees.tre
perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' < intron-indel-genetrees.tre > shuffled-intron-indel-genetrees.tre
perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' < intron-nt-genetrees.tre > shuffled-intron-nt-genetrees.tre

compareTrees genetrees/intron-indel-genetrees.tre shuffled-intron-indel-genetrees.tre -simplify 1>/dev/null  2> RF-indel-indel.stat
compareTrees genetrees/intron-nt-genetrees.tre shuffled-intron-nt-genetrees.tre -simplify 1>/dev/null  2> RF-nt-nt.stat
compareTrees genetrees/intron-indel-genetrees.tre genetrees/intron-nt-genetrees.tre 1>/dev/null  2> RF.stat  ## Sorry for the bad name !

compareTrees genetrees/uce-indel-genetrees.tre genetrees/shuffled-uce-indel-genetrees.tre -simplify  1>/dev/null 2> uce-RF-indel-indel.stat 
compareTrees genetrees/uce-nt-genetrees.tre genetrees/shuffled-uce-nt-genetrees.tre -simplify  1>/dev/null 2> uce-RF-nt-nt.stat 
compareTrees genetrees/uce-indel-genetrees.tre genetrees/uce-nt-genetrees.tre 1>/dev/null  2> uce-RF.stat 


################################################### Quartet scores

# See "RF distance" for how to get allastraltrees.txt
#
cat allastraltrees.txt|sed -e "s/.tre$/.log/" |xargs -I@ grep Normalized @|tee quartetscores.stat
cat quartetscores.stat |sed -e "s:.*/::g" -e "s/-/ /g" -e "s/astraltree//g" -e "s/contract.log.*: //g" -e "s/.log.*: /no-contract /g" -e "s/  */ /g" > qscore.stat



################################################### Draw trees

# These are not figures used in the paper. In the paper, we use FigTree. These just didn't look nice enough. 
echo *tree-polytometest.tre  *tree-?-contract-polytometest.tre *tree-??-contract-polytometest.tre  |grep -v nexus |xargs cat|nw_topology -|nw_order -cn - | nw_ed - 'i & b>0.05' o| nw_rename - names-common-2.csv   | draw-trees.sh - poly-contracted.pdf 1200
echo *tree.tre  *tree-?-contract.tre *tree-??-contract.tre  |grep -v nexus |xargs cat|nw_topology -|nw_order -cn - | nw_rename - names-common-2.csv   | draw-trees.sh - allastraltrees.pdf 1200
echo *tree-quartetscore.tre  *tree-?-contract-quartetscore.tre *tree-??-contract-quartetscore.tre  |grep -v nexus |xargs cat |nw_reroot -s - ANAPL|nw_reroot -s - STRCA TINMA| nw_topology -|nw_order -cn - | nw_rename - names-common-2.csv   | draw-trees.sh - $(basename `pwd`)-quartetscore.pdf 1200


################################################### Compute contraction stats
for x in 0 3 5 10 20 33; do cat <( cat uce-indel-genetrees-$x-contract.tre|grep STRCA|nw_reroot - STRCA| nw_stats -fl - ) <( cat uce-indel-genetrees-$x-contract.tre|grep -v STRCA|nw_reroot - TINMA| nw_stats -fl - )  |tee uce-indel-genetrees-$x-contract.nw.stats|wc -l; done
cat uce-indel-genetrees.tre |grep STRCA|nw_reroot - STRCA| nw_stats -fl - |tee uce-indel-genetrees.nw.stats|wc
for x in 0 3 5 10 20 33; do cat <( cat uce-nt-genetrees-$x-contract.tre|grep STRCA|nw_reroot - STRCA| nw_stats -fl - ) <( cat uce-nt-genetrees-$x-contract.tre|grep -v STRCA|nw_reroot - TINMA| nw_stats -fl - )  |tee uce-nt-genetrees-$x-contract.nw.stats|wc -l; done
cat uce-nt-genetrees.tre |grep STRCA|nw_reroot - STRCA| nw_stats -fl - |tee uce-nt-genetrees.nw.stats|wc
grep "gram" *nw.stats|gsed -e "s/:/\t/g"|gsed -e "s/-genetrees./\t/g" -e "s/-contr.*stats//g" -e "s/nw.stats/no-contract/g" |tee nwstats.txt|wc -l



################################################### Bootstrap summary
bootstrap_summary.sh genetrees/intron-nt-genetrees.tre 2>/dev/null > intron-nt-genetrees-bs.stat
bootstrap_summary.sh genetrees/intron-indel-genetrees.tre 2>/dev/null > intron-indel-genetrees-bs.stat
bootstrap_summary.sh genetrees/uce-indel-genetrees.tre 2>/dev/null > uce-indel-genetrees-bs.stat
bootstrap_summary.sh genetrees/uce-nt-genetrees.tre 2>/dev/null > uce-nt-genetrees-bs.stat
grep " " *bs.stat|sed -e "s/.*\///g" -e "s/-genetrees.tre//g" > bootstrap.stats



################################################### DiscoVista
docker run -v `pwd`/:/data esayyari/discovista discoVista.py -c clade-defs.txt -p species/ -t 0.95 -y newModel.txt -w newTaxa.txt -m 0 -o results/



################################################### Draw support values on RAxML indel gene trees
# TODO: this will not run unless we provide original RAxML files. Too big currently for github (bu available on google drive). 
for x in `cat loci`; do ./raxmlHPC -fb -m GTRCAT -t RAxML_bestTree.$x.indels.phy.bestml.out -z RAxML_bootstrap.$x.indels.phy.bs.out  -n $x.indels.phy.bp.out ; done
