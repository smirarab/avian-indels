# RAxML commands for inferring gene trees

##  UCE  Gene Trees


```
for file in *.phy
do
        raxmlHPC-PTHREADS-SSE3_WA -T 4 -s $file -m ASC_BINGAMMA -n $file.bestml.out --asc-corr=lewis --no-seq-check -p 12345 -N 10&

        raxmlHPC-PTHREADS-SSE3_WA -T 8 -s $file -m ASC_BINGAMMA -n $file.bs.out --asc-corr=lewis --no-seq-check -p 12345 -x 12345 -N 200&
done
```

 
## 2) Intron Gene Trees


```
for file in *.phy
do
        raxmlHPC-PTHREADS-SSE3_WA -T 4 -s $file -m ASC_BINGAMMA -n $file.bestml.out --asc-corr=lewis --no-seq-check -p 12345 -N 10&

        raxmlHPC-PTHREADS-SSE3_WA -T 8 -s $file -m ASC_BINGAMMA -n $file.bs.out --asc-corr=lewis --no-seq-check -p 12345 -x 12345 -N 200&

done
```
