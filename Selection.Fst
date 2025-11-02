##variant call
for POP in AR_B MO_B NJ_B NY_B; do echo $POP;/local/opt/angsd/angsd -b "$POP".25m.imb.pco3.txt -nThreads 20 -anc /local/home/robertk/sodalis/sodalisa1/genome/mylu.genome.rm.gmap.k100.fa -GL 1 -dosaf 1 
-minMapQ 30 -minQ 20 -rf /local/home/robertk/sodalis/sodalisa1/sites/chrs.imb.txt -sites /local/home/robertk/sodalis/sodalisa1/sites/skno.new.imb.sites1.sorted.txt 
-out /local/home/robertk/sodalis/sodalisa1/safnok/25m/"$POP"_25m_fst; done

/local/opt/angsd/misc/realSFS -sites /local/home/robertk/sodalis/sodalisa1/sites/skno.new.imb.sites1.sorted.txt AR_A_75m_fst.saf.idx AR_B_75m_fst.saf.idx -nSites 6000000 -P 20 > AR_AB_75m.fst.ml

for POP in AR MO NJ NY; do echo $POP;
awk 'BEGIN{FS=OFS=" "}
     NR>=1{for (i=1;i<=NF;i++) a[i]+=$i}
     END{for (i=1;i<=NF;i++) printf a[i] OFS; printf "\n"}' "$POP"_AB_75m.fst.ml > "$POP"_AB_75m.fst.sum.ml
done 

########Genome-wide Fst between pre- and post-WNS groups
/local/opt/angsd/misc/realSFS fst index AR_A_75m_fst.saf.idx AR_B_75m_fst.saf.idx -sfs AR_AB_75m.fst.sum.ml -fstout AR_AB_75m.fstd

/local/opt/angsd/misc/realSFS fst stats AR_AB_75m.fstd.fst.idx

##########10,000 and 50,0000 bp windows with Fst > 5 standard deviations of the mean
#Sliding window for FST outlier regions
#50000, 10000 windows
for POP in AR MO NJ NY; do echo $POP;/local/opt/angsd/misc/realSFS fst stats2 "$POP"_AB_75m.fstd.fst.idx -win 50000 -step 10000 -type 0 > "$POP".75m.win.50.sw2
done

for POP in AR MO NJ NY; do echo $POP;awk '{print $2" "$3" "$4" "$5}' "$POP".75m.win.50.sw2 > "$POP".75m.win.50.cut
done 

for POP in AR MO NJ NY; do echo $POP;tail -n+2 "$POP".75m.win.50.cut > "$POP".75m.win.50.cut.nh
done

#average window fst
awk '{ sum += $4 } END { if (NR > 0) print sum/NR }' AR.25m.win.50.cut.nh

#standard dev window fst
awk '{ sum += $4; sumsq += ($4)^2 } END { if (NR > 1) print sqrt((sumsq - sum^2/NR) / (NR-1)) }' AR.25m.win.10.cut.nh

awk '($4 > 0.035007) ' AR.25m.win.10.cut.nh > 5std/AR.fst.10.5sd

awk '($3 > 10) ' AR.25m.fst.10.5sd > AR.25m.fst.10.5sd.n

#look for same regions
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.25m.fst.10.5sd.n MO.25m.fst.10.5sd.n > ARMO.10.fst

#10000 2000 windows
for POP in AR MO NJ NY; do echo $POP;/local/opt/angsd/misc/realSFS fst stats2 "$POP"_AB_25m.fstd.fst.idx -win 10000 -step 2000 -type 0 > "$POP".25m.win.10.sw2
done


for POP in AR MO NJ NY; do echo $POP;awk '{print $2" "$3" "$4" "$5}' "$POP".25m.win.10.sw2 > "$POP".25m.win.10.cut
done 

for POP in AR MO NJ NY; do echo $POP;tail -n+2 "$POP".25m.win.10.cut > "$POP".25m.win.10.cut.nh
done

#average window fst
awk '{ sum += $4 } END { if (NR > 0) print sum/NR }' AR.25m.win.10.cut.nh

#standard dev window fst
awk '{ sum += $4; sumsq += ($4)^2 } END { if (NR > 1) print sqrt((sumsq - sum^2/NR) / (NR-1)) }' AR.25m.win.10.cut.nh

awk '($4 > 0.02173534) ' NY.25m.win.10.cut.nh > NY.25m.10.5sd
awk '($3 > 10) ' AR.25m.10.5sd > AR.25m.10.5sd.n

#look for same regions
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.25m.fst.10.5sd.n MO.25m.fst.10.5sd.n > ARMO.10.mat


