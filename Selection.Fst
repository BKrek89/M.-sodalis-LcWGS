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


##########################################################99th percentile and hypergeometric test
##in R. fst windows. find 99th percentile
# Read input data
data.ar2510 <- read.table("AR.25m.win.10.cut.nh.n", header = FALSE)
data.mo2510 <- read.table("MO.25m.win.10.cut.nh.n", header = FALSE)
data.nj2510 <- read.table("NJ.25m.win.10.cut.nh.n", header = FALSE)
data.ny2510 <- read.table("NY.25m.win.10.cut.nh.n", header = FALSE)
data.ar2550 <- read.table("AR.25m.win.50.cut.nh.n", header = FALSE)
data.mo2550 <- read.table("MO.25m.win.50.cut.nh.n", header = FALSE)
data.nj2550 <- read.table("NJ.25m.win.50.cut.nh.n", header = FALSE)
data.ny2550 <- read.table("NY.25m.win.50.cut.nh.n", header = FALSE)

data.ar7510 <- read.table("AR.75m.win.10.cut.nh.n", header = FALSE)
data.mo7510 <- read.table("MO.75m.win.10.cut.nh.n", header = FALSE)
data.nj7510 <- read.table("NJ.75m.win.10.cut.nh.n", header = FALSE)
data.ny7510 <- read.table("NY.75m.win.10.cut.nh.n", header = FALSE)
data.ar7550 <- read.table("AR.75m.win.50.cut.nh.n", header = FALSE)
data.mo7550 <- read.table("MO.75m.win.50.cut.nh.n", header = FALSE)
data.nj7550 <- read.table("NJ.75m.win.50.cut.nh.n", header = FALSE)
data.ny7550 <- read.table("NY.75m.win.50.cut.nh.n", header = FALSE)

# Calculate the 99th percentile of column 4
threshold.ar2510 <- quantile(data.ar2510$V4, 0.99)
threshold.mo2510 <- quantile(data.mo2510$V4, 0.99)
threshold.nj2510 <- quantile(data.nj2510$V4, 0.99)
threshold.ny2510 <- quantile(data.ny2510$V4, 0.99)
threshold.ar2550 <- quantile(data.ar2550$V4, 0.99)
threshold.mo2550 <- quantile(data.mo2550$V4, 0.99)
threshold.nj2550 <- quantile(data.nj2550$V4, 0.99)
threshold.ny2550 <- quantile(data.ny2550$V4, 0.99)

threshold.ar7510 <- quantile(data.ar7510$V4, 0.99)
threshold.mo7510 <- quantile(data.mo7510$V4, 0.99)
threshold.nj7510 <- quantile(data.nj7510$V4, 0.99)
threshold.ny7510 <- quantile(data.ny7510$V4, 0.99)
threshold.ar7550 <- quantile(data.ar7550$V4, 0.99)
threshold.mo7550 <- quantile(data.mo7550$V4, 0.99)
threshold.nj7550 <- quantile(data.nj7550$V4, 0.99)
threshold.ny7550 <- quantile(data.ny7550$V4, 0.99)

# Filter rows based on the threshold
filtered_data.ar2510 <- data.ar2510[data.ar2510$V4 >= threshold.ar2510, ]
filtered_data.mo2510 <- data.mo2510[data.mo2510$V4 >= threshold.mo2510, ]
filtered_data.nj2510 <- data.nj2510[data.nj2510$V4 >= threshold.nj2510, ]
filtered_data.ny2510 <- data.ny2510[data.ny2510$V4 >= threshold.ny2510, ]
filtered_data.ar2550 <- data.ar2550[data.ar2550$V4 >= threshold.ar2550, ]
filtered_data.mo2550 <- data.mo2550[data.mo2550$V4 >= threshold.mo2550, ]
filtered_data.nj2550 <- data.nj2550[data.nj2550$V4 >= threshold.nj2550, ]
filtered_data.ny2550 <- data.ny2550[data.ny2550$V4 >= threshold.ny2550, ]

filtered_data.ar7510 <- data.ar7510[data.ar7510$V4 >= threshold.ar7510, ]
filtered_data.mo7510 <- data.mo7510[data.mo7510$V4 >= threshold.mo7510, ]
filtered_data.nj7510 <- data.nj7510[data.nj7510$V4 >= threshold.nj7510, ]
filtered_data.ny7510 <- data.ny7510[data.ny7510$V4 >= threshold.ny7510, ]
filtered_data.ar7550 <- data.ar7550[data.ar7550$V4 >= threshold.ar7550, ]
filtered_data.mo7550 <- data.mo7550[data.mo7550$V4 >= threshold.mo7550, ]
filtered_data.nj7550 <- data.nj7550[data.nj7550$V4 >= threshold.nj7550, ]
filtered_data.ny7550 <- data.ny7550[data.ny7550$V4 >= threshold.ny7550, ]

# Write the filtered data to a new file
write.table(filtered_data.ar2510, "AR.25m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.mo2510, "MO.25m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.nj2510, "NJ.25m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.ny2510, "NY.25m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.ar2550, "AR.25m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.mo2550, "MO.25m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.nj2550, "NJ.25m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.ny2550, "NY.25m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(filtered_data.ar7510, "AR.75m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.mo7510, "MO.75m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.nj7510, "NJ.75m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.ny7510, "NY.75m.10.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.ar7550, "AR.75m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.mo7550, "MO.75m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.nj7550, "NJ.75m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_data.ny7550, "NY.75m.50.99p.fst", sep = "\t", quote = FALSE, row.names = FALSE)

###in linux. find number shared
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.25m.10.99p.fst MO.25m.10.99p.fst > ARMO.25m.10.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.25m.10.99p.fst.match NJ.25m.10.99p.fst > ARMONJ.25m.10.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.25m.10.99p.fst.match NY.25m.10.99p.fst > ARMONJNY.25m.10.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.25m.50.99p.fst MO.25m.50.99p.fst > ARMO.25m.50.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.25m.50.99p.fst.match NJ.25m.50.99p.fst > ARMONJ.25m.50.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.25m.50.99p.fst.match NY.25m.50.99p.fst > ARMONJNY.25m.50.99p.fst.match

awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.75m.10.99p.fst MO.75m.10.99p.fst > ARMO.75m.10.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.75m.10.99p.fst.match NJ.75m.10.99p.fst > ARMONJ.75m.10.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.75m.10.99p.fst.match NY.75m.10.99p.fst > ARMONJNY.75m.10.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.75m.50.99p.fst MO.75m.50.99p.fst > ARMO.75m.50.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.75m.50.99p.fst.match NJ.75m.50.99p.fst > ARMONJ.75m.50.99p.fst.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.75m.50.99p.fst.match NY.75m.50.99p.fst > ARMONJNY.75m.50.99p.fst.match


#find total number of shared snps
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.25m.win.10.cut.nh.n MO.25m.win.10.cut.nh.n > ARMO.25m.10.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.25m.10.fst.tot.match NJ.25m.win.10.cut.nh.n > ARMONJ.25m.10.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.25m.10.fst.tot.match NY.25m.win.10.cut.nh.n > ARMONJNY.25m.10.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.25m.win.50.cut.nh.n MO.25m.win.50.cut.nh.n > ARMO.25m.50.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.25m.50.fst.tot.match NJ.25m.win.50.cut.nh.n > ARMONJ.25m.50.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.25m.50.fst.tot.match NY.25m.win.50.cut.nh.n > ARMONJNY.25m.50.fst.tot.match

awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.75m.win.10.cut.nh.n MO.75m.win.10.cut.nh.n > ARMO.75m.10.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.75m.10.fst.tot.match NJ.75m.win.10.cut.nh.n > ARMONJ.75m.10.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.75m.10.fst.tot.match NY.75m.win.10.cut.nh.n > ARMONJNY.75m.10.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' AR.75m.win.50.cut.nh.n MO.75m.win.50.cut.nh.n > ARMO.75m.50.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMO.75m.50.fst.tot.match NJ.75m.win.50.cut.nh.n > ARMONJ.75m.50.fst.tot.match
awk 'FNR==NR{arr[$1,$2];next}(($1,$2) in arr)' ARMONJ.75m.50.fst.tot.match NY.75m.win.50.cut.nh.n > ARMONJNY.75m.50.fst.tot.match


### usually hypergeometric test next but no shared regions. 



##############################################################################higher Fst than expected given null model. All code adapted from Pinsky et al., 2021
for POP in AR MO NJ NY; do echo $POP;/local/opt/angsd/misc/realSFS fst print "$POP"_AB_75m.fstd.fst.idx | gzip > model/"$POP"_AB_75m.fst.AB.gz
done

#for construction of null model and comparison, see attached R files angsd_fst_siteshuffle_null_stats_25m_10_BK and angsd_fst_siteshuffle_null_stats_25m_10_BK in that order
