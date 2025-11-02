## examine results from site reshuffled FST from ANGSD genotypes
## run after angsd_fst_siteshuffle_null.sh

#################
# parameters
#################
minloci <- 10 # should match angsd_fst_siteshuffle_null.r
winsz <- 10000 # window size
winstp <- 2000 # window step

###########################
# load functions
###########################
require(data.table)
#require(plyr)
require(ggplot2)
require(RColorBrewer)

calcp <- function(fst, null) return((sum(null > fst)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen

#####################
# read in and prep data
#####################

# max FST per genome from reshuffling (all sites)
nullAR <- fread('AR_AB.25m.fst.siteshuffle.10.csv.gz')
nullMO <- fread('MO_AB.25m.fst.siteshuffle.10.csv.gz')
nullNJ <- fread('NJ_AB.25m.fst.siteshuffle.10.csv.gz')
nullNY <- fread('NY_AB.25m.fst.siteshuffle.10.csv.gz')



# sliding window FST A and B components from ANGSD, after collapsing to unlinked loci
# header is missing the fst column, so have to skip and make our own
# need to make the all loci AB files
AR <- fread('AR_win_sw2', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) # all sites
KY <- fread('KY_win_50_sw2', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst'))
MO <- fread('MO_win_sw2', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
NJ <- fread('NJ_win_sw2', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
NY <- fread('NY_win_sw2', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 

ARab <- fread('AR_AB_25m.fst.AB.gz') # gatk sites
colnames(ARab) <- c("chr", "POS", "A", "B")
MOab <- fread('MO_AB_25m.fst.AB.gz') 
colnames(MOab) <- c("chr", "POS", "A", "B")
NJab <- fread('NJ_AB_25m.fst.AB.gz')
colnames(NJab) <- c("chr", "POS", "A", "B")
NYab <- fread('NY_AB_25m.fst.AB.gz') 
colnames(NYab) <- c("chr", "POS", "A", "B")

# nucleotide position for the whole genome (start position for each chr)
#in linux
cut -f1,2 mylu.genome.rm.gmap.k100.fa.fai > scaffold_lengths.txt
#back in R
chrmax <- fread('scaffold_lengths.txt')
colnames(chrmax) <- c("chr", "len")
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
chrmax[, chr := sub("^>", "", chr)]

#my alterated way of adding posgen
ARabr <- merge(ARab, unique(chrmax), by = "chr", all.x = TRUE)
ARabr[, posgen := POS + start]

MOabr <- merge(MOab, unique(chrmax), by = "chr", all.x = TRUE)
MOabr[, posgen := POS + start]
NJabr <- merge(NJab, unique(chrmax), by = "chr", all.x = TRUE)
NJabr[, posgen := POS + start]
NYabr <- merge(NYab, unique(chrmax), by = "chr", all.x = TRUE)
NYabr[, posgen := POS + start]
######################
# Calc genome-wide FST
######################

# already removes unplaced
ARabr[!is.na(A), sum(A)/sum(B)]
MOabr[!is.na(A), sum(A)/sum(B)]
NJabr[!is.na(A), sum(A)/sum(B)]
NYabr[!is.na(A), sum(A)/sum(B)]


######################
# Calc windowed FST
######################

# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
  ARabr[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  MOabr[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  NJabr[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
  NYabr[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# calc fst and # snps per window
for(j in 1:(winsz/winstp)){
  if(j ==1){
    ARabrfsts <- ARabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))]
    MOabrfsts <- MOabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))]
    NJabrfsts <- NJabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))]
    NYabrfsts <- NYabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))]
  } 
  if(j > 1){
    ARabrfsts <- rbind(ARabrfsts, ARabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))])
    MOabrfsts <- rbind(MOabrfsts, MOabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))])
    NJabrfsts <- rbind(NJabrfsts, NJabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))])
    NYabrfsts <- rbind(NYabrfsts, NYabr[, .(fst = sum(A)/sum(B), nloci = .N), by = .(chr, midPos = get(paste0('win', j)))])
  } 
}

#remove low nloci my way
AR.10 <- ARabrfsts[nloci >= 10]
MO.10 <- MOabrfsts[nloci >= 10]
NJ.10 <- NJabrfsts[nloci >= 10]
NY.10 <- NYabrfsts[nloci >= 10]

#######################
## null model stats
#######################

# nullcan[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
# nulllof0711[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
# nulllof0714[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
# nulllof1114[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

nullAR[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nullMO[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nullNJ[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nullNY[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

###########################
# calc p-values per window
###########################

# can[, p := calcp(fst, nullcan$x), by = .(chr, midPos)]
# lof0711[, p := calcp(fst, nulllof0711$x), by = .(chr, midPos)]
# lof0714[, p := calcp(fst, nulllof0714$x), by = .(chr, midPos)]
# lof1114[, p := calcp(fst, nulllof1114$x), by = .(chr, midPos)]

AR.10[, p := calcp(fst, nullAR$x), by = .(chr, midPos)]
MO.10[, p := calcp(fst, nullMO$x), by = .(chr, midPos)]
NJ.10[, p := calcp(fst, nullNJ$x), by = .(chr, midPos)]
NY.10[, p := calcp(fst, nullNY$x), by = .(chr, midPos)]


######################
# combine the datasets
######################
# all loci
# lof0711[, pop := 'lof0711']
# lof0714[, pop := 'lof0714']
# lof1114[, pop := 'lof1114']
# can[, pop := 'can']
# 
# dat <- rbind(lof0711, lof0714, lof1114, can)
# nrow(dat)
# dat
# 
# dat[, pop := factor(pop, levels = c('can', 'lof0711', 'lof0714', 'lof1114'))]


# gatk loci
AR.10[, pop := 'AR']
MO.10[, pop := 'MO']
NJ.10[, pop := 'NJ']
NY.10[, pop := 'NY']

datall <- rbind(AR.10, MO.10, NJ.10, NY.10)
nrow(datall)
datall

datall[, pop := factor(pop, levels = c('AR', 'MO', 'NJ', 'NY'))]

# remove NAs and negative windows
datall <- datall[!is.na(fst) & midPos > 0, ]

# sort by window
setkey(datall, pop, chr, midPos)

##############
# Write out
##############

# write.csv(dat, file = gzfile('output/fst_siteshuffle.angsd.csv.gz'), row.names = FALSE)
write.csv(datall, file = gzfile('all.25m.fst_siteshuffle.10.csv.gz'), row.names = FALSE)

##############
# plots
##############

# plot p-value vs. position (all loci)
# cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
# p1 <- ggplot(dat, aes(posgen, -log10(p), color = chr)) + 
#   geom_point(size = 0.2, alpha = 0.3) +
#   facet_wrap(~pop, ncol = 1) +
#   scale_color_manual(values = cols) +
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
# p1
# 
# ggsave(plot = p1, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)
# 


# plot p-value vs. position (gatk loci)
# only where nloci >= minloci
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datall[nloci >= minloci, ], aes(midPos, -log10(p))) + 
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~pop, ncol = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')

ggsave(plot = p2, device = 'png', filename = 'fst.siteshuffle.p_vs_pos.50.png', width = 7.5, height = 6, units = 'in')


# plot p-value vs. nloci (gatk loci)
pn <- ggplot(datall, aes(nloci, -log10(p), color = pop)) + 
  geom_point(size = 0.2, alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') + 
  scale_x_log10()
ggsave(plot = pn, device = 'png', filename = 'fst.siteshuffle.p_vs_nloci.50.png', width = 7.5, height = 6, units = 'in')

#################
# print outliers
#################

library(dplyr)
datall[p < 0.05, ]

datall.out.10 <- filter(datall, p <= 0.9)

write.table(datall.out.10, file = gzfile('all.wind.outliers.loc.10.txt.gz'), row.names = FALSE)

