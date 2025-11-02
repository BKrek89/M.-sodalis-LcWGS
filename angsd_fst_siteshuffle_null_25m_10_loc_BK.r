# shuffle ANGSD FST A and B values across sites and calculate windowed FST to get a null distribution of max genome-wide FST
# run last part of angsd_fst.sh first to get the *.fst.AB.gz files
# if using GATK nodam2 loci, groups into linkage blocks based on ngsLD output. Need to run ngsLD_find_blocks.sh/.r first

# to run on saga

# read command line arguments



# parameters
winsz <- 10000 # window size
winstp <- 2000 # window step
nrep <- 1000 # number of reshuffles
minloci <- 10 # minimum number of loci per window to consider

outfileAR <- 'AR_AB.25m.fst.siteshuffle.10.csv.gz' # used if all loci are used
outfileMO <- 'MO_AB.25m.fst.siteshuffle.10.csv.gz'
outfileNJ <- 'NJ_AB.25m.fst.siteshuffle.10.csv.gz'
outfileNY <- 'NY_AB.25m.fst.siteshuffle.10.csv.gz'


# load functions
require(data.table)

setwd('/local/home/robertk/sodalis/sodalisa1/safnok/25m/model')

#############
# Prep data
#############
# load fst A/B data
AR <- fread('AR_AB_25m.fst.AB.gz')
setnames(AR, c('CHROM', 'POS', 'A', 'B'))

MO <- fread('MO_AB_25m.fst.AB.gz')
setnames(MO, c('CHROM', 'POS', 'A', 'B'))

NJ <- fread('NJ_AB_25m.fst.AB.gz')
setnames(NJ, c('CHROM', 'POS', 'A', 'B'))

NY <- fread('NY_AB_25m.fst.AB.gz')
setnames(NY, c('CHROM', 'POS', 'A', 'B'))








# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
	AR[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	MO[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	NJ[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
	NY[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# mark windows with < minloci for removal
rem <- rep(0,4) # number of windows removed for each of the 4 comparisons
for(j in 1:(winsz/winstp)){
	ARwin <- AR[, .(nsnps = length(POS)), by = .(win = get(paste0('win', j)))] # calc num snps per window
	rem[1] <- rem[1] + ARwin[, sum(nsnps < minloci)] # record number to be removed
	ARwin[, (paste0('win', j, 'keep')) := 1] # create col to mark which windows to keep
	ARwin[nsnps < minloci, (paste0('win', j, 'keep')) := 0] # mark windows to remove
	ARwin[, nsnps := NULL] # drop column
	setnames(ARwin, "win", paste0('win', j)) # change col name
	AR <- merge(AR, ARwin, by = paste0('win', j), all.x = TRUE) # merge keeper col back to full dataset

	MOwin <- MO[, .(nsnps = length(POS)), by = .(win = get(paste0('win', j)))]
	rem[2] <- rem[2] + MOwin[, sum(nsnps < minloci)]
	MOwin[, (paste0('win', j, 'keep')) := 1]
	MOwin[nsnps < minloci, (paste0('win', j, 'keep')) := 0]
	MOwin[, nsnps := NULL]
	setnames(MOwin, "win", paste0('win', j))
	MO <- merge(MO, MOwin, by = paste0('win', j), all.x = TRUE)

	NJwin <- NJ[, .(nsnps = length(POS)), by = .(win = get(paste0('win', j)))]
	rem[3] <- rem[3] + NJwin[, sum(nsnps < minloci)]
	NJwin[, (paste0('win', j, 'keep')) := 1]
	NJwin[nsnps < minloci, (paste0('win', j, 'keep')) := 0]
	NJwin[, nsnps := NULL]
	setnames(NJwin, "win", paste0('win', j))
	NJ <- merge(NJ, NJwin, by = paste0('win', j), all.x = TRUE)

	NYwin <- NY[, .(nsnps = length(POS)), by = .(win = get(paste0('win', j)))]
	rem[4] <- rem[4] + NYwin[, sum(nsnps < minloci)]
	NYwin[, (paste0('win', j, 'keep')) := 1]
	NYwin[nsnps < minloci, (paste0('win', j, 'keep')) := 0]
	NYwin[, nsnps := NULL]
	setnames(NYwin, "win", paste0('win', j))
	NY <- merge(NY, NYwin, by = paste0('win', j), all.x = TRUE)
	
}

rem # number of windows removed for each comparison




####################################
# shuffle and recalc windowed FST
####################################
colnms <- c('CHROM', 'POS', paste0('win', 1:(winsz/winstp)), paste0('win', 1:(winsz/winstp), 'keep')) # list of column names we want out of the base data.table

# CAN
print('Starting Can')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(AR), nrow(AR), replace = FALSE)
	temp <- cbind(AR[, ..colnms], AR[inds, .(A, B)]) # shuffle FSTs across positions
		
	# calc fst for each window to keep
	for(j in 1:(winsz/winstp)){
		temp2 <- temp[get(paste0('win', j, 'keep')) == 1, ] # trim to windows to keep. can't combine with next line for some reason.
		if(j ==1) tempfsts <- temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	# exclude windows with negative midpoints
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfileAR), row.names = FALSE)

rm(maxfst)



# Lof0711
print('Starting Lof0711')
for(i in 1:nrep){
	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(MO), nrow(MO), replace = FALSE)
	temp <- cbind(MO[, ..colnms], MO[inds, .(A, B)]) # shuffle FSTs across positions
	
	# calc fst for each window
	for(j in 1:(winsz/winstp)){
		temp2 <- temp[get(paste0('win', j, 'keep')) == 1, ]
		if(j ==1) tempfsts <- temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfileMO), row.names = FALSE)

rm(maxfst)


# Lof0714
print('Starting Lof0714')
for(i in 1:nrep){
 	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(NJ), nrow(NJ), replace = FALSE)
	temp <- cbind(NJ[, ..colnms], NJ[inds, .(A, B)]) # shuffle FSTs across positions
		
	# calc fst for each window
	for(j in 1:(winsz/winstp)){
		temp2 <- temp[get(paste0('win', j, 'keep')) == 1, ]
		if(j ==1) tempfsts <- temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfileNJ), row.names = FALSE)

rm(maxfst)


# Lof1114
print('Starting Lof1114')
for(i in 1:nrep){
 	cat(i); cat(' ')
	# create new dataset
	inds <- sample(1:nrow(NY), nrow(NY), replace = FALSE)
	temp <- cbind(NY[, ..colnms], NY[inds, .(A, B)]) # shuffle FSTs across positions
		
	# calc fst for each window
	for(j in 1:(winsz/winstp)){
		temp2 <- temp[get(paste0('win', j, 'keep')) == 1, ]
		if(j ==1) tempfsts <- temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
		if(j > 1) tempfsts <- rbind(tempfsts, temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
	}

	# save the max windowed fst
	if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
	if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfileNY), row.names = FALSE)

