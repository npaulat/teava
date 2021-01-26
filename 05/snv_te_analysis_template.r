library(data.table)
library(tictoc)

tic("Read input to memory")

## Read in the species' variant vcf as a tab-separated table
snps <- read.table("<VCF_FILE>", stringsAsFactors=F, sep="\t")

## Read in reference index as a tab-separated table
x <- read.table("<REFINDEX>", sep="\t")
## Make list of scaffolds that are greater than/equal to 1Mb
x.list <- as.character(x[x[,2] >= 1000000,1]) #scaffold list of those > 1 Mbp

## Read in TEs in bed format (Chrom, Start, Stop) as tab-separated table
TEs <- read.table("<BED_FILE>", sep="\t", stringsAsFactors=F)
## Force all data into single data object type
TEs <- cbind(TEs)

## Log statement
log1 <- "Input data read in."
print(log1, quote=FALSE)
toc()

#### Observed values below

tic("Observed SNPs loop")

### Loop for each scaffold to get distance to TE for every SNP
## Initialize output list
output <- list()
## For each scaffold (a) in range of scaffold list
for(a in 1:length(x.list)) {
	## Subset TE file by current scaffold
	a_TEs <- TEs[TEs[,1] %in% x.list[a], ]
	## Set up output TE list
	TE_list <- list()
	## Run the following if there are TEs on this scaffold
	if(nrow(a_TEs) > 0) {
		## Obtain all the bp locations of TEs and remove overlaps
		## For each TE in scaffold TE list
		for(b in 1:nrow(a_TEs)) {
			## Make list of lists by TE index # > each base of the TE
			TE_list[[b]] <- seq(a_TEs[b,2], a_TEs[b,3])
		}
		## Set up data table (a_TEs) of this scaffold's TEs (unique only=remove overlaps)
		a_TEs <- data.table(TEs=unique(unlist(TE_list)))
		## a_TEs with new column 'merge', that contains the info from TEs data table? not scaffold ID, but either TE start or stop (which?)
		a_TEs[,merge:=TEs]
		## Sort the rows of a data.table (a_TEs) by the column (c) 'merge' (now the key)
		setkeyv(a_TEs, c('merge'))
		## Set up this scaffold's SNPs
		## Create data table of SNP positions for each SNP on the scaffold
		a_SNPs <- data.table(SNPs=snps[snps[,1] %in% x.list[a], 2])
		## a_SNPs with new column 'merge', that contains the info from the SNP data table?
		a_SNPs[,merge:=SNPs]
		### Sort the rows of a data.table (a_SNPs) by the column (c) 'merge' (now the key)
		setkeyv(a_SNPs, c('merge'))
		## Identify the nearest TE base pair for each SNP
		merge_SNPs_TEs <- a_TEs[a_SNPs, roll='nearest']
		## Measure the distance from each SNP to the TE
		## A value of zero means the snp is in the TE
		## Output = List of lists by scaffold index # > list of distances
		output[[a]] <- abs(merge_SNPs_TEs$SNPs - merge_SNPs_TEs$TEs)
		print(a)
	}
}
## Vector of all SNP to TE distances
output2 <- unlist(output)
## Table of distances (counts of SNPs at each distance from TE (in list) per scaffold)
output3 <- table(output2)

test1 <- as.numeric(names(output3))
test2 <- as.vector(output3)

## Log statement
log2 <- "Observed SNP to TE data tables complete."
print(log2, quote=FALSE)
toc()


#### Expected values below

tic("Calculate all possible SNP:TE distances loop")

### Loop for each scaffold to get distance to TE for every position in the genome
output_r <- list()
for(a in 1:length(x.list)) {
	## Subset TE file by this scaffold
	a_TEs <- TEs[TEs[,1] %in% x.list[a], ]
	## Set up output file list
	TE_list <- list()
	## Run the following if there are TEs on this scaffold
	if(nrow(a_TEs) > 0) {
		## Obtain all the bp locations of TEs and remove overlaps
		for(b in 1:nrow(a_TEs)) {
			TE_list[[b]] <- seq(a_TEs[b,2], a_TEs[b,3])
		}
		## Set up data table and properties of this scaffold's TEs
		a_TEs <- data.table(TEs=unique(unlist(TE_list)))
		a_TEs[,merge:=TEs]
		setkeyv(a_TEs, c('merge'))
		## Set up data table of every location on the scaffold
		a_SNPs <- data.table(SNPs=seq(from=1, to=x[x[,1] == x.list[a], 2], by=1))
		a_SNPs[,merge:=SNPs]
		setkeyv(a_SNPs, c('merge'))
		## Identify the nearest TE base pair for each position on the scaffold
		merge_SNPs_TEs <- a_TEs[a_SNPs, roll='nearest']
		## Measure the distance from each position to the TE
		## A value of zero means the position is in the TE
		output_r[[a]] <- abs(merge_SNPs_TEs$SNPs - merge_SNPs_TEs$TEs)
		print(a)
	}
}
## Vector of all positions to TE distances
output_r2 <- unlist(output_r)

## Log statement
log3 <- "All possible SNP to TE distances data table generated for expected values loop."
print(log3, quote=FALSE)
toc()

tic("Randomized SNP distribution sampling loop")

## One hundred sampling replicates from the whole genome values with summary
## Also more detailed check of the first five hundred bp
first500 <- list()
replicate_output <- c()
## Modified August 13 2020, to be 500 reps instead of 100 reps ##
#for(a in 1:100) {
for(a in 1:500) {
	## Sample the whole genome by the number of SNPs in the observed dataset
	a_rep <- sample(output_r2, length(output2))
	## Count the proportion that are inside of TEs
	in_TEs <- length(a_rep[a_rep == 0]) / length(a_rep)
	## Summarize the output
	summary_a_rep <- as.vector(summary(a_rep[a_rep > 0]))
	replicate_output <- rbind(replicate_output, c(in_TEs, summary_a_rep[2:5]))
	
	## Table the results
	a_rep <- table(a_rep)
	## Measure the proportion of SNPs at each position (SNPs that are not in TEs)
	## If you want to measure including SNPs in TEs, use the next commented line instead
	#a_rep <- as.vector(a_rep[names(a_rep) %in% seq(from=1, to=500, by=1)]) / sum(as.vector(a_rep[1:length(a_rep)]))
	#a_rep <- as.vector(a_rep[names(a_rep) %in% seq(from=1, to=500, by=1)]) / sum(as.vector(a_rep[2:length(a_rep)]))
	## New line below, line above commented out May 2, 2019 ##
	a_rep <- cbind(as.vector(names(a_rep)[names(a_rep) %in% seq(from=1, to=500, by=1)]), as.vector(a_rep[names(a_rep) %in% seq(from=1, to=500, by=1)]) / sum(as.vector(a_rep[2:length(a_rep)])))
	## Modified the line below May 2, 2019 ##
	first500[[a]] <- data.frame(position=as.numeric(a_rep[,1]), proportion=as.numeric(a_rep[,2]))
	print(a)
}
replicate_output <- data.frame(in_TEs=as.numeric(replicate_output[,1]), Q1=as.numeric(replicate_output[,2]), Median=as.numeric(replicate_output[,3]), Mean=as.numeric(replicate_output[,4]), Q3=as.numeric(replicate_output[,5]))

## Log statement
log4 <- "Expected SNP to TE data tables complete, Q1, Q3, median and mean calculated."
print(log4, quote=FALSE)
toc()

####### Write tables

tic("Write tables")

## Write table of expected values (summary)
write.table(replicate_output, file="<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_expected_values.txt", sep="\t", quote=F, row.names=F, col.names=T)

## Write table of observed values (summary)
write.table(rbind(names(summary(output2)), summary(output2)), file="<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_obs_values.txt", sep="\t",quote=F, row.names=F, col.names=F)

## Log statement
log5 <- "Expected and observed values table files created."
print(log5, quote=FALSE)
toc()

##### Plots

tic("Calculate basic stats for average expected SNP distribution")

## Get the mean for those 500 positions
#output_r3 <- rep(0,500)
#for(a in 1:500) {
#	a_rep <- mean(sapply(first500, "[[", a)) # Can change mean to min or median
#	output_r3[a] <- a_rep
#}

output_r3 <- rep(0,500)
output_r3_min <- c()
output_r3_max <- c()
for(a in 1:500) {
	a_rep <- c()
	for(b in 1:length(first500)) {
		b_rep <- first500[[b]][first500[[b]][,1] %in% a, 2]
		if(length(b_rep) == 1) {
			a_rep <- c(a_rep, b_rep)
		}
	}
	## If some randomizations didn't sample this site, add some zeros to represent those reps
	## Modified  line below August 13 2020, to be 500 reps, not 100 ##
	#a_rep <- c(a_rep, rep(0, (100 - length(a_rep))))
	a_rep <- c(a_rep, rep(0, (500 - length(a_rep))))
	## Added the following two lines June 3, 2019 ##
	output_r3_min <- c(output_r3_min, min(a_rep))
	output_r3_max <- c(output_r3_max, max(a_rep))
	a_rep <- mean(a_rep) # Can change mean to min or median
	output_r3[a] <- a_rep
}

output_close <- rep(0,500)
output3_names <- names(output3)
output3_values <- as.numeric(output3)
output3_sum <- sum(output3_values[2:length(output3_values)])
for(a in 1:500) {
	if(length(output3_values[output3_names %in% a]) > 0) {
		output_close[a] <- output3_values[output3_names %in% a] / output3_sum
	}
}

toc()
tic("Save workspace")

## Discard input tables from memory (SNP table is huge)
rm(snps)
rm(TEs)
rm(x)

## Save R workspace (SNP-TE positions and basic stats
## Important for ability to replot runs
save.image("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>.RData")

tic("Generate plots")

### Plot the <VARIANT> enrichment from TE distance (modify the second abline as necessary) ###
output_close_r <- output_r3
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
pdf("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_enrichment_scaled_plot.pdf")
plot((output_close - output_close_r) * 100, pch=19, cex=0.5, xlab="Distance to TE (bp)", ylab="% Enrichment <VARIANT> Relative to Random", ylim=c(-0.0005, 0.002))
abline(h=0)
dev.off()

### Plot the <VARIANT> enrichment from TE distance, include exp ranges ###
output_close_r <- output_r3
pdf("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_enrichment_ranges_scaled_plot.pdf")
plot((output_close - output_close_r) * 100, pch=19, cex=0.2, xlab="Distance to TE (bp)", ylab="% Enrichment <VARIANT> Relative to Random", bty="l", ylim=c(-0.0005, 0.002))
abline(h=0)
points((output_close - output_close_r) * 100, pch=19, cex=0.2)
rep_range <- c()
for(a in 1:500) {
expected <- c(output_r3_min[a], output_r3_max[a])
expected_y <- (output_close[a] * 100) - (expected * 100)
expected_x <- c(a, a)
lines(expected_x, expected_y, lwd=0.3, col=rgb(0.3,0.3,0.3))
rep_range <- rbind(rep_range, c(a, expected_y[1], expected_y[2]))
}
dev.off()

rep_range <- data.frame(base_distance=as.numeric(rep_range[,1]), Expected_Max=as.numeric(rep_range[,2]), Expected_Min=as.numeric(rep_range[,3]))

write.table(rep_range, file="<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_plot_min_max_values.txt", sep="\t", quote=F, row.names=F, col.names=T)

toc()
tic.clearlog()
