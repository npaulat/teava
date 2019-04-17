# Create list of input files from files in the directory that end in "depth.txt"
x_files <- list.files(pattern="*depth.txt")
# Create list of sample names from by reading the substring of characters [3:len-10) of each file name
sample_names <- substr(x_files, 3, nchar(x_files)-10)

# Define limits of plot, max_coverage evaluated is x, plot axis limit is x+2
coverage_max_to_check <- 70
plot_max <- coverage_max_to_check + 2

# Initialize output lists: 
# 1 = list of lists of coverage depth from ranging 1-50, 
# 2 = list of lists of quotients=(each value in list 1) divided by (sum of the values in list 1)
# 3 = list of mean coverage values and the mean value of the file(?)
output1 <- list()
output2 <- list()
mean_coverage <- c()
# For each file in the x_files list:
for(a in 1:length(x_files)) {
	## Make a list of all data values in the file (DATA_VALUES)
	a_rep <- scan(x_files[a])
	## Initialize local output list (OUTPUT_VALUES)
	a_output <- c()
	## For each coverage value within range of 1 to plot_max:
	for(b in 1:plot_max) {
		### If the value is equal to the plot_max:
		if(b == plot_max) {
			#### Append to local output list the count of all data values >= (plot_max - 1)
			a_output <- c(a_output, length(a_rep[a_rep >= (b - 1)]))
		} else {
			#### Else append to local output list the count of all data values == coverage value
			a_output <- c(a_output, length(a_rep[a_rep == (b - 1)]))
		}
	}
	## Local output #2 is each(?) output #1 value divided by the sum of all values in output #1
	a_output2 <- a_output / sum(a_output)
	## Append local output #1 to output list #1
	output1[[a]] <- a_output
	## Append local output #2 to output list #2
	output2[[a]] <- a_output2
	## Append mean coverage list with the mean of the current file
	mean_coverage <- c(mean_coverage, mean(a_rep))
}

# Remove the last file's list of data values from memory
rm(a_rep)
# Save work (lists) as RData file
save.image("myotis_coverage.RData")


# Make distribution plot of mean coverage vs proportion of genome

# Load work from RData file
load("myotis_coverage.RData")

# Initialize pdf for output plots
pdf("Myotis_coverage_plots.pdf")

# Plot the distribution of coverage vs proportion of ref genome for each individual/sample
par(mfrow=c(3,4))
for(a in 1:length(output2)) {
	plot(0:71, output2[[a]], pch=19, cex=0.1, xlab="Coverage", ylab="Proportion M. lucifugus genome", main=sample_names[a], ylim=c(0,0.10))
	poly.plot <- rbind(cbind(0:71, output2[[a]]), c(71, 0), c(0,0))
	polygon(poly.plot, col="gray")
	abline(v=mean_coverage[a], col="red")
}
dev.off()