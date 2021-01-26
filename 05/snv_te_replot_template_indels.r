library(data.table)

load(file = "<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>.RData")
ls()

########

### Plot the <VARIANT> enrichment from TE distance (modify the second abline as necessary) ###
output_close_r <- output_r3
# next line commented out (June 3)
#output_close <- as.vector(output3[names(output3) %in% seq(from=1, to=500, by=1)]) / sum(as.vector(output3[2:length(output3)]))
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
pdf("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_enrichment_scaled_plot.pdf")
plot((output_close - output_close_r) * 100, pch=19, cex=0.5, xlab="Distance to TE (bp)", ylab="% Enrichment <VARIANT> Relative to Random", ylim=c(-0.0005, 0.001), yaxt="n")
marks <- c(-0.0005,0.0000,0.0005,0.0010)
axis(2,at=marks,labels=format(marks,scientific=FALSE))
abline(h=0)
dev.off()

svg("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_enrichment_scaled_plot.svg")
plot((output_close - output_close_r) * 100, pch=19, cex=0.5, xlab="Distance to TE (bp)", ylab="% Enrichment <VARIANT> Relative to Random", ylim=c(-0.0005, 0.001), yaxt="n")
marks <- c(-0.0005,0.0000,0.0005,0.0010)
axis(2,at=marks,labels=format(marks,scientific=FALSE))
abline(h=0)
dev.off()

### Plot the <VARIANT> enrichment from TE distance, include exp ranges ###
output_close_r <- output_r3
# next line commented out (June 3)
#output_close <- as.vector(output3[names(output3) %in% seq(from=1, to=500, by=1)]) / sum(as.vector(output3[2:length(output3)]))
pdf("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_enrichment_ranges_scaled_plot.pdf")
plot((output_close - output_close_r) * 100, pch=19, cex=0.2, xlab="Distance to TE (bp)", ylab="% Enrichment <VARIANT> Relative to Random", bty="l", ylim=c(-0.0005, 0.001), yaxt="n")
marks <- c(-0.0005,0.0000,0.0005,0.0010)
axis(2,at=marks,labels=format(marks,scientific=FALSE))
abline(h=0)
points((output_close - output_close_r) * 100, pch=19, cex=0.2)
rep_range <- c()
for(a in 1:500) {
# changed the next line on June 3
expected <- c(output_r3_min[a], output_r3_max[a])
expected_y <- (output_close[a] * 100) - (expected * 100)
expected_x <- c(a, a)
lines(expected_x, expected_y, lwd=0.3, col=rgb(0.3,0.3,0.3))
rep_range <- rbind(rep_range, c(a, expected_y[1], expected_y[2]))
}
dev.off()

svg("<OUTDIR>/<DATASET>_<SPECIES>_<TE_TYPE>_<IN_TYPE>_<VTYPE>_enrichment_ranges_scaled_plot.svg")
plot((output_close - output_close_r) * 100, pch=19, cex=0.2, xlab="Distance to TE (bp)", ylab="% Enrichment <VARIANT> Relative to Random", bty="l", ylim=c(-0.0005, 0.001), yaxt="n")
marks <- c(-0.0005,0.0000,0.0005,0.0010)
axis(2,at=marks,labels=format(marks,scientific=FALSE))
abline(h=0)
points((output_close - output_close_r) * 100, pch=19, cex=0.2)
rep_range <- c()
for(a in 1:500) {
# changed the next line on June 3
expected <- c(output_r3_min[a], output_r3_max[a])
expected_y <- (output_close[a] * 100) - (expected * 100)
expected_x <- c(a, a)
lines(expected_x, expected_y, lwd=0.3, col=rgb(0.3,0.3,0.3))
rep_range <- rbind(rep_range, c(a, expected_y[1], expected_y[2]))
}
dev.off()
