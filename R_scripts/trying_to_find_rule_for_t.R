filt <- read.table("filtered_files/filtered_47.txt",header=F)
superheat((filt), pretty.order.cols=F, pretty.order.rows=F, bottom.label.text.size=2)

delta <- matrix(data=0, ncol=80, nrow=9)
for (i in 1:9) {
    delta[i,] <- as.numeric(filt[i,] - filt[i+1,])
}
superheat((delta), pretty.order.cols=F, pretty.order.rows=F, bottom.label.text.size=2)

segm <- delta
for (i in 1:9) {
 dmed <- median(as.numeric(delta[i,]))
 dmad <- mad(as.numeric(delta[i,]))
 segm[i,] <- ifelse(abs(as.numeric(delta[i,])) > (dmed+dmad/2.0), 1, 0)
}

superheat((segm), pretty.order.cols=F, pretty.order.rows=F, bottom.label.text.size=2)