require('ggplot2')
d<-read.csv('b.txt',sep='\t')
wrongBranchPost  <- d[is.na(d$Node2)  , 2 ]
trueBranchesPost <- d[!is.na(d$Node2) , 2 ]


CountFreq<-function(x,v) {sum(v==x)}
trueBranchesPost <- seq(0,1,0.05)
wrongBranchPost <- seq(0,1,0.02)
edges <- seq(0,1,0.1)
r<- 1:(length(edges)-1)


trueIndexc    <- findInterval(trueBranchesPost, edges, rightmost.closed = FALSE, all.inside = TRUE)
trueHistc <- sapply(r,function(rr) CountFreq(rr,v=trueIndexc))
wrongIndexc   <- findInterval(wrongBranchPost,edges,rightmost.closed = FALSE, all.inside = TRUE)
wrongHistc <- sapply(r,function(rr) CountFreq(rr,v=wrongIndexc));
TotalHistc<-wrongHistc+trueHistc
freq <- trueHistc/TotalHistc
M <- matrix(c(edges[1:length(edges)-1],edges[2:length(edges)]),2,length(edges)-1)
Mid<-apply(M,2,function(x) (x[1]+x[2])/2)
cbind(Mid,freq)
qplot(Mid,freq,geom="step")
