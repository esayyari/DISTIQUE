qplot(V8,V12,data=d,color=interaction(V1,V2,V3),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(V7~V9)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)



qplot(X..sites,FN.rate,data=d,color=interaction(method),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(model~X..genes)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)

qplot(X..sites,FN.rate,data=d,color=interaction(method),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(model~X..genes)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)


qplot(X..sites,FN.rate,data=d,color=interaction(method),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(model~X..genes,scales="free_y")+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=5)
ggsave("newres.pdf") 

require(reshape2)
install.packages(c('reshape'),dep=TRUE)
d<-read.csv('filename',header=TRUE/FALSE,sep=',')
The command is recast from the reshape package. 

write.csv(data_wide, file = "wideFN_rates_11-taxon.csv",row.names=FALSE)

data_wide <- spread(d, method,FN_rate )
require(tidyr)
install.packages('tidyr',dep=TRUE)

require(reshape2)
newdata<-dcast(d, model + X_sites + X_taxa + X_genes ~ method, value.var="FN_rate")

d<-read.csv('FN_rates_11-taxon.csv',header=TRUE,sep=",",col.names = c('model','Xtaxa','method','Xsites','Xgenes','replicate','FNrate'))

qplot(Xsites,FNrate,data=d,color=interaction(method),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(model~Xgenes,scales="free_y")+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=5)