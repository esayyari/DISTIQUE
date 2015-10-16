d<-read.csv('mammalian-avian-res.csv',sep=",",header=FALSE)
require('ggplot2')
levels(d$V6)<-list("500bp"="500","True gene trees"="true","1000bp"="1000","1500bp"="1500","250bp"="250")
c<-d$V4 %in% c("0.2X")
b<-d$V1 %in% c("mammalian")
a<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm","distique-cons-prod","distique-cons-min")
g<-d$V5 %in% c(200)
h<-d$V6 %in% c("500bp","True gene trees")
e<-d$V2 %in% c("distique-2-prod","astrid","astral")
f<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm")
i<-d$V1 %in% c("avian")
j<-d$V4 %in% c("1X")
l<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm","distique-cons-prod","distique-cons-min")
m<-d$V5 %in% c(1000)
n<-d$V6 %in% c("500bp","True gene trees")
o<-d$V2 %in% c("distique-2-prod","astrid","astral")
p<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm")
 
 #d$V5<-as.factor(d$V5)
qplot(V5,V8,data=d[a & b & c & h,],color=interaction(V2),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","Distique-rndsample-min","Distique-rndsample-all-pairs","all-pairs"))
ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique-vt.pdf')

qplot(V4,V8,data=d[a & b & g & h,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","Distique-rndsample-min","Distique-rndsample-all-pairs","all-pairs"))
ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-vt.pdf')
qplot(V4,V8,data=d[e & b & g & h,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique-all-pairs"))
ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-Vs-Other.pdf')

qplot(V5,V8,data=d[e & b & c,],color=interaction(V2),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique-all-pairs"))
ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique-Vs-others.pdf')
 

qplot(V5,V8,data=d[f & b & c,],color=interaction(V2),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","all-pairs"))
ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique.pdf')

qplot(V4,V8,data=d[f & b & g & h,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","all-pairs"))

ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique.pdf')

qplot(V5,V8,data=d[l & i & j & n,],color=interaction(V2),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","Distique-rndsample-min","Distique-rndsample-all-pairs","all-pairs"))
ggsave('Avian-1ILS-genetreesVSFNrates-Distique-vt.pdf')


qplot(V4,V8,data=d[l & m & i & n,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","Distique-rndsample-min","Distique-rndsample-all-pairs","all-pairs"))
ggsave('Avian-1000-genes-ILS-Vs-FNrate-Distique-vt.pdf')

qplot(V5,V8,data=d[j & i & o & n,],color=interaction(V2),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique-all-pairs"))
ggsave('Avian-1ILS-genes-ILS-Vs-FNrate-Distique-Vs-Other.pdf')

qplot(V4,V8,data=d[o & m & i & n,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique-all-pairs"))
ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Other.pdf')

qplot(V5,V8,data=d[j & i & p & n,],color=interaction(V2),geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.1)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","all-pairs"))
ggsave('Avian-1ILS-genetreesVSFNrates-Distique.pdf',width=4.5,height=4.5)

qplot(V4,V8,data=d[i & p & n & m,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Distique-min","Distique-all-pairs","all-pairs"))

ggsave('Avian-1000-genes-ILS-Vs-FNrate-Distique.pdf')