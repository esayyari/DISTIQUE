d<-read.csv('mammalian-avian-res.csv',sep=",",header=FALSE)
require('ggplot2')
levels(d$V6)<-list("500bp"="500","True gene trees"="true","1000bp"="1000","1500bp"="1500","250bp"="250")
c<-d$V4 %in% c("0.2X")
b<-d$V1 %in% c("mammalian")
a<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme","distique-cons-prod","distique-cons-min")
g<-d$V5 %in% c(200)
h<-d$V6 %in% c("500bp","True gene trees")
e<-d$V2 %in% c("distique-2-prod","astrid","astral")
f<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme")
i<-d$V1 %in% c("avian")
j<-d$V4 %in% c("1X")
l<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme","distique-cons-prod","distique-cons-min")
m<-d$V5 %in% c(1000)
n<-d$V6 %in% c("500bp","True gene trees")
o<-d$V2 %in% c("distique-2-prod","astrid","astral")
p<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme")
 
 #d$V5<-as.factor(d$V5)
qplot(V5,V8,data=d[a & b & c & h,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Average (Cons.)","Min (Cons.)","Min (Cons.-sampling)","Average (Cons.-sampling)","Average"))+coord_cartesian(ylim=c(-0.003,0.303))
ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique-vt.pdf',width=4.5,height=4.5)

qplot(V4,V8,data=d[a & b & g & h,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Average (Cons.)","Min (Cons.-sampling)","Average (Cons.-sampling)","Average"))+coord_cartesian(ylim=c(-0.003,0.203))
ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-vt.pdf',width=4.5,height=4.5)
qplot(V4,V8,data=d[e & b & g & h,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique (Avg. Cons.)"))+coord_cartesian(ylim=c(-0.002,0.152))
ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-Vs-Other.pdf',width=4.5,height=4.5)

qplot(V5,V8,data=d[e & b & c,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique-all-pairs"))+coord_cartesian(ylim=c(-0.003,0.203))
ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique-Vs-others.pdf',width=4.5,height=4.5)
 

qplot(V5,V8,data=d[f & b & c,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Average (Cons.)","Average"))+coord_cartesian(ylim=c(-0.003,0.303))
ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique.pdf',width=4.5,height=4.5)

qplot(V4,V8,data=d[f & b & g & h,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Average (Cons.)","Average"))+coord_cartesian(ylim=c(-0.003,0.203))

ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique.pdf',width=4.5,height=4.5)

qplot(V5,V8,data=d[l & i & j & n,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Average (Cons.)","Min (Cons.)","Min (Cons.-sampling)","Average (Cons.-sampling)","Average"))+coord_cartesian(ylim=c(-0.004,0.403))
ggsave('Avian-1ILS-genetreesVSFNrates-Distique-vt.pdf',width=4.5,height=4.5)


qplot(V4,V8,data=d[l & m & i & n,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Average (Cons.)","Min (Cons.)","Min (Cons.-sampling)","Average (Cons.-sampling)","Average"))+coord_cartesian(ylim=c(-0.003,0.303))
ggsave('Avian-1000-genes-ILS-Vs-FNrate-Distique-vt.pdf',width=4.5,height=4.5)

qplot(V5,V8,data=d[j & i & o & n,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique (Avg. Cons.)"))+coord_cartesian(ylim=c(-0.003,0.253))
ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other.pdf',width=4.5,height=4.5)

qplot(V4,V8,data=d[o & m & i & n,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Astrid (njst)","Distique (Avg. Cons.)"))+coord_cartesian(ylim=c(-0.002,0.152))
ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Other.pdf',width=4.5,height=4.5)

qplot(V5,V8,data=d[j & i & p & n,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.1)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Avgerage (Cons.)","Average"))+coord_cartesian(ylim=c(-0.003,0.403))
ggsave('Avian-1ILS-genetreesVSFNrates-Distique.pdf',width=4.5,height=4.5)

qplot(V4,V8,data=d[i & p & n & m,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Avgerage (Cons.)","Average"))+coord_cartesian(ylim=c(-0.002,0.252))
ggsave('Avian-1000-ILS-genetreesVSFNrates-Distique.pdf',width=4.5,height=4.5)




w<-read.csv('11-taxon-res.csv',sep=",",header=FALSE)
levels(w$V4)<-list("M1"="model.10.5400000.0.000000037","M2"="model.10.1800000.0.000000111","M3"="model.10.600000.0.000000333","M4"="model.10.200000.0.000001000")
a<-w$V2 %in% c('astral','RAxML','njst','distique-2-prod')
b<-w$V2 %in% c('distique-2-min','distique-2-prod','distique-cons-min','distique-cons-prod','prod_fm','min_fm')

qplot(V5,V8,data=w[a,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(V4~V6,scales="free_y")+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Astral","Distique (Avg. Cons.)","Astrid","RAxML"))
ggsave('11-taxon-Distique-vs-others.pdf')

qplot(V5,V8,data=w[b,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(V4~V6,scales="free_y")+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Average (Cons.)","Min (Cons.-sample)","Average (Cons.-sample)","Min","Average"))
ggsave('different-methods-tmp-version.pdf')


c<-w$V2 %in% c('distique-2-min','distique-2-prod','prod_fm','min_fm')
qplot(V5,V8,data=w[c,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(V4~V6,scales="free_y")+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("Min (Cons.)","Average (Cons.)","Min","Average"))
ggsave('different-methods.pdf')
