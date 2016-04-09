d<-read.csv('../mammalian-avian-res.csv',sep=",",header=FALSE)
require('ggplot2')
levels(d$V6)<-list("500bp"="500","True gene trees"="true","1000bp"="1000","1500bp"="1500","250bp"="250")
c<-d$V4 %in% c("0.2X")
b<-d$V1 %in% c("mammalian")
g<-d$V5 %in% c(200)
h<-d$V6 %in% c("500bp","True gene trees")
e<-d$V2 %in% c("distique-2-prod","distique-anchoring-log","distique-anchoring-freq","distique-anchoring-freq-2",
               "distique-anchoring-freq-4","distique-anchoring-freq-8",
               "distique-anchoring-log-2","distique-anchoring-log-4",
               "distique-anchoring-log-8")
#e<-d$V2 %in% c("distique-2-prod","distique-2-log","astral")

f<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme")
i<-d$V1 %in% c("avian")
j<-d$V4 %in% c("1X")
l<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme","distique-cons-prod","distique-cons-min")
m<-d$V5 %in% c(1000)
n<-d$V6 %in% c("500bp","True gene trees")
o<-d$V2 %in% c("distique-2-prod","distique-anchoring-log","distique-anchoring-freq","distique-anchoring-freq-2",
               "distique-anchoring-freq-4","distique-anchoring-freq-8",
               "distique-anchoring-log-2","distique-anchoring-log-4",
               "distique-anchoring-log-8")
#o<-d$V2 %in% c("distique-2-prod","distique-2-log","astral")

p<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme")

qplot(V4,V8,data=d[e & b & g & h ,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+
  facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",palette = "Set2",labels=c("DISTIQUE (O(n^4) sum of distance matrices)","DISTIQUE-anchoring (log of sums)","DISTIQUE-anchoring (log of sums, 2 rounds, outlier removed)","DISTIQUE-anchoring (sum of logs)"))


ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-Vs-Other-reg-vs-log2.pdf',width=14.5,height=7.5)

qplot(V5,V8,data=d[e & b & c & d$V2 %in% c("distique-2-prod","distique-anchoring-log","distique-anchoring-log-2","distique-anchoring-log-4","distique-anchoring-log-8"),],
      color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+
  theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+
  scale_color_brewer(name="",palette = "Set2",
                     labels=c("DISTIQUE (O(n^4) sum of distance matrices)","DISTIQUE-anchoring (log of sums)","DISTIQUE-anchoring (log of sums, 2 rounds, outlier removed)","DISTIQUE-anchoring (sum of logs)"))

ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique-Vs-others-reg-vs-log2.pdf',width=14.5,height=7.5)

qplot(V4,V8,data=d[e & b & g & h & d$V2 %in% c("distique-2-prod","distique-anchoring-log","distique-anchoring-log-2","distique-anchoring-log-4","distique-anchoring-log-8") ,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+
  facet_grid(.~V6)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",palette = "Set2",
                     labels=c("DISTIQUE (O(n^4) sum of distance matrices)","DISTIQUE-anchoring (log of sums)","DISTIQUE-anchoring (log of sums, 2 rounds, outlier removed)","DISTIQUE-anchoring (sum of logs)"))

ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-Vs-Other-reg-vs-log-anchors-20.pdf',height=8)

qplot(V5,V8,data=d[e & b & c & d$V2 %in% c("anchores-all","anchores-mrl-20-log",
                                           "anchores-mrl-20-freq","distique-2-prod"),],
      color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+
  theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+
  scale_color_brewer(name="",palette = "Set2",
                     labels=c("anchors (O(n^2) of samples), sum of log","anchors (20 Samp.) MRL log of sum","anchors (20 Samp.) MRL sum of log",
                              "DISTIQUE (O(n^4) sum of distance matrices)"))

ggsave('Mammalian-0.2ILS-genetreesVSFNrates-Distique-Vs-others-reg-vs-log-20.pdf',height=8)




qplot(V5,V8,data=d[j & i & o & n ,],color=V2,group = V2,
      geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) 
    {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+
scale_color_brewer(name="",palette = "Set2",
                   labels=c("DISTIQUE (O(n^4) sum of distance matrices)","DISTIQUE-anchoring (log of sums)","DISTIQUE-anchoring (log of sums, 2 rounds, outlier removed)","DISTIQUE-anchoring (sum of logs)"))

ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other-reg-vs-log2.pdf',width=14.5,height=7.5)


qplot(V4,V8,data=d[o & m & i & n,],
      color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
scale_color_brewer(name="",palette = "Set2",
                   labels=c("DISTIQUE (O(n^4) sum of distance matrices)","DISTIQUE-anchoring (log of sums)","DISTIQUE-anchoring (log of sums, 2 rounds, outlier removed)","DISTIQUE-anchoring (sum of logs)"))

ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Other-reg-vs-log2.pdf',width=14.5,height=7.5)




qplot(V5,V8,data=d[j & i & o & n & d$V2 %in% c("anchores-all",
                                               "anchores-mrl-20-log","anchores-mrl-20-freq",
                                               "distique-2-prod"),],color=V2,group = V2,
      geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) 
  {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+
  scale_color_brewer(name="",palette = "Set2",
                     labels=c("anchors (O(n^2) of samples), sum of log",
                              "anchors (20 Samp.) MRL log of sum","anchors (20 Samp.) MRL sum of log",
                              "DISTIQUE (O(n^4) sum of distance matrices)"))

ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other-reg-vs-log-20.pdf',height=8)


qplot(V4,V8,data=d[o & m & i & n & d$V2 %in% 
                     c("anchores-all",
                       "anchores-mrl-20-log","anchores-mrl-20-freq",
                       "distique-2-prod"),],
      color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",palette = "Set2",
                     labels=c("anchors (O(n^2) of samples), sum of log",
                              "anchors (20 Samp.) MRL log of sum",
                              "anchors (20 Samp.) MRL sum of log",
                              "DISTIQUE (O(n^4) sum of distance matrices)"))

ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Other-reg-vs-log-20.pdf',height=8)






w<-read.csv('../11-taxon-res.csv',sep=",",header=FALSE)
levels(w$V4)<-list("M1"="model.10.5400000.0.000000037","M2"="model.10.1800000.0.000000111","M3"="model.10.600000.0.000000333","M4"="model.10.200000.0.000001000")
#a<-w$V2 %in% c('astral','RAxML','njst_prev','distique-2-prod','mrl','anchores')
a<-w$V2 %in% c('distique-2-prod','distique-2-log','astral')

qplot(V5,V8,data=w[a,],color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(V4~V6,scales="free_y")+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+scale_color_brewer(name="",palette = "Paired",labels=c("ASTRAL","Distique (Avg. Dist. Cons.)","Distique (Avg. Frq. Cons.)"))
ggsave('11-taxon-Distique-vs-others-reg-vs-log.pdf',width=8,height=9)
e<-d$V2 %in% c("distique-2-prod","astrid","astral","mrl-all","anchores-all","anchores-mrl-5-log","anchores-mrl-10-log","anchores-mrl-20-log","anchores-mrl-40-log","anchores-mrl-100-log")



qplot(V7,V8,data=d[e & b & g & h & d$V2 %in% 
                     c( "anchores-all","distique-2-prod","mrl-all","anchores-mrl-10-freq") & 
                     d$V6 =="500bp",],color=V2,group = interaction(V7,V2),geom=c("violin"))+
  facet_wrap(~V4,ncol=1)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+
  ylab("FN rates")+stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
                                fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)


qplot(V7,V8,data=d[e & b & g & h & d$V2 %in% c("anchores-all","mrl-all","anchores-mrl-20-freq") &
                     d$V6 =="500bp",],color=V2,group = interaction(V7,V2),
      geom=c("jitter"))+facet_wrap(~V4,ncol=1)+theme_bw()+theme(legend.position="bottom")+
  xlab("Number of genes")+ylab("FN rates")+
  scale_color_brewer(name="",palette = 2,
                     labels=c("anchors (O(n^2) of samples), sum of log",
                              "anchors (20 Samp.) MRL log of sum",
                              "anchors (O(n^2) Samp.) MRL log of sum"))
ggsave('Mammalian-200-genes-ILS-Vs-FNrate-jitter-500bp.pdf',height=8)


qplot(V7,V8,data=d[o & m & i & n & d$V2 %in% 
                     c("anchores-all","mrl-all",
                       "anchores-mrl-20-freq") & d$V6 =="500bp",],color=V2,group = interaction(V7,V2),
      geom=c("jitter"))+facet_wrap(~V4,ncol=1)+theme_bw()+theme(legend.position="bottom")+
  xlab("Number of genes")+ylab("FN rates")+
  scale_color_brewer(name="",palette = 2,
                     labels=c("anchors (O(n^2) of samples), sum of log",
                              "anchors (20 Samp.) MRL log of sum",
                              "anchors (O(n^2) Samp.) MRL log of sum"))

ggsave('Avian-1000-genetreesVSFNrates-jitter-500bp.pdf',height=8)



