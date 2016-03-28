d<-read.csv('../../mammalian-avian-res.csv',sep=",",header=FALSE)
require('ggplot2')
levels(d$V6)<-list("500bp"="500","True gene trees"="true","1000bp"="1000","1500bp"="1500","250bp"="250")
c<-d$V4 %in% c("0.2X")
b<-d$V1 %in% c("mammalian")
g<-d$V5 %in% c(200)
h<-d$V6 %in% c("500bp","True gene trees")
e1<-d$V2 %in% c("distique-2-prod","distique-Anchoring-SUM-2-rounds","distique-Anchoring-SUM-4-rounds",
               "distique-Anchoring-SUM-1-rounds","distique-Anchoring-SUM-8-rounds","distique-Anchoring-SUM-16-rounds")

e2<-d$V2 %in% c("distique-2-prod","distique-Anchoring-MRL-2-rounds","distique-Anchoring-MRL-4-rounds",
                "distique-Anchoring-MRL-1-rounds","distique-Anchoring-MRL-8-rounds")

e3<-d$V2 %in% c("distique-2-prod","distique-Anchoring-SUM-2-rounds","distique-Anchoring-SUM-4-rounds",
                "distique-Anchoring-MRL-2-rounds","distique-Anchoring-MRL-4-rounds")

e4<-d$V2 %in% c("distique-2-prod","distique-Anchoring-SUM-4-rounds","distique-Anchoring-SUM-8-rounds",
                "distique-Anchoring-MRL-4-rounds","distique-Anchoring-MRL-8-rounds")


#e<-d$V2 %in% c("distique-2-prod","distique-2-log","astral")

f<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme")
i<-d$V1 %in% c("avian")
j<-d$V4 %in% c("1X")
l<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme","distique-cons-prod","distique-cons-min")
m<-d$V5 %in% c(1000)
n<-d$V6 %in% c("500bp","True gene trees")
o1<-d$V2 %in% c("distique-2-prod","distique-Anchoring-SUM-2-rounds","distique-Anchoring-SUM-4-rounds",
                "distique-Anchoring-SUM-1-rounds","distique-Anchoring-SUM-8-rounds","distique-Anchoring-SUM-16-rounds")

o2<-d$V2 %in% c("distique-2-prod","distique-Anchoring-MRL-2-rounds","distique-Anchoring-MRL-4-rounds",
                "distique-Anchoring-MRL-1-rounds","distique-Anchoring-MRL-8-rounds")

o3<-d$V2 %in% c("distique-2-prod","distique-Anchoring-SUM-2-rounds","distique-Anchoring-SUM-4-rounds",
                "distique-Anchoring-MRL-2-rounds","distique-Anchoring-MRL-4-rounds")

o4<-d$V2 %in% c("distique-2-prod","distique-Anchoring-SUM-4-rounds","distique-Anchoring-SUM-8-rounds",
                "distique-Anchoring-MRL-4-rounds","distique-Anchoring-MRL-8-rounds")

#o<-d$V2 %in% c("distique-2-prod","distique-2-log","astral")

p<-d$V2 %in% c("distique-2-prod","distique-2-min","prod_fm-fastme")

qplot(V4,V8,data=d[e1 & b & g & h ,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+
  facet_grid(.~V6)+theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Distance-SUM (1 round)","Distance-SUM (16 rounds)","Distance-SUM (2 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 round)"))


ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-vs-Anchoring-Distance-Sum_different-rounds.pdf',width=9.5,height=4)

qplot(V5,V8,data=d[e1 & b & c ,],
      color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+
  theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Distance-SUM (1 round)","Distance-SUM (16 rounds)","Distance-SUM (2 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))

ggsave('Mammalian-0.2ILS-genetrees-VS-FNrates-Distique-vs-Anchoring-Distance-SUM_different-rounds.pdf',width=9.5,height=4)





qplot(V4,V8,data=d[e2 & b & g & h ,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+
  facet_grid(.~V6)+theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (1 round)","Tree-SUM (2 rounds)",
                                               "Tree-SUM (4 rounds)","Tree-SUM (8 rounds)"))

ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-vs-Anchoring-Tree-Sum_different-rounds.pdf',width=9.5,height=4)

qplot(V5,V8,data=d[e2 & b & c ,],
      color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+
  theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=10)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (1 round)","Tree-SUM (2 rounds)",
                                               "Tree-SUM (4 rounds)","Tree-SUM (8 rounds)"))

ggsave('Mammalian-0.2ILS-genetrees-VS-FNrates-Distique-vs-Anchoring-Tree-Sum_different-rounds.pdf',width=9.5,height=4)





qplot(V4,V8,data=d[e3 & b & g & h ,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+
  facet_grid(.~V6)+theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (2 rounds)","Tree-SUM (4 rounds)",
                                               "Distance-SUM (2 rounds)","Distance-SUM (4 rounds)"))


ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-vs-Anchoring-Tree-Sum-VS-Distance-Sum-different-rounds.pdf',width=9.5,height=4)

qplot(V5,V8,data=d[e3 & b & c ,],
      color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+
  theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (2 rounds)","Tree-SUM (4 rounds)",
                                               "Distance-SUM (2 rounds)","Distance-SUM (4 rounds)"))

ggsave('Mammalian-0.2ILS-genetrees-VS-FNrates-Distique-vs-Anchoring-Tree-Sum-VS-Distance-Sum-different-rounds.pdf',width=9.5,height=4)


qplot(V4,V8,data=d[e4 & b & g & h ,],color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+
  facet_grid(.~V6)+theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (4 rounds)","Tree-SUM (8 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))

ggsave('Mammalian-200-genes-ILS-Vs-FNrate-Distique-vs-Anchoring-Tree-Sum-VS-Distance-Sum-different-rounds-4-vs-8.pdf',width=9.5,height=4)




qplot(V5,V8,data=d[e4 & b & c ,],
      color=V2,group = V2,geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+
  theme_bw()+theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (4 rounds)","Tree-SUM (8 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))

ggsave('Mammalian-0.2ILS-genetrees-VS-FNrates-Distique-vs-Anchoring-Tree-Sum-VS-Distance-Sum-different-rounds-2-vs-4.pdf',width=9.5,height=4)





qplot(V5,V8,data=d[j & i & o1 & n ,],color=V2,group = V2,
      geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) 
  {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.5)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Distance-SUM (1 round)","Distance-SUM (16 rounds)","Distance-SUM (2 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))

ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other-reg-vs-Anchoring-Distance-Sum-different_rounds.pdf',width=9.5,height=4)


qplot(V4,V8,data=d[o1 & m & i & n,],
      color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Distance-SUM (1 round)","Distance-SUM (16 rounds)","Distance-SUM (2 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))
ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Anchoring-Distance-Sum-different_rounds.pdf',width=9.5,height=4)





qplot(V5,V8,data=d[j & i & o2 & n ,],color=V2,group = V2,
      geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) 
  {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.5)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (1 round))","Tree-SUM (2 rounds)",
                                               "Tree-SUM (4 rounds)","Tree-SUM (8 rounds)"))
ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other-reg-vs-Anchoring-Tree-SUM-different_rounds.pdf',width=9.5,height=4)


  qplot(V4,V8,data=d[o2 & m & i & n,],
        color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (1 round))","Tree-SUM (2 rounds)",
                                               "Tree-SUM (4 rounds)","Tree-SUM (8 rounds)"))
ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Anchoring-Tree-Sum-different_rounds.pdf',width=9.5,height=4)






qplot(V5,V8,data=d[j & i & o3 & n ,],color=V2,group = V2,
      geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) 
  {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.5)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (2 rounds))","Tree-SUM (4 rounds)",
                                               "Distance-SUM (2 rounds)","Distance-SUM (4 rounds)"))
ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other-reg-vs-Anchoring-Tree-Sum-Vs-Distance-sum-different_rounds-2-4.pdf',width=9.5,height=4)


qplot(V4,V8,data=d[o3 & m & i & n,],
      color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (2 rounds))","Tree-SUM (4 rounds)",
                                               "Distance-SUM (2 rounds)","Distance-SUM (4 rounds)"))
ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Anchoring-Tree-Sum-Vs-Distance-Sum-different_rounds-2-4.pdf',width=9.5,height=4)



qplot(V5,V8,data=d[j & i & o4 & n ,],color=V2,group = V2,
      geom=c("point","line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) 
  {mean(x)-sd(x)/sqrt(length(x))},fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.5)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (4 rounds))","Tree-SUM (8 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))
ggsave('Avian-1ILS-genes-ILSVsFNrate-Distique-Vs-Other-reg-vs-Anchoring-Tree-Sum-Vs-Distance-Sum-different_rounds-4_VS-8.pdf',width=9.5,height=4)


qplot(V4,V8,data=d[o4 & m & i & n,],
      color=V2,group = V2,geom=c("line"),fun.y=mean,stat='summary')+facet_grid(.~V6)+theme_bw()+
  theme(legend.position="right")+xlab("Number of genes")+ylab("FN rates")+
  stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.05)+
  scale_color_brewer(name="",
                     palette = "Set2",labels=c("DISTIQUE (All Pairs)",
                                               "Tree-SUM (4 rounds))","Tree-SUM (8 rounds)",
                                               "Distance-SUM (4 rounds)","Distance-SUM (8 rounds)"))
ggsave('Avian-1000-genetreesVSFNrates-Distique-Vs-Anchoring-Tree-Sum-Vs-Distance-Sum-different_rounds-4_Vs-8.pdf',width=9.5,height=4)





