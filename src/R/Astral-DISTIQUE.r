d<-read.csv('../../Astral-res.csv',sep=",",header=FALSE)
require('ggplot2')
require('scales')
levels(d$V7)<-list("Estimated gene trees"="half","True gene trees"="true")
d$V2<-as.factor(d$V2)
d$V5<-factor(d$V5,levels(d$V5)[c(1,4,2,3,5)])
levels(d$V2)<-list("10M"="10000000","2M"="2000000","500K"="500000")
levels(d$V6)<-list("50"="50gt","200"="200gt","1000"="1000gt")
levels(d$V6)<-as.numeric(as.character(levels(d$V6)))
ggplot(data=d[d$V1 %in% c(200) & d$V7 %in% "Estimated gene trees"  ,])+
  geom_boxplot(aes(x=V6,y=V8,fill=V5),outlier.size=1.2,outlier.colour=rgb(0,0,0,0.7))+
  facet_grid(V3~V2)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  scale_fill_brewer(name="",
                    palette = "BuPu",labels=c("ASTRAL","NJST","Distance-SUM (4 rounds)","Tree-SUM (2 rounds)","CA-ML"))

ggsave("ASTRALII-200taxa.pdf",width=9.5,height=7)

ggplot(data=d[d$V2 %in% c("2M") & d$V3 %in% c(1e-06) & d$V7 %in% "Estimated gene trees"  & d$V1 %in% c(10,50,100,200) ,])+
  geom_boxplot(aes(x=V6,y=V8,fill=V5),outlier.size=1.2,outlier.colour=rgb(0,0,0,0.7))+
  facet_wrap(~V1)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  scale_fill_brewer(name="",
                    palette = "BuPu",labels=c("ASTRAL","NJST","Distance-SUM (4 rounds)","Tree-SUM (2 rounds)","CA-ML"))

ggsave("ASTRALII-different-taxa.pdf",width=6.5,height=7)


ggplot(data=d[d$V1 %in% c(200) & d$V7 %in% "True gene trees"  ,])+
  geom_boxplot(aes(x=V6,y=V8,fill=V5),outlier.size=1.2,outlier.colour=rgb(0,0,0,0.7))+
  facet_grid(V3~V2)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  scale_fill_brewer(name="",
                    palette = "BuPu",labels=c("ASTRAL","Distance-SUM (4 rounds)","Tree-SUM (2 rounds)"))

ggsave("ASTRALII-true-genetrees-200taxa.pdf",width=9.5,height=7)

ggplot(data=d[d$V2 %in% c("2M") & d$V3 %in% c(1e-06) & d$V7 %in% "True gene trees"  & d$V1 %in% c(10,50,100,200) ,])+
  geom_boxplot(aes(x=V6,y=V8,fill=V5),outlier.size=1.2,outlier.colour=rgb(0,0,0,0.7))+
  facet_wrap(~V1)+theme_bw()+theme(legend.position="bottom")+xlab("Number of genes")+ylab("FN rates")+
  scale_fill_brewer(name="",
                    palette = "BuPu",labels=c("ASTRAL","Distance-SUM (4 rounds)","Tree-SUM (2 rounds)"))

ggsave("ASTRALII-true-genetrees-different-taxa.pdf",width=6.5,height=7)