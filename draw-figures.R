require(ggplot2)
require(scales)


################### RF among gene trees
 inc=rbind(data.frame(read.csv("RF.stat",sep="\t",h=F),t="nt-indel (same locus)"),data.frame(read.csv("RF-indel-indel.stat",sep="\t",h=F),t="indel-indel (two random loci)"),data.frame(read.csv("RF-nt-nt.stat",sep="\t",h=F),t="nt-nt (two random loci)"))
qplot(1-V7,data=inc,xlab="Normalized RF distance",geom="density",color=t)+theme_classic()+coord_cartesian(xlim=c(0.34,1))+scale_color_manual(name="",values=c("#e31a1c","#1f78b4","#a6cee3"))+theme(legend.position=c(.25,.86))
ggsave("genetreediscordance-RF-density.png",width=5,height=3.2)
ggsave("genetreediscordance-RF-density.pdf",width=5,height=3.2)
qplot(1-V7,data=inc,xlab="RF",ylab="Number of loci")+theme_bw()+facet_wrap(~t,ncol=1)+coord_cartesian(xlim=c(0.34,1))
ggsave("RF-all.pdf",width=3,height=5)

inc$g="intron"
inc2=rbind(data.frame(read.csv("uce-RF.stat",sep="\t",h=F),t="nt-indel (same locus)"),data.frame(read.csv("uce-RF-indel-indel.stat",sep="\t",h=F),t="indel-indel (two random loci)"),data.frame(read.csv("uce-RF-nt-nt.stat",sep="\t",h=F),t="nt-nt (two random loci)"))
inc2$g="uce"
inc=rbind(inc,inc2)

qplot(1-V7,data=inc,xlab="Normalized RF distance",geom="density",linetype=interaction(t,g,sep=": "), color=interaction(t,g,sep=": "),adjust=1.6)+theme_classic()+coord_cartesian(xlim=c(0.34,1))+scale_color_manual(name="",values=c("#63acce","#1f78b4","#a9ceef","#73c98c","#33a02c","#b2df8a"))+theme(legend.position=c(.25,.78))+scale_linetype_manual(values=c(2,1,1,2,1,1),name="")
ggsave("RF-intron-uce.pdf",width=6,height=3.2)
ggsave("RF-intron-uce.png",width=6,height=3.2)

inc$RF=1-inc$V7
write.table(recast(inc,t~g,measure.var="RF",fun=function(x) paste(format(summary(x)[c(2,3,5)],digits=3),collapse=" ") ),file="rf-random.csv")


################## RI

 qplot(V1,color=V2,data=read.csv("RI.csv",head=F),geom="density")+theme_classic()+scale_color_brewer(palette="Paired",name="")+xlab("Retention Index (RI)")+theme(legend.position=c(.7,.8))
ggsave("RI-genetrees.pdf",width=4,height=2.7)

ggplot(aes(V1,stat(scaled), color=V2),data=read.csv("RI.csv",head=F))+geom_density()+theme_classic()+scale_color_brewer(palette="Paired",name="")+xlab("Retention Index (RI)")+theme(legend.position=c(.88,.8))+ylab("Normalized frequency")
ggsave("RI-genetrees-scaled.png",width=4.2,height=2.9)


################### Compute contraction stats
nw=rbind(read.csv("nwstats.txt",hea=F,sep="\t"))
ggplot(aes(y=(V6)/(V5-1),x=factor(V2,levels=levels(V2)[c(1,4,6,2,3,5)],labels=c("0%","3%","5%","10%","20%","33%")),fill=factor(V1,levels=levels(V1)[c(2,1,4,3)])),data=nw[nw$V2!="no-contract",])+geom_boxplot(outlier.alpha=0.33333,outlier.size=.4)+theme_classic()+scale_fill_brewer(name="",palette="Paired")+xlab("Contraction level (BS \u2264 x contracted)")+ylab("binary nodes (normalized by max possible)")+theme(legend.position="bottom")
ggsave("genetree-resolution-combined.png",width=4.3*1.2, height=3.4*1.2)


################# Bootstrap summary

b=read.csv("bootstrap.stats",header=F,sep=" ")

qplot(V3/100,color=V1,linetype=V2,geom="density",data=b)+theme_classic()+scale_color_brewer(name="",palette="Dark2")+scale_x_continuous(label=percent,name="mean gene tree bootstrap support")+theme(legend.position=c(.8,.75))+scale_linetype(name="")
ggsave("genetree-bootstrap-lt.pdf",width=4.3, height=3.4)
qplot(V3/100,color=factor(paste(V1,V2,sep="-"),levels=c("intron-nt","intron-indel","uce-nt","uce-indel")),geom="density",data=b)+theme_classic()+scale_color_brewer(name="",palette="Paired")+scale_x_continuous(label=percent,name="mean gene tree bootstrap support")+theme(legend.position=c(.8,.75))+scale_linetype(name="")
ggsave("genetree-bootstrap.pdf",width=4.3, height=3.4)
ggsave("genetree-bootstrap.png",width=4.3, height=3.4)

qplot(V3/100,data=b)+theme_classic()+scale_fill_brewer(name="",palette="Dark2")+scale_x_continuous(label=percent,name="mean gene tree bootstrap support")+theme(legend.position="bottom")+facet_grid(V2~V1)
ggsave("genetree-bootstrap-hist.pdf",width=4.3, height=3.4)

qplot(V2,V3/100,geom="boxplot",data=b,fill=V1)+theme_classic()+scale_x_discrete(name="")+scale_y_continuous(label=percent,name="mean gene tree bootstrap support")+theme(legend.position=c(.4,.85))+scale_fill_brewer(name="",palette="Dark2")
ggsave("genetree-bootstrap-boxplot.pdf",width=4.3, height=3.4)

write.table(recast(b,V1~V2,measure.var="V3",fun=function(x) paste(format(summary(x)[c(2,3,5)],digits=3),collapse="/")),file="bootstrap.csv")

#################### Heatmap


 d=read.csv('RF.stats',sep=" ",header=T)
 o=levels(d$V1)[c(3,18:24,2,11:17,5,32:38,4,25:31,1,7:10,6)]; 
 ggplot(aes(x=factor(V1,levels=rev(o)),y=factor(V2,levels=o),fill=V4),data=d)+geom_tile()+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position=c(0.45,-.25),legend.direction="horizontal")+scale_fill_continuous(name="RF",low="white",high="red",breaks=c(0,2,4,6,8,10,12,14))+xlab("")+ylab("")+geom_text(aes(label=V4),color="blue",size=2.8)
 ggsave("heatmap.pdf",width=7.4,height=7.4)
 ggsave("heatmap.png",width=7.4,height=7.4)


#################  Quartet score figure
q=(read.csv("qscore.stat",sep=" ",hea=F))
ggplot(aes(x=factor(V3,levels=levels(q$V3)[c(7,1,4,6,2,3,5)]),color=V1,linetype=V2,group=interaction(V1,V2),y=V4),data=q[q$V1!="combined",])+geom_line()+geom_point()+theme_classic()+theme(legend.position=c(.2,.8),legend.spacing.y=unit(-.3,"cm"))+scale_color_brewer(name="",palette="Dark2",labels=c("intron","UCE"))+xlab("gene tree contraction")+ylab("Normalized ASTRAL quartet score")+scale_linetype(name="")
 ggsave("quartet-score.png",width=4.2,height=3.7)
 ggsave("quartet-score.pdf",width=4.2,height=3.7)



