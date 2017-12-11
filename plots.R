intervention=alpha.nine(1,1,1,1,1,1,1,1,1)
females2017=matrix(0,nrow=80, ncol=8)
for(i in 1:79){
  females2017[1,]<-Prevnum.f(31,2017,intervention)[[1]]
  females2017[i+1,]=Prevnum.f(30+i,2017,intervention)[[i+1]]
}
males2017=matrix(0,nrow=80, ncol=8)
for(i in 1:79){
  males2017[1,]<-Prevnum.m(31,2017,intervention)[[1]]
  males2017[i+1,]=Prevnum.m(30+i,2017,intervention)[[i+1]]
}

total2017=males2017+females2017
plot2017<-data.frame(cbind(age=c(seq(30,109,1)),total2017/1000000)) 
plot2017.b<-data.frame(cbind(age=rep(c(seq(30,109,1)),8),group=c(rep(1,80),rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80),rep(7,80),rep(8,80)),c(total2017[,1:8])/1000000))
write.csv(plot2017,"/Users/n_a_abdallah/Desktop/GSR/prev2017.csv" )
plot2017.bno1<-data.frame(cbind(age=rep(c(seq(30,109,1)),7),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                    rep(5,80), rep(6,80),rep(7,80),rep(8,80)),c(total2017[,2:8])/1000000))
head(plot2017.b) 
#mci pooled
total2017pooled=cbind(total2017[,1:4],total2017[,5]+total2017[,6],total2017[,7:8])
plot2017pool<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                   rep(5,80), rep(6,80), rep(7,80)),c(total2017pooled[,2:7])/1000000))
dpool=data.frame(lt=c(rep("solid",5),rep( "dashed", 1)))
jpeg("/Users/n_a_abdallah/Desktop/GSR/prevbyagenew.jpg",width = 10, height = 6, units = 'in',res = 300)
#quartz()
ggplot(plot2017pool, aes(age,V3+0.01, colour=factor(group),linetype=factor(group))) +
  geom_smooth(method="gam", fill=NA, formula=y~s(x), size=0.92)+
  #scale_y_log10(breaks =c(50,100,200,400,1000),
  #            labels = c(50,100,200,400,1000))+
  xlab("Age")+
  ylab("U.S. Prevalence (millions)")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  scale_color_manual(name="",values=c("#33CC33","#0000CC","#3399CC","#000000","#FF66CC","#CC0033"),
                     labels=c("Amyloidosis","Neurodegeneration","Amyloidosis+Neurodegen.", "MCI due to AD","Early AD", "Late AD"))+
  scale_linetype_manual(values=c(rep("solid",3),rep( "dotted", 1),rep("solid",2)),name="",
                        labels=c("Amyloidosis","Neurodegeneration","Amyloidosis+Neurodegen.", "MCI due to AD","Early AD", "Late AD"))+
  theme(legend.key = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.key.size = unit(1.5, 'lines'),
        legend.text=element_text(size=13),text = element_text(size=12),axis.text=element_text(color = "black", size=12))



dev.off()
quartz()
ggplot(plot2017.bno1, aes(age,V3, colour=group)) +
  geom_path() +
  xlab("Age")+
  ylab("Population in millions")+
  scale_colour_gradientn(colours=rainbow(4))
write.csv(plot2017.b,"/Users/n_a_abdallah/Desktop/GSR/2017.csv")
y=2060
intervention=alpha(1,1,1,1,1,1)
females2060.a=matrix(unlist(Prevnum.f(109,y,intervention)),
                     nrow=80,ncol=6, byrow=T)
females2060.a[1,]<-c(f.census[y-2014+1,36],rep(0,5))
males2060.a=matrix(unlist(Prevnum.m(109,y,intervention)),
                   nrow=80,ncol=6, byrow=T)
males2060.a[1,]<-c(m.census[y-2014+1,36],rep(0,5))
total2060.a=males2060.a+females2060.a

#plot2060<-data.frame(cbind(age=c(seq(30,109,1)),total2060.a,dem.a)) 
plot2060.a<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(1,80),rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80)),c(total2060.a[,1:6])/1000000))
plot2060.ano1<-data.frame(cbind(age=rep(c(seq(30,109,1)),5),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                    rep(5,80), rep(6,80)),c(total2060.a[,2:6])/1000000))
dem.a=total2060.a[,5]+total2060.a[,6]
plot2060.ano1.56<-data.frame(cbind(age=rep(c(seq(30,109,1)),4),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                       rep(5,80)),c(total2060.a[,2:4],dem.a)/1000000))
head(plot2060.ano1.56)  


intervention=alpha(0.5,1,1,1,1,1)
females2060.b=matrix(unlist(Prevnum.f(109,y,intervention)),
                     nrow=80,ncol=6, byrow=T)
females2060.b[1,]<-c(f.census[y-2014+1,36],rep(0,5))
males2060.b=matrix(unlist(Prevnum.m(109,y,intervention)),
                   nrow=80,ncol=6, byrow=T)
males2060.b[1,]<-c(m.census[y-2014+1,36],rep(0,5))
total2060.b=males2060.b+females2060.b
plot2060.b<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(1,80),rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80)),c(total2060.b[,1:6])/1000000))
plot2060.bno1<-data.frame(cbind(age=rep(c(seq(30,109,1)),5),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                    rep(5,80), rep(6,80)),c(total2060.b[,2:6])/1000000))
dem.b=total2060.b[,5]+total2060.b[,6]
plot2060.bno1.56<-data.frame(cbind(age=rep(c(seq(30,109,1)),4),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                       rep(5,80)),c(total2060.b[,2:4],dem.b)/1000000))
head(plot2060.b)  

intervention=alpha(1,1,1,0.5,0.5,1)
females2060.c=matrix(unlist(Prevnum.f(109,y,intervention)),
                     nrow=80,ncol=6, byrow=T)
females2060.c[1,]<-c(f.census[y-2014+1,36],rep(0,5))
males2060.c=matrix(unlist(Prevnum.m(109,y,intervention)),
                   nrow=80,ncol=6, byrow=T)
males2060.c[1,]<-c(m.census[y-2014+1,36],rep(0,5))
total2060.c=males2060.c+females2060.c
plot2060.c<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(1,80),rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80)),c(total2060.c[,1:6])/1000000))
plot2060.cno1<-data.frame(cbind(age=rep(c(seq(30,109,1)),5),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                    rep(5,80), rep(6,80)),c(total2060.c[,2:6])/1000000))
dem.c=total2060.c[,5]+total2060.c[,6]
plot2060.cno1.56<-data.frame(cbind(age=rep(c(seq(30,109,1)),4),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                       rep(5,80)),c(total2060.c[,2:4],dem.c)/1000000))
head(plot2060.c)  

intervention=alpha(0.5,1,1,0.5,0.5,1)
females2060.d=matrix(unlist(Prevnum.f(109,y,intervention)),
                     nrow=80,ncol=6, byrow=T)
females2060.d[1,]<-c(f.census[y-2014+1,36],rep(0,5))
males2060.d=matrix(unlist(Prevnum.m(109,y,intervention)),
                   nrow=80,ncol=6, byrow=T)
males2060.d[1,]<-c(m.census[y-2014+1,36],rep(0,5))
total2060.d=males2060.d+females2060.d
plot2060.d<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(1,80),rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80)),c(total2060.d[,1:6])/1000000))
plot2060.dno1<-data.frame(cbind(age=rep(c(seq(30,109,1)),5),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                    rep(5,80), rep(6,80)),c(total2060.d[,2:6])/1000000))
dem.d=total2060.d[,5]+total2060.d[,6]
plot2060.dno1.56<-data.frame(cbind(age=rep(c(seq(30,109,1)),4),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                       rep(5,80)),c(total2060.d[,2:4],dem.d)/1000000))
head(plot2060.d) 
p1<-ggplot(plot2060.ano1, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  ggtitle("No Intervention")+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#FF66CC","#CC0033"),
                     guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2<-ggplot(plot2060.bno1, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  ggtitle("Alpha12=0.5")+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#FF66CC","#CC0033"),
                     guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3<-ggplot(plot2060.cno1, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  ggtitle("Alpha35=Alpha45=0.5")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#FF66CC","#CC0033"),
                     guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4<-ggplot(plot2060.dno1, aes(age,V3, colour=factor(group))) +
  geom_smooth() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  ggtitle("Alpha12=Alpha35=Alpha45=0.5")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#FF66CC","#CC0033"),
                     guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
quartz()
multiplot(p1, p2, p3, p4, cols=2)
#combining 5+6
p1.56<-ggplot(plot2060.ano1.56, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  ggtitle("No Intervention")+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#CC0033"),
                     labels=c("2","3","4","5+6"),guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2.56<-ggplot(plot2060.bno1.56, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  ggtitle("Alpha12=0.5")+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#CC0033"),
                     labels=c("2","3","4","5+6"),  guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3.56<-ggplot(plot2060.cno1.56, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  ggtitle("Alpha35=Alpha45=0.5")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#CC0033"),
                     labels=c("2","3","4","5+6"),guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4.56<-ggplot(plot2060.dno1.56, aes(age,V3, colour=factor(group))) +
  geom_path() +
  xlab("Age")+
  ylab("Prevalence (million)")+
  ggtitle("Alpha12=Alpha35=Alpha45=0.5")+
  scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
  scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
  scale_color_manual(values=c("#0000CC","#3399CC","#99CCFF","#CC0033"),
                     labels=c("2","3","4","5+6"),guide=guide_legend(title="Stage"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

quartz()
multiplot(p1.56, p2.56, p3.56, p4.56, cols=2)
#prev 5+6 year2018-2060
table<-read.csv("/Users/n_a_abdallah/Desktop/GSR/table.9statesnew.june16.csv")
table1.1<-read.csv("/Users/n_a_abdallah/Desktop/GSR/table1.1.9statesnew.june16.csv")
table2.2<-read.csv("/Users/n_a_abdallah/Desktop/GSR/table2.2.9statesnew.june16.csv")
table3.2<-read.csv("/Users/n_a_abdallah/Desktop/GSR/table3.2.9statesnew.june16.csv")

head(table3.2)
prev.by.year<-data.frame(cbind(group=c(rep(1,44),rep(2,44),rep(3,44),rep(4,44)),
                               year=rep(seq(2017,2060,1),4),dem.census=c((table[3:46,8]+table[3:46,9]),(table1.1[3:46,8]+table1.1[3:46,9]),
                                                                         (table2.2[3:46,8]+table2.2[3:46,9]),(table3.2[3:46,8]+table3.2[3:46,9]))) )        
prev.by.year$group2[prev.by.year$group==1]<-"No Intervention"
prev.by.year$group2[prev.by.year$group==2]<-"I - Lower risk of amyloidosis  "
prev.by.year$group2[prev.by.year$group==3]<-"II- Lower risk of MCI due to AD"
prev.by.year$group2[prev.by.year$group==4]<-"III -Lower risk of  MCI progression to AD"
prev.by.year$group2<-factor(prev.by.year$group2,levels=unique(prev.by.year$group2))
#prev.by.year$group2<-sub("-","-\n",prev.by.year$group2)
jpeg("/Users/n_a_abdallah/Desktop/GSR/prevbyyearnew.jpg",width = 8, height = 4, units = 'in',res = 250)
#quartz()
ggplot(prev.by.year, aes(year,dem.census, group=group2,colour=group2)) +
  geom_line(size=0.7) +
  xlab("Year")+
  ylab("U.S. Prevalence of Alzheimer's Disease (millions)")+
  scale_x_continuous(breaks = c(seq(2015,2060,5)), labels=c("",2020, "",2030,"",2040,"",2050,"",2060))+
  scale_y_continuous(breaks = c(seq(-1,12,1)),labels=c("",0,"",2,"",4,"",6,"",8,"",10,"",12),expand=c(0,0),limits=c(0,12))+
  scale_color_manual(values=c("#330000","#CC0000","green","#0066CC"),
                     guide=guide_legend(title="") )+
  theme(legend.key = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(color = "black",size=12),legend.key.size = unit(1.5, 'lines'),
        legend.text=element_text(size=12))


dev.off()
#incidence plots
incidencemales<-matrix(0, nrow=31, ncol=1)
for ( i in 1:31){
  incidencemales[i,]<-Incidence.m(64+i,y,alpha=intalpha)
}
incidencefemales<-matrix(0, nrow=31, ncol=1)
for ( i in 1:31){
  incidencefemales[i,]<-Incidence.f(64+i,y,alpha=intalpha)
}
write.csv(cbind(incidencemales, incidencefemales),"/Users/n_a_abdallah/Desktop/GSR/incforplots.csv")
incidencemales<-(read.csv("/Users/n_a_abdallah/Desktop/GSR/incforplots.csv"))[,2]
incidencefemales<-(read.csv("/Users/n_a_abdallah/Desktop/GSR/incforplots.csv"))[,3]
av.inc<-100*(incidencefemales+incidencemales)/2
av.inc
#inc.low
incidencemaleslow<-matrix(0, nrow=31, ncol=1)
for ( i in 1:31){
  incidencemaleslow[i,]<-Incidence.m(64+i,y,alpha=intalpha)
}
incidencefemaleslow<-matrix(0, nrow=31, ncol=1)
for ( i in 1:31){
  incidencefemaleslow[i,]<-Incidence.f(64+i,y,alpha=intalpha)
}
write.csv(cbind(incidencemaleslow, incidencefemaleslow),"/Users/n_a_abdallah/Desktop/GSR/incforplotslow.csv")
av.inclow<-100*(incidencefemaleslow+incidencemaleslow)/2
incidencemaleshigh<-matrix(0, nrow=31, ncol=1)
for ( i in 1:31){
  incidencemaleshigh[i,]<-Incidence.m(64+i,y,alpha=intalpha)
}
incidencefemaleshigh<-matrix(0, nrow=31, ncol=1)
for ( i in 1:31){
  incidencefemaleshigh[i,]<-Incidence.f(64+i,y,alpha=intalpha)
}
write.csv(cbind(incidencemaleshigh, incidencefemaleshigh),"/Users/n_a_abdallah/Desktop/GSR/incforplotshigh.csv")
av.inchigh<-100*(incidencefemaleshigh+incidencemaleshigh)/2

incform<-matrix(0,nrow=31,ncol=1)
for(i in 1:31){
  incform[i,]<-0.00117*exp(0.126*(4+i))*100
}
incform
quartz()
jpeg("/Users/n_a_abdallah/Desktop/GSR/incnewlowhigh.jpg",width = 7.5, height = 6, units = 'in',res = 300)
par(las=1)
plot(seq(65,95,1),av.inc,log="y",type="l",xlab="Age",ylab="Incidence (% per year)",bty="L",
     ylim=c(0.05,20))
lines(seq(65,95,1),incform,log="y",col="red")
lines(seq(65,95,1),av.inclow,log="y",col="black",type="l",lty=3)
lines(seq(65,95,1),av.inchigh,log="y",col="black",type="l",lty=3)
legend(85,0.5,c("Multistate Model","Systematic Review"),lwd=c(2.5,2.5),col=c("black","red"),
       bty="n") 
dev.off()
