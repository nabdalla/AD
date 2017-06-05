
rates<-read.csv("/Users/n_a_abdallah/Desktop/GSR/mortality rates all.csv")
#females
dr.f<-function(a,y){
  rates.sub<-subset(rates, rates[,1]==y)
  d.f<-rates.sub[a+1,3]
  return(d.f)
}
#males
dr.m<-function(a,y){
  rates.sub<-subset(rates, rates[,1]==y)
  d.f<-rates.sub[a+1,4]
  return(d.f)
}
#posible interventions
alpha.nine<-function(one.two,  two.four, four.five, one.three, three.four,three.six, 
                five.seven, six.seven, seven.eight){
  alpha<-matrix(1,nrow=8,ncol=9)
  alpha[1,2]<-one.two
  alpha[2,4]<-two.four
  alpha[4,5]<-four.five
  alpha[1,3]<-one.three
  alpha[3,4]<-three.four
  alpha[3,6]<-three.six
  alpha[5,7]<-five.seven
  alpha[6,7]<-six.seven
  alpha[7,8]<-seven.eight
  return(alpha)
}
intalpha=alpha.nine(1,1,1,1,1,1,1,1,1)
#constants of the probabilities of transition from state i to j
k0<-matrix(0,6,7)
k1<-matrix(0,6,7)
k0[1,2]<-0.000149
k1[1,2]<-0.086
k0[1,3]<-7.531246e-06
k1[1,3]<-0.117866
k0[3,4]<-0.00012646
k1[3,4]<-0.081938
k0[2,4]<-0.001109085
k1[2,4]<-0.0616514
k0[3,6]<-0.00327305
k1[3,6]<-0.01983114
k0[6,7]<-0.09
k1[6,7]<-0
k0[4,5]<-2.595831e-05
k1[4,5]<-0.09623778
k0[5,7]<-0.3
k1[5,7]<-0
#The one step transition matrix at age a in year y  
#females
#define intervention year
int.year=2017
pr.y.alpha<-function(a,y,alpha){
  if(y<int.year){
    alpha=matrix(1,nrow=8,ncol=9)
  }else if(y>=int.year){
    alpha=alpha
  }
  if(y>2014){
    y=2014
  } else if (y<1933){
    y<-1933
  }
  else {y=y}
  if(a>109){
    a=110
  } else {
    a=a
  }
  p<-matrix(0,8,9)
  p[7,8]<-alpha[7,8]*0.167*(1-dr.f(a,y))
  p[7,9]<-alpha[7,9]*dr.f(a,y)
  p[8,9]<-alpha[8,9]*(dr.f(a,y)+0.078)
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.f(a,y))}
        else if(65 <= a & a <= 75){
          p[1,2]<-alpha[1,2]*0.00105*exp(0.05596*a)*(1-dr.f(a,y))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.f(a,y))
        }
        else if(a>75){
          p[1,2]<-alpha[1,2]*0.07*(1-dr.f(a,y))
          p[3,6]<-alpha[3,6]*k0[3,6]*exp(k1[3,6]*75)*(1-dr.f(a,y))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.f(a,y))
        }
        p[i,9]<-alpha[i,9]*dr.f(a,y)}
    }
  }
  for(i in 1:8){
    if(a<100){
      p[i,i]<-1-sum(p[i, 1:9]) }
    else if (a>99){
      p[1,1]<-0
      p[i,i]<-1-sum(p[i, 1:9]) }
  }
  
  return(p)}
pr.y.alpha(65,2017,intalpha)
#males
pr.y.males.alpha<-function(a,y,alpha){
  if(y<int.year){
    alpha=matrix(1,nrow=8,ncol=9)
  }else if(y>=int.year){
    alpha=alpha
  }
  if(y>2014){
    y=2014
  } else if (y<1933){
    y<-1933
  }
  else {y=y}
  if(a>109){
    a=110
  } else {
    a=a
  }
  p<-matrix(0,8,9)
  p[7,8]<-alpha[7,8]*0.167*(1-dr.m(a,y))
  p[7,9]<-alpha[7,9]*dr.m(a,y)
  p[8,9]<-alpha[8,9]*(dr.m(a,y)+0.078)
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.m(a,y))}
        else if(65 <= a & a <= 75){
          p[1,2]<-alpha[1,2]*0.00105*exp(0.05596*a)*(1-dr.m(a,y))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.m(a,y))
        }
        else if(a>75){
          p[1,2]<-alpha[1,2]*0.07*(1-dr.m(a,y))
          p[3,6]<-alpha[3,6]*k0[3,6]*exp(k1[3,6]*75)*(1-dr.m(a,y))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.m(a,y))
        }
        p[i,9]<-alpha[i,9]*dr.m(a,y)}
    }
  }
  for(i in 1:8){
    if(a<100){
      p[i,i]<-1-sum(p[i, 1:9]) }
    else if (a>99){
      p[1,1]<-0
      p[i,i]<-1-sum(p[i, 1:9]) }
  }
  return(p)}

#The absolute risk of dementia is the probability given you are in state I at age a in 
#year y, what is the probability of reaching state j in t years
#females
prob.adj<-prob.adj.m<-NULL
pr.rec.adj<-function(a,y){
  prob.adj<-rbind(pr.y.alpha(a,y,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))[1:6,],
                  c(rep(0,8),1),pr.y.alpha(a,y,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))[8,],c(rep(0,8),1))
  return(prob.adj)
}
prod<-sum<-vector("list")
AR.f<-function(a,t,y){
  for(i in 1:(t-1)){
    prod[[1]]<-pr.rec.adj(a,y)
    prod[[i+1]]<-(prod[[i]])%*%pr.rec.adj(a+i, y+i)
  }
  for(k in 1:(t-1)){
    sum[[1]]<-prod[[1]]
    sum[[k+1]]<-sum[[k]]+prod[[k+1]]
  }
  
  return(sum[[t]])
}
y=2014
trans.rec5<-matrix(0,6,6 )
trans.rec10<-matrix(0,6,6 )
trans.rec15<-matrix(0,6,6 )

for(i in 1:6){
  trans.rec5[i,]<-AR.f((65+((i-1)*5)),5,y)[1:6,7]
  trans.rec10[i,]<-AR.f((65+((i-1)*5)),10,y)[1:6,7]
  trans.rec15[i,]<-AR.f((65+((i-1)*5)),15,y)[1:6,7]
}

#males
pr.rec.adj.males<-function(a,y){
  prob.adj.m<-rbind(pr.y.males.alpha(a,y,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))[1:6,],c(rep(0,8),1),
                    pr.y.males.alpha(a,y,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))[8,],c(rep(0,8),1))
  return(prob.adj.m)
}
prod.m<-sum.m<-vector("list")
AR.m<-function(a,t,y){
  for(i in 1:(t-1)){
    prod.m[[1]]<-pr.rec.adj.males(a,y)
    prod.m[[i+1]]<-(prod.m[[i]])%*%pr.rec.adj.males(a+i, y+i)
  }
  for(k in 1:(t-1)){
    sum.m[[1]]<-prod.m[[1]]
    sum.m[[k+1]]<-sum.m[[k]]+prod.m[[k+1]]
  }
  
  return(sum.m[[t]])
}
y=2014
trans.rec5.males<-matrix(0,6,6 )
trans.rec10.males<-matrix(0,6,6 )
trans.rec15.males<-matrix(0,6,6 )

for(i in 1:6){
  trans.rec5.males[i,]<-AR.m((65+((i-1)*5)),5,y)[1:6,7]
  trans.rec10.males[i,]<-AR.m((65+((i-1)*5)),10,y)[1:6,7]
  trans.rec15.males[i,]<-AR.m((65+((i-1)*5)),15,y)[1:6,7]
}


#The SurvAD survival for a person who has Alzheimer's at age a in year y
prob.adj.exp<-NULL
pr.rec.adj.exp<-function(a,y){
  prob.adj.exp<-rbind(cbind(pr.y.alpha(a,y,alpha=alpha.nine(1,1,1,1,1,1,1,1,1)),rep(0,8)),c(rep(0,9),1),c(rep(0,9),1))
  return(prob.adj.exp)
}
SurvAD.f<-function(a,n=(108-a),y){
  exp.surv<-NULL
  for(i in 1:(n-1)){
    prod[[1]]<-pr.rec.adj.exp(a,y)
    prod[[i+1]]<-(prod[[i]])%*%pr.rec.adj.exp(a+(i+1)-1,y)
  }
  for(j in 1:(107-a)){
    exp.surv[1]<-prod[[1]][7,9]
    exp.surv[j+1]<-(exp.surv[j])+((j+1)*prod[[j+1]][7,9])
  }
  return(exp.surv[n])
}
exp.survival.f<-matrix(0,nrow=8, ncol=1)
for(i in 1:8){
  a[i]=(60+(i-1)*5)
  exp.survival.f[i,1]<-SurvAD.f(a[i],(108-a[i]),2014)
}
#males
prob.adj.exp.m<-NULL
pr.rec.adj.exp.males<-function(a,y){
  prob.adj.exp.m<-rbind(cbind(pr.y.males.alpha(a,y,alpha=alpha.nine(1,1,1,1,1,1,1,1,1)),rep(0,8)),c(rep(0,9),1),c(rep(0,9),1))
  return(prob.adj.exp.m)
}

SurvAD.m<-function(a,n=(108-a),y){
  exp.surv.m<-NULL
  for(i in 1:(n-1)){
    prod.m[[1]]<-pr.rec.adj.exp.males(a,y)
    prod.m[[i+1]]<-(prod.m[[i]])%*%pr.rec.adj.exp.males(a+(i+1)-1,y)
  }
  for(j in 1:(107-a)){
    exp.surv.m[1]<-prod.m[[1]][7,9]
    exp.surv.m[j+1]<-(exp.surv.m[j])+((j+1)*prod.m[[j+1]][7,9])
  }
  return(exp.surv.m[n])
}
exp.survival.m<-matrix(0,nrow=8, ncol=1)
for(i in 1:8){
  a[i]=(60+(i-1)*5)
  exp.survival.m[i,1]<-SurvAD.m(a[i],(108-a[i]),2014)
}

#The prevalence rate in state I for persons at age a in year y is the proportion of persons alive at age a in calendar year y who are in state I 
#females
prob<-NULL

TP.f<-function(a,y, alpha){
  prob<-rbind(pr.y.alpha(a,y,alpha),c(rep(0,8),1))
  return(prob)
}
Prevrate.f<-function(a,y, alpha){
  n=a-30  
  y2<-y-n-1
  prevalence<-NULL
  #The probability of being in state I at age a in year y given that you started in state 1 at age (30) is 
  #\phi_i(a,y)
  for(i in 1:(n)){
    prod[[1]]<-TP.f(30,y2,alpha)
    prod[[i+1]]<-(prod[[i]])%*%TP.f(30+i, y2+i,alpha)}
  for(i in 1:8){
    prevalence[i]<-prod[[n+1]][1,i]/sum(prod[[n+1]][1,1:8])
  }
  return(c(prevalence))
}
1-dr.f(89,2014)
#males
TP.m<-function(a,y, alpha){
  prob<-rbind(pr.y.males.alpha(a,y,alpha),c(rep(0,8),1))
  return(prob)
}
Prevrate.m<-function(a,y, alpha){
  n=a-30
  y2<-y-n-1
  prevalence<-NULL
  for(i in 1:(n)){
    prod[[1]]<-TP.m(30,y2,alpha)
    prod[[i+1]]<-(prod[[i]])%*%TP.m(30+i, y2+i,alpha)
  }
  for(i in 1:8){
    prevalence[i]<-prod[[n+1]][1,i]/sum(prod[[n+1]][1,1:8])
  }
  return(c(prevalence))
}

#The incidence rate of dementia at age a in year y:  I (a,y). This is a conditional probability.
#females
incidence.females<-NULL
Incidence.f<-function(a,y,alpha){
  incidence.females<-((Prevrate.f(a-1,y, alpha)[6]*TP.f(a-1,y-1,alpha)[6,7])+
                        (Prevrate.f( a-1,y, alpha)[5]*TP.f(a-1,y-1,alpha)[5,7]))/
    sum(Prevrate.f(a-1,y,alpha)[1:6])
  return(incidence.females)
}
incidencefemales<-matrix(0, nrow=9, ncol=1)
for ( i in 1:9){
  incidencefemales[i,]<-Incidence.f(55+(5*i),y,alpha=intalpha)
}

#without death adjustment in transition rates
incidence.females<-NULL
Incidence.f<-function(a,y,alpha){
  incidence.females<-(((Prevrate.f(a-1,y, alpha)[6]*TP.f(a-1,y-1,alpha)[6,7])+
                        (Prevrate.f( a-1,y, alpha)[5]*TP.f(a-1,y-1,alpha)[5,7]))/
                        (1-dr.f(a-1,2014)))/
    sum(Prevrate.f(a-1,y,alpha)[1:6])
  return(incidence.females)
}
incidencefemales<-matrix(0, nrow=9, ncol=1)
for ( i in 1:9){
  incidencefemales[i,]<-Incidence.f(55+(5*i),y,alpha=intalpha)
}
#males
incidence.males<-NULL
Incidence.m<-function(a,y,alpha){
  incidence.males<-((Prevrate.m(a-1,y,alpha)[6]*TP.m(a-1,y-1,alpha)[6,7])+
                      (Prevrate.m( a-1,y,alpha)[5]*TP.m(a-1,y-1,alpha)[5,7]))/
    sum(Prevrate.m(a-1,y,alpha)[1:6])
  return(incidence.males)
}

incidencemales<-matrix(0, nrow=9, ncol=1)
for ( i in 1:9){
  incidencemales[i,]<-Incidence.m(55+(5*i),y,alpha=intalpha)
}

#no death adjustment
incidence.males<-NULL
Incidence.m<-function(a,y,alpha){
  incidence.males<-(((Prevrate.m(a-1,y,alpha)[6]*TP.m(a-1,y-1,alpha)[6,7])+
                      (Prevrate.m( a-1,y,alpha)[5]*TP.m(a-1,y-1,alpha)[5,7]))/
    (1-dr.m(a-1,2014)))/
    sum(Prevrate.m(a-1,y,alpha)[1:6])
  return(incidence.males)
}

incidencemales<-matrix(0, nrow=9, ncol=1)
for ( i in 1:9){
  incidencemales[i,]<-Incidence.m(55+(5*i),y,alpha=intalpha)
}
write.csv(cbind(incidencemales, incidencefemales), "/Users/n_a_abdallah/Desktop/GSR/inc9statesnodeath.April24.csv")

#projections

#Projections of numbers of people in the United States in year y state i are given by Mi (y)

#The n-step transition matrix. 
#The elements represent the probability that given you are in state I at age a in calendar year y, that you are in state j n steps (years later)
#function used to calculate prevalence matrix R for the first year to be
#used in recursive function males and females
TPN.m<-function(a,n, y,alpha){
  prod<-vector("list")
  for(i in 1:(n)){
    prod[[1]]<-TP.m(a,y-n-1,alpha)
    prod[[i+1]]<-(prod[[i]])%*%TP.m(a+i, y-n+i-1,alpha)
  }
  return(prod)
}
TPN.f<-function(a,n,y,alpha){
  prod<-vector("list")
  for(i in 1:(n)){
    prod[[1]]<-TP.f(a,(y-n-1),alpha)
    prod[[i+1]]<-(prod[[i]])%*%TP.f(a+i, (y-n+i-1),alpha)
  }
  return(prod)
}

#The number of people in the U.S at age a, in year y
Prevnum.f<-function(a, y, alpha){
  prob<-pf<-num.females<-prevalencef<-vector("list") 
  prob<-TPN.f(a=30,n=a-30,y-1,alpha=alpha)
  for(k in 1:(a-30)){
    pf[[1]]<-((TP.f(30,y,alpha)))
    pf[[k+1]]<-(prob[[k]]%*%(TP.f(30+k,(y-1),alpha)))
    prevalencef[[1]]<-pf[[1]][1,1:8]/sum(pf[[(1)]][1,1:8])
    prevalencef[[k+1]]<-pf[[k+1]][1,1:8]/sum(pf[[(k+1)]][1,1:8])
    num.females[[1]]<-prevalencef[[1]]*f.census[(y-2014+1),(6+(30))]
    num.females[[k+1]]<-prevalencef[[k]]*f.census[(y-2014+1),(6+(30+k))]}
  return(num.females)
}
Prevnum.m<-function(a, y, alpha){
  prob<-pm<-num.males<-prevalencem<-vector("list") 
  prob<-TPN.m(a=30,n=a-30,y-1,alpha=alpha)
  for(k in 1:(a-30)){
    pm[[1]]<-((TP.m(30,y,alpha)))
    pm[[k+1]]<-(prob[[k]]%*%(TP.m(30+k,(y-1),alpha)))
    prevalencem[[1]]<-pm[[1]][1,1:8]/sum(pm[[(1)]][1,1:8])
    prevalencem[[k+1]]<-pm[[k+1]][1,1:8]/sum(pm[[(k+1)]][1,1:8])
    num.males[[1]]<-prevalencem[[1]]*m.census[(y-2014+1),(6+(30))]
    num.males[[k+1]]<-prevalencem[[k]]*m.census[(y-2014+1),(6+(30+k))]}
  return(num.males)
}
#Projections of numbers of people in the United States in year y state i are given by Mi (y)
Prevnum.Total<-function(prob1,prob2, y, alpha){
  a=109
  p<-pf<-num.males<-num.females<-num.tot<-prevalencef<- prevalence<-inc.males<-sum.inc.males<-
    mci4<-inc.females<-sum.inc.females<-inc.tot<-inc.mci4<-vector("list")   
  for(k in 1:(a-30)){
    p[[1]]<-((TP.m(30,y-1,alpha)))
    prevalence[[1]]<-p[[1]][1,1:8]/sum(p[[(1)]][1,1:8])
    #num.males[[1]]<-prevalence[[1]]*m.census[(i+1),(6+(30))]
    num.males[[1]]<-prevalence[[1]]*m.census[(i+1),(6+(30))]
    inc.males[[1]]<-0
    p[[k+1]]<-(prob1[[k]]%*%(TP.m(30+k,(y-1),alpha)))
    prevalence[[k+1]]<-p[[k+1]][1,1:8]/sum(p[[(k+1)]][1,1:8])
    #The number of persons in the United States by disease state at age a in year y, gender g
    num.males[[k+1]]<-prevalence[[k]]*m.census[(i+1),(6+(30+k))]
    inc.males[[k+1]]<-(num.males[[k+1]][6]*trans6to7.m[i,k+1])+
      (num.males[[k+1]][5]*trans5to7.m[i,k+1])
    pf[[1]]<-((TP.f(30,y-1,alpha)))
    prevalencef[[1]]<-pf[[1]][1,1:8]/sum(pf[[(1)]][1,1:8])
    #num.females[[1]]<-prevalencef[[1]]*f.census[(i+1),(6+(30))]
    num.females[[1]]<-prevalencef[[1]]*f.census[(i+1),(6+(30))]
    inc.females[[1]]<-0
    pf[[k+1]]<-(prob2[[k]]%*%(TP.f(30+k,(y-1),alpha)))
    prevalencef[[k+1]]<-pf[[k+1]][1,1:8]/sum(pf[[(k+1)]][1,1:8])
    num.females[[k+1]]<-prevalencef[[k]]*f.census[(i+1),(6+(30+k))]
    #The incidence of  Alzheimer's disease at age a in year y (numbers of people, in millions)   
    inc.females[[k+1]]<-(num.females[[k+1]][6]*trans6to7.f[i,k+1])+
      (num.females[[k+1]][5]*trans5to7.f[i,k+1])
    mci4[[1]]<-0
    mci4[[k+1]]<-((num.females[[k+1]][5]*trans5to7.f[i,k+1])+
      (num.males[[k+1]][5]*trans5to7.m[i,k+1]))
  }
  sum.females<-Reduce('+', num.females)
  sum.males<-Reduce('+', num.males) 
  sum.tot<-sum.males+sum.females
  sum.inc.males<-Reduce('+', inc.males)
  sum.inc.females<-Reduce('+', inc.females)
  inc.tot<-sum.inc.males+sum.inc.females
  inc.mci4<-Reduce('+', mci4)
  return(c(p,pf,sum.tot,inc.tot, inc.mci4))
}  
#define starting year and ending year for projections
year=2014
year2=2060
#define the intervention
intalpha<-alpha.nine(1,1,1,1,1,1,1,1,1)

trans6to7.m<-matrix(0,year2-year,(109-29))
trans5to7.m<-matrix(0,year2-year,(109-29))
trans6to7.f<-matrix(0,year2-year,(109-29))
trans5to7.f<-matrix(0,year2-year,(109-29))
prob1<-prob2<-fit<-sum.year<-inc.year<-vector("list")
for(l in 1:(109-29)){
  trans6to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha)[6,7]
  trans5to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha)[5,7]
  trans6to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha)[6,7]
  trans5to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha)[5,7]
}
for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum<-unlist(sum.year)
table<-matrix((sum)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table, "/Users/n_a_abdallah/Desktop/GSR/table.9states.may2.csv")

#Primary intervention
#Amyloid
intalpha1.1<-alpha.nine(0.25,1,1,1,0.25,1,1,1,1)
for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha1.1)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum1.1<-unlist(sum.year)
table1.1<-matrix((sum1.1)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table1.1, "/Users/n_a_abdallah/Desktop/GSR/table1.1.9states.April24.csv")

intalpha1.2<-alpha.nine(0.5,1,1,1,0.5,1,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha1.2)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum1.2<-unlist(sum.year)
table1.2<-matrix((sum1.2)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table1.2, "/Users/n_a_abdallah/Desktop/GSR/table1.2.9states.April24.csv")

intalpha1.3<-alpha.nine(0.75,1,1,1,0.75,1,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha1.3)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum1.3<-unlist(sum.year)
table1.3<-matrix((sum1.3)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table1.3, "/Users/n_a_abdallah/Desktop/GSR/table1.3.9states.April24.csv")

intalpha1.4<-alpha.nine(0.9,1,1,1,0.9,1,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha1.4)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum1.4<-unlist(sum.year)
table1.4<-matrix((sum1.4)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table1.4, "/Users/n_a_abdallah/Desktop/GSR/table1.4.9states.April24.csv")

#secondary
#Target progression to MCI
intalpha2.1<-alpha.nine(1,1,0.25,1,1,0.25,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha2.1)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum2.1<-unlist(sum.year)
table2.1<-matrix((sum2.1)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table2.1, "/Users/n_a_abdallah/Desktop/GSR/table2.1.9states.April24.csv")

intalpha2.2<-alpha.nine(1,1,0.5,1,1,0.5,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha2.2)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum2.2<-unlist(sum.year)
table2.2<-matrix((sum2.2)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table2.2, "/Users/n_a_abdallah/Desktop/GSR/table2.2.9states.April24.csv")

intalpha2.3<-alpha.nine(1,1,0.75,1,1,0.75,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha2.3)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum2.3<-unlist(sum.year)
table2.3<-matrix((sum2.3)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table2.3, "/Users/n_a_abdallah/Desktop/GSR/table2.3.9states.April24.csv")

intalpha2.4<-alpha.nine(1,1,0.9,1,1,0.9,1,1,1)

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha2.4)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum2.4<-unlist(sum.year)
table2.4<-matrix((sum2.4)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table2.4, "/Users/n_a_abdallah/Desktop/GSR/table2.4.9states.April24.csv")

#secondary prevention
#Target progression MCI to AD 5-7 and 6-7
intalpha3.1<-alpha.nine(1,1,1,1,1,1,0.5,0.5,1)

trans6to7.m<-matrix(0,year2-year,(109-29))
trans5to7.m<-matrix(0,year2-year,(109-29))
trans6to7.f<-matrix(0,year2-year,(109-29))
trans5to7.f<-matrix(0,year2-year,(109-29))
prob1<-prob2<-fit<-sum.year<-inc.year<-vector("list")
for(l in 1:(109-29)){
  trans6to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.1)[6,7]
  trans5to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.1)[5,7]
  trans6to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.1)[6,7]
  trans5to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.1)[5,7]
}

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha3.1)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum3.1<-unlist(sum.year)
table3.1<-matrix((sum3.1)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table3.1, "/Users/n_a_abdallah/Desktop/GSR/table3.1.9states.April24.csv")

intalpha3.2<-alpha.nine(1,1,1,1,1,1,0.75,0.75,1)
trans6to7.m<-matrix(0,year2-year,(109-29))
trans5to7.m<-matrix(0,year2-year,(109-29))
trans6to7.f<-matrix(0,year2-year,(109-29))
trans5to7.f<-matrix(0,year2-year,(109-29))
prob1<-prob2<-fit<-sum.year<-inc.year<-vector("list")
for(l in 1:(109-29)){
  trans6to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.2)[6,7]
  trans5to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.2)[5,7]
  trans6to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.2)[6,7]
  trans5to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.2)[5,7]
}

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha3.2)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum3.2<-unlist(sum.year)
table3.2<-matrix((sum3.2)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table3.2, "/Users/n_a_abdallah/Desktop/GSR/table3.2.9states.April24.csv")

intalpha3.3<-alpha.nine(1,1,1,1,1,1,0.90,0.90,1)
trans6to7.m<-matrix(0,year2-year,(109-29))
trans5to7.m<-matrix(0,year2-year,(109-29))
trans6to7.f<-matrix(0,year2-year,(109-29))
trans5to7.f<-matrix(0,year2-year,(109-29))
prob1<-prob2<-fit<-sum.year<-inc.year<-vector("list")
for(l in 1:(109-29)){
  trans6to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.3)[6,7]
  trans5to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.3)[5,7]
  trans6to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.3)[6,7]
  trans5to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.3)[5,7]
}

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha3.3)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum3.3<-unlist(sum.year)
table3.3<-matrix((sum3.3)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table3.3, "/Users/n_a_abdallah/Desktop/GSR/table3.3.9states.April24.csv")


intalpha3.4<-alpha.nine(1,1,1,1,1,1,0.95,0.95,1)
trans6to7.m<-matrix(0,year2-year,(109-29))
trans5to7.m<-matrix(0,year2-year,(109-29))
trans6to7.f<-matrix(0,year2-year,(109-29))
trans5to7.f<-matrix(0,year2-year,(109-29))
prob1<-prob2<-fit<-sum.year<-inc.year<-vector("list")
for(l in 1:(109-29)){
  trans6to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.4)[6,7]
  trans5to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha3.4)[5,7]
  trans6to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.4)[6,7]
  trans5to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha3.4)[5,7]
}

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha3.4)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum3.4<-unlist(sum.year)
table3.4<-matrix((sum3.4)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table3.4, "/Users/n_a_abdallah/Desktop/GSR/table3.4.9states.April24.csv")
#combination prevention
intalpha4.1<-alpha.nine(0.50,1,0.75,1,0.50,0.75,0.90,0.90,1)

trans6to7.m<-matrix(0,year2-year,(109-29))
trans5to7.m<-matrix(0,year2-year,(109-29))
trans6to7.f<-matrix(0,year2-year,(109-29))
trans5to7.f<-matrix(0,year2-year,(109-29))
prob1<-prob2<-fit<-sum.year<-inc.year<-vector("list")
for(l in 1:(109-29)){
  trans6to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.f[1:3,l]<-pr.y.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[6,7]
  trans5to7.m[1:3,l]<-pr.y.males.alpha(28+l, 2014, intalpha)[5,7]
  trans6to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha4.1)[6,7]
  trans5to7.f[4:(year2-year),l]<-pr.y.alpha(28+l, 2018, intalpha4.1)[5,7]
  trans6to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha4.1)[6,7]
  trans5to7.m[4:(year2-year),l]<-pr.y.males.alpha(28+l, 2018, intalpha4.1)[5,7]
}

for(i in 1:(year2-year)){
  prob1[[1]]<-TPN.m(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  prob2[[1]]<-TPN.f(a=30,n=109-30,year,alpha=alpha.nine(1,1,1,1,1,1,1,1,1))
  fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha4.1)
  prob1[[i+1]]<-fit[[i]][1:80]
  prob2[[i+1]]<-fit[[i]][81:160]
  sum.year[[i]]<-fit[[i]][161:170]
}

sum4.1<-unlist(sum.year)
table4.1<-matrix((sum4.1)/1000000, nrow=46, ncol=10,byrow=T)
write.csv(table4.1, "/Users/n_a_abdallah/Desktop/GSR/table4.1.9states.may2.csv")

#intervention in terms of delays

pr.y.alpha2<-function(a,y,alpha){
  if(y<int.year){
    alpha=matrix(1,nrow=8,ncol=9)
  }else if(y>=int.year){
    alpha=alpha
  }
  if(y>2014){
    y=2014
  } else if (y<1933){
    y<-1933
  }
  else {y=y}
  p<-matrix(0,8,9)
  p[7,8]<-alpha[7,8]*0.167
  p[7,9]<-alpha[7,9]
  p[8,9]<-alpha[8,9]*(0.078)
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)}
        else if(65 <= a & a <= 75){
          p[1,2]<-alpha[1,2]*0.00105*exp(0.05596*a)
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)
        }
        else if(a>75){
          p[1,2]<-alpha[1,2]*0.07
          p[3,6]<-alpha[3,6]*k0[3,6]*exp(k1[3,6]*75)
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)
        }
        p[i,9]<-alpha[i,9]}
    }
  }
  for(i in 1:8){
    
      p[i,i]<-1-sum(p[i, 1:8]) 
    
  }
  return(p)}
prob.adj.am<-exp.am<-all.prob<-NULL
pr.rec.adj.am<-function(a,y,alpha){
  prob.adj.am<-rbind(c(pr.y.alpha2(a,y,alpha)[1,1:8],0),c(rep(0,8),1),
                        c(pr.y.alpha2(a,y,alpha)[3,1:8],0),c(rep(0,8),1),
                       c(rep(0,9)),c(rep(0,9)),c(rep(0,9)),c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.am)
}
prod.am<-vector("list")
lateage<-120
Amonset<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.am[[1]]<-pr.rec.adj.am(a,y,alpha)
    prod.am[[i+1]]<-(prod.am[[i]])%*%pr.rec.adj.am(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.am[1]<-a*(prod.am[[1]][1,2]+prod.am[[1]][1,4])
    exp.am[j+1]<-(exp.am[j])+((a+j)*(prod.am[[j+1]][1,2]+prod.am[[j+1]][1,4]))
    all.prob[1]<-(prod.am[[1]][1,2]+prod.am[[1]][1,4])
    all.prob[j+1]<-(all.prob[j])+((prod.am[[j+1]][1,2]+prod.am[[j+1]][1,4]))
  }
  return(c(all.prob[n],exp.am[n]))
}
age=30
Amonset(age,(lateage-age),2018,intalpha1.3)
#age of onset of MCI

prob.adj.mci<-exp.mci.m2<-all.prob.mci<-NULL
pr.rec.adj.mci<-function(a,y,alpha){
  prob.adj.mci<-rbind(cbind(pr.y.alpha2(a,y,alpha)[1:4,1:8],rep(0,4)),c(rep(0,8),1),c(rep(0,8),1),
                       c(rep(0,9)),c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.mci)
}
prod.mci<-vector("list")
lateage<-120
MCIonset<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.mci[[1]]<-pr.rec.adj.mci(a,y,alpha)
    prod.mci[[i+1]]<-(prod.mci[[i]])%*%pr.rec.adj.mci(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.mci.m2[1]<-a*(prod.mci[[1]][1,5]+prod.mci[[1]][1,6])
    exp.mci.m2[j+1]<-(exp.mci.m2[j])+((a+j)*(prod.mci[[j+1]][1,5]+prod.mci[[j+1]][1,6]))
    all.prob.mci[1]<-(prod.mci[[1]][1,5]+prod.mci[[1]][1,6])
    all.prob.mci[j+1]<-(all.prob.mci[j])+((prod.mci[[j+1]][1,5]+prod.mci[[j+1]][1,6]))
  }
  return(c(all.prob.mci[n],exp.mci.m2[n]))
}
age=30
MCIonset(age,(lateage-age),2018,alpha=intalpha2.4)

#age of onset of AD

prob.adj.ad<-exp.ad.m2<-all.prob.ad<-NULL
pr.rec.adj.ad<-function(a,y,alpha){
  prob.adj.ad<-rbind(cbind(pr.y.alpha2(a,y,alpha)[1:6,1:8],rep(0,6)),c(rep(0,8),1),
                      c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.ad)
}
prod.ad<-vector("list")
lateage<-120
ADonset<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.ad[[1]]<-pr.rec.adj.ad(a,y,alpha)
    prod.ad[[i+1]]<-(prod.ad[[i]])%*%pr.rec.adj.ad(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.ad.m2[1]<-a*(prod.ad[[1]][1,7])
    exp.ad.m2[j+1]<-(exp.ad.m2[j])+((a+j)*(prod.ad[[j+1]][1,7]))
    all.prob.ad[1]<-(prod.ad[[1]][1,7])
    all.prob.ad[j+1]<-(all.prob.ad[j])+((prod.ad[[j+1]][1,7]))
  }
  return(c(all.prob.ad[n],exp.ad.m2[n]))
}
age=30
ADonset(age,(lateage-age),2018,alpha=intalpha)

#using death not =0
#females
prob.adj.am.f<-exp.am<-all.prob<-NULL
pr.rec.adj.am.females<-function(a,y,alpha){
  prob.adj.am.f<-rbind(c(pr.y.alpha(a,y,alpha)[1,]),c(rep(0,8),1),
                     c(pr.y.alpha(a,y,alpha)[3,]),c(rep(0,8),1),
                     c(rep(0,9)),c(rep(0,9)),c(rep(0,9)),c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.am.f)
}
prod.am.f<-vector("list")
lateage<-130
Amonset.f<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.am.f[[1]]<-pr.rec.adj.am.females(a,y,alpha)
    prod.am.f[[i+1]]<-(prod.am.f[[i]])%*%pr.rec.adj.am.females(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.am[1]<-a*(prod.am.f[[1]][1,2]+prod.am.f[[1]][1,4])
    exp.am[j+1]<-(exp.am[j])+((a+j)*(prod.am.f[[j+1]][1,2]+prod.am.f[[j+1]][1,4]))
    all.prob[1]<-(prod.am.f[[1]][1,2]+prod.am.f[[1]][1,4])
    all.prob[j+1]<-(all.prob[j])+((prod.am.f[[j+1]][1,2]+prod.am.f[[j+1]][1,4]))
  }
  return(c(all.prob[n],exp.am[n]))
}
age=30
exp.amf1<-Amonset.f(age,(lateage-age),2018,intalpha)
exp.amf1.1<-Amonset.f(age,(lateage-age),2018,intalpha1.1)
exp.amf1.2<-Amonset.f(age,(lateage-age),2018,intalpha1.2)
exp.amf1.3<-Amonset.f(age,(lateage-age),2018,intalpha1.3)
exp.amf1.4<-Amonset.f(age,(lateage-age),2018,intalpha1.4)
write.csv(rbind(t(exp.amf1),t(exp.amf1.1),t(exp.amf1.2),t(exp.amf1.3),t(exp.amf1.4)), "/Users/n_a_abdallah/Desktop/GSR/exp.amf.may11.csv")
#males
prob.adj.am.m<-exp.am<-all.prob<-NULL
pr.rec.adj.am.males<-function(a,y,alpha){
  prob.adj.am.m<-rbind(c(pr.y.males.alpha(a,y,alpha)[1,]),c(rep(0,8),1),
                       c(pr.y.males.alpha(a,y,alpha)[3,]),c(rep(0,8),1),
                       c(rep(0,9)),c(rep(0,9)),c(rep(0,9)),c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.am.m)
}
prod.am.m<-vector("list")
lateage<-109
Amonset.m<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.am.m[[1]]<-pr.rec.adj.am.males(a,y,alpha)
    prod.am.m[[i+1]]<-(prod.am.m[[i]])%*%pr.rec.adj.am.males(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.am[1]<-a*(prod.am.m[[1]][1,2]+prod.am.m[[1]][1,4])
    exp.am[j+1]<-(exp.am[j])+((a+j)*(prod.am.m[[j+1]][1,2]+prod.am.m[[j+1]][1,4]))
    all.prob[1]<-(prod.am.m[[1]][1,2]+prod.am.m[[1]][1,4])
    all.prob[j+1]<-(all.prob[j])+((prod.am.m[[j+1]][1,2]+prod.am.m[[j+1]][1,4]))
  }
  return(c(all.prob[n],exp.am[n]))
}
age=30
exp.amm1<-Amonset.m(age,(lateage-age),2018,intalpha)
exp.amm1.1<-Amonset.m(age,(lateage-age),2018,intalpha1.1)
exp.amm1.2<-Amonset.m(age,(lateage-age),2018,intalpha1.2)
exp.amm1.3<-Amonset.m(age,(lateage-age),2018,intalpha1.3)
exp.amm1.4<-Amonset.m(age,(lateage-age),2018,intalpha1.4)
write.csv(rbind(t(exp.amm1),t(exp.amm1.1),t(exp.amm1.2),t(exp.amm1.3),t(exp.amm1.4)), "/Users/n_a_abdallah/Desktop/GSR/exp.amm.may11.csv")
#age of onset of MCI
#females
prob.adj.mci.f<-exp.mci.m2<-all.prob.mci<-NULL
pr.rec.adj.mci.f<-function(a,y,alpha){
  prob.adj.mci.f<-rbind(cbind(pr.y.alpha(a,y,alpha)[1:4,]),c(rep(0,8),1),c(rep(0,8),1),
                      c(rep(0,9)),c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.mci.f)
}
prod.mci.f<-vector("list")
MCIonset.f<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.mci.f[[1]]<-pr.rec.adj.mci.f(a,y,alpha)
    prod.mci.f[[i+1]]<-(prod.mci.f[[i]])%*%pr.rec.adj.mci.f(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.mci.m2[1]<-a*(prod.mci.f[[1]][1,5]+prod.mci.f[[1]][1,6])
    exp.mci.m2[j+1]<-(exp.mci.m2[j])+((a+j)*(prod.mci.f[[j+1]][1,5]+prod.mci.f[[j+1]][1,6]))
    all.prob.mci[1]<-(prod.mci.f[[1]][1,5]+prod.mci.f[[1]][1,6])
    all.prob.mci[j+1]<-(all.prob.mci[j])+((prod.mci.f[[j+1]][1,5]+prod.mci.f[[j+1]][1,6]))
  }
  return(c(all.prob.mci[n],exp.mci.m2[n]))
}
age=30
exp.mcif1<-MCIonset.f(age,(lateage-age),2018,intalpha)
exp.mcif2.1<-MCIonset.f(age,(lateage-age),2018,intalpha2.1)
exp.mcif2.2<-MCIonset.f(age,(lateage-age),2018,intalpha2.2)
exp.mcif2.3<-MCIonset.f(age,(lateage-age),2018,intalpha2.3)
exp.mcif2.4<-MCIonset.f(age,(lateage-age),2018,intalpha2.4)
write.csv(rbind(t(exp.mcif1),t(exp.mcif2.1),t(exp.mcif2.2),t(exp.mcif2.3),t(exp.mcif2.4)), "/Users/n_a_abdallah/Desktop/GSR/exp.mcif.may11.csv")
#males
prob.adj.mci.m<-exp.mci.m2<-all.prob.mci<-NULL
pr.rec.adj.mci.m<-function(a,y,alpha){
  prob.adj.mci.m<-rbind(cbind(pr.y.males.alpha(a,y,alpha)[1:4,]),c(rep(0,8),1),c(rep(0,8),1),
                        c(rep(0,9)),c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.mci.m)
}
prod.mci.m<-vector("list")
MCIonset.m<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.mci.m[[1]]<-pr.rec.adj.mci.m(a,y,alpha)
    prod.mci.m[[i+1]]<-(prod.mci.m[[i]])%*%pr.rec.adj.mci.m(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.mci.m2[1]<-a*(prod.mci.m[[1]][1,5]+prod.mci.m[[1]][1,6])
    exp.mci.m2[j+1]<-(exp.mci.m2[j])+((a+j)*(prod.mci.m[[j+1]][1,5]+prod.mci.m[[j+1]][1,6]))
    all.prob.mci[1]<-(prod.mci.m[[1]][1,5]+prod.mci.m[[1]][1,6])
    all.prob.mci[j+1]<-(all.prob.mci[j])+((prod.mci.m[[j+1]][1,5]+prod.mci.m[[j+1]][1,6]))
  }
  return(c(all.prob.mci[n],exp.mci.m2[n]))
}
age=30
exp.mcim1<-MCIonset.m(age,(lateage-age),2018,intalpha)
exp.mcim2.1<-MCIonset.m(age,(lateage-age),2018,intalpha2.1)
exp.mcim2.2<-MCIonset.m(age,(lateage-age),2018,intalpha2.2)
exp.mcim2.3<-MCIonset.m(age,(lateage-age),2018,intalpha2.3)
exp.mcim2.4<-MCIonset.m(age,(lateage-age),2018,intalpha2.4)
write.csv(rbind(t(exp.mcim1),t(exp.mcim2.1),t(exp.mcim2.2),t(exp.mcim2.3),t(exp.mcim2.4)), "/Users/n_a_abdallah/Desktop/GSR/exp.mcim.may11.csv")
#age of onset of AD
#females
prob.adj.ad.f<-exp.ad.m2<-all.prob.ad<-NULL
pr.rec.adj.ad.f<-function(a,y,alpha){
  prob.adj.ad.f<-rbind(cbind(pr.y.alpha(a,y,alpha)[1:6,]),c(rep(0,8),1),
                     c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.ad.f)
}
prod.ad.f<-vector("list")
ADonset.f<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.ad.f[[1]]<-pr.rec.adj.ad.f(a,y,alpha)
    prod.ad.f[[i+1]]<-(prod.ad.f[[i]])%*%pr.rec.adj.ad.f(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.ad.m2[1]<-a*(prod.ad.f[[1]][1,7])
    exp.ad.m2[j+1]<-(exp.ad.m2[j])+((a+j)*(prod.ad.f[[j+1]][1,7]))
    all.prob.ad[1]<-(prod.ad.f[[1]][1,7])
    all.prob.ad[j+1]<-(all.prob.ad[j])+((prod.ad.f[[j+1]][1,7]))
  }
  return(c(all.prob.ad[n],exp.ad.m2[n]))
}
age=30
ADonset.f(age,(lateage-age),2018,alpha=intalpha3.3)
exp.adf1<-ADonset.f(age,(lateage-age),2018,intalpha)
exp.adf3.1<-ADonset.f(age,(lateage-age),2018,intalpha3.1)
exp.adf3.2<-ADonset.f(age,(lateage-age),2018,intalpha3.2)
exp.adf3.3<-ADonset.f(age,(lateage-age),2018,intalpha3.3)
exp.adf3.4<-ADonset.f(age,(lateage-age),2018,intalpha3.4)
write.csv(rbind(t(exp.adf1),t(exp.adf3.1),t(exp.adf3.2),t(exp.adf3.3),t(exp.adf3.4)), "/Users/n_a_abdallah/Desktop/GSR/exp.adf.may11.csv")
#males
prob.adj.ad.m<-exp.ad.m2<-all.prob.ad<-NULL
pr.rec.adj.ad.m<-function(a,y,alpha){
  prob.adj.ad.m<-rbind(cbind(pr.y.males.alpha(a,y,alpha)[1:6,]),c(rep(0,8),1),
                       c(rep(0,9)),c(rep(0,9)))
  return(prob.adj.ad.m)
}
prod.ad.m<-vector("list")
ADonset.m<-function(a,n=(lateage-a),y,alpha){
  ad.m<-NULL
  for(i in 1:(n-1)){
    prod.ad.m[[1]]<-pr.rec.adj.ad.m(a,y,alpha)
    prod.ad.m[[i+1]]<-(prod.ad.m[[i]])%*%pr.rec.adj.ad.m(a+(i+1)-1,y,alpha)
  }
  for(j in 1:(lateage-a-1)){
    exp.ad.m2[1]<-a*(prod.ad.m[[1]][1,7])
    exp.ad.m2[j+1]<-(exp.ad.m2[j])+((a+j)*(prod.ad.m[[j+1]][1,7]))
    all.prob.ad[1]<-(prod.ad.m[[1]][1,7])
    all.prob.ad[j+1]<-(all.prob.ad[j])+((prod.ad.m[[j+1]][1,7]))
  }
  return(c(all.prob.ad[n],exp.ad.m2[n]))
}
age=30
ADonset.m(age,(lateage-age),2018,alpha=intalpha)
exp.adm1<-ADonset.m(age,(lateage-age),2018,intalpha)
exp.adm3.1<-ADonset.m(age,(lateage-age),2018,intalpha3.1)
exp.adm3.2<-ADonset.m(age,(lateage-age),2018,intalpha3.2)
exp.adm3.3<-ADonset.m(age,(lateage-age),2018,intalpha3.3)
exp.adm3.4<-ADonset.m(age,(lateage-age),2018,intalpha3.4)
write.csv(rbind(t(exp.adm1),t(exp.adm3.1),t(exp.adm3.2),t(exp.adm3.3),t(exp.adm3.4)), "/Users/n_a_abdallah/Desktop/GSR/exp.adm.may11.csv")

#lifetime risk 
#males
lifetimerisk<-matrix(0,nrow=8,ncol=6)
lifetimerisk[1,]<-AR.m(30,(110-30),2017)[1:6,7]
lifetimerisk[2,]<-AR.m(60,(110-60),2017)[1:6,7]
lifetimerisk[3,]<-AR.m(65,(110-65),2017)[1:6,7]
lifetimerisk[4,]<-AR.m(70,(110-70),2017)[1:6,7]
lifetimerisk[5,]<-AR.m(75,(110-75),2017)[1:6,7]
lifetimerisk[6,]<-AR.m(80,(110-80),2017)[1:6,7]
lifetimerisk[7,]<-AR.m(85,(110-85),2017)[1:6,7]
lifetimerisk[8,]<-AR.m(90,(110-90),2017)[1:6,7]
write.csv(lifetimerisk,"/Users/n_a_abdallah/Desktop/GSR/lifetime.csv")
a=90
y=2017
n=a-30  
y2<-y-n-1
prevalence<-NULL
prod<-vector("list")
for(i in 1:(n)){
  prod[[1]]<-TP.m(30,y2,intalpha)
  prod[[i+1]]<-(prod[[i]])%*%TP.m(30+i, y2+i,intalpha)}
for(i in 1:6){
  prevalence[i]<-(lifetimerisk[8,i]*prod[[n+1]][1,i])/sum(prod[[n+1]][1,1:6])
}
sum(prevalence)

#females
lifetimerisk.f<-matrix(0,nrow=8,ncol=6)
lifetimerisk.f[1,]<-AR.f(30,(110-30),2017)[1:6,7]
lifetimerisk.f[2,]<-AR.f(60,(110-60),2017)[1:6,7]
lifetimerisk.f[3,]<-AR.f(65,(110-65),2017)[1:6,7]
lifetimerisk.f[4,]<-AR.f(70,(110-70),2017)[1:6,7]
lifetimerisk.f[5,]<-AR.f(75,(110-75),2017)[1:6,7]
lifetimerisk.f[6,]<-AR.f(80,(110-80),2017)[1:6,7]
lifetimerisk.f[7,]<-AR.f(85,(110-85),2017)[1:6,7]
lifetimerisk.f[8,]<-AR.f(90,(110-90),2017)[1:6,7]
write.csv(lifetimerisk.f,"/Users/n_a_abdallah/Desktop/GSR/lifetime.f.csv")
a=90
y=2017
n=a-30  
y2<-y-n-1
prevalence<-NULL
prod<-vector("list")
for(i in 1:(n)){
  prod[[1]]<-TP.f(30,y2,intalpha)
  prod[[i+1]]<-(prod[[i]])%*%TP.f(30+i, y2+i,intalpha)}
for(i in 1:6){
  prevalence[i]<-(lifetimerisk.f[8,i]*prod[[n+1]][1,i])/sum(prod[[n+1]][1,1:6])
}
sum(prevalence)

#10 year risk 
#males
tenyearrisk<-matrix(0,nrow=8,ncol=6)
tenyearrisk[1,]<-AR.m(30,10,2017)[1:6,7]
tenyearrisk[2,]<-AR.m(60,10,2017)[1:6,7]
tenyearrisk[3,]<-AR.m(65,10,2017)[1:6,7]
tenyearrisk[4,]<-AR.m(70,10,2017)[1:6,7]
tenyearrisk[5,]<-AR.m(75,10,2017)[1:6,7]
tenyearrisk[6,]<-AR.m(80,10,2017)[1:6,7]
tenyearrisk[7,]<-AR.m(85,10,2017)[1:6,7]
tenyearrisk[8,]<-AR.m(90,10,2017)[1:6,7]
write.csv(tenyearrisk,"/Users/n_a_abdallah/Desktop/GSR/tenyearrisk.csv")
a=90
y=2017
n=a-30  
y2<-y-n-1
prevalence<-NULL
prod<-vector("list")
for(i in 1:(n)){
  prod[[1]]<-TP.m(30,y2,intalpha)
  prod[[i+1]]<-(prod[[i]])%*%TP.m(30+i, y2+i,intalpha)}
for(i in 1:6){
  prevalence[i]<-(tenyearrisk[8,i]*prod[[n+1]][1,i])/sum(prod[[n+1]][1,1:6])
}
sum(prevalence)

#females
tenyearrisk.f<-matrix(0,nrow=8,ncol=6)
tenyearrisk.f[1,]<-AR.f(30,10,2017)[1:6,7]
tenyearrisk.f[2,]<-AR.f(60,10,2017)[1:6,7]
tenyearrisk.f[3,]<-AR.f(65,10,2017)[1:6,7]
tenyearrisk.f[4,]<-AR.f(70,10,2017)[1:6,7]
tenyearrisk.f[5,]<-AR.f(75,10,2017)[1:6,7]
tenyearrisk.f[6,]<-AR.f(80,10,2017)[1:6,7]
tenyearrisk.f[7,]<-AR.f(85,10,2017)[1:6,7]
tenyearrisk.f[8,]<-AR.f(90,10,2017)[1:6,7]
write.csv(tenyearrisk.f,"/Users/n_a_abdallah/Desktop/GSR/tenyearriskf.csv")
a=90
y=2017
n=a-30  
y2<-y-n-1
prevalence<-NULL
prod<-vector("list")
for(i in 1:(n)){
  prod[[1]]<-TP.f(30,y2,intalpha)
  prod[[i+1]]<-(prod[[i]])%*%TP.f(30+i, y2+i,intalpha)}
for(i in 1:6){
  prevalence[i]<-(tenyearrisk.f[8,i]*prod[[n+1]][1,i])/sum(prod[[n+1]][1,1:6])
}
sum(prevalence)
