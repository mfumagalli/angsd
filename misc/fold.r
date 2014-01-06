source("http://www.popgen.dk/thorfinn/src/supersrc.R")


n <- as.numeric(system("wc -l taj/norm.sfs.pos|cut -f1 -d \" \"",intern=T))
a<-exp(scan("taj/norm.sfs.em.ml"))
b <- c(0,read.table("taj/msout.sfs")[,1],0)
barplot(rbind(a*n,b),be=T)

pair(b,20)/100
sum(b)/a1(40)/100
