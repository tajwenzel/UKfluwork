

#data 
proposed.contact.ids <<- current.contact.ids
if (runif(1,0,1) < 0.1) {
  rs <<- round(runif(2,1,length(proposed.contact.ids)))
  proposed.contact.ids[rs[1]] <<- rs[2]
}
# Resample contact ids. T

p.contacts<-NULL
for (jj in 1:out)
{
  p.contacts[[jj]]<-contact.matrix(as.matrix(polymod[contact.ids[[jj]],]),
                             age_sizes[,1], agegrouplimits)
}

trial<-array(unlist(p.contacts), c(7,7,out))
iter.plot.set<-as.data.frame(t(combn(1:7, m=2, FUN = NULL, simplify = TRUE)))
#diagnol<-c(1,1,2,2,3,3,4,4,5,5,6,67,7)
#iter.plot.set<-rbind(iter.plot.set,diagnol)

set.mat<-matrix(
            data=c(1, 1, 1, 1, 1, 1, 1,
                   0, 1, 1, 1, 1, 1, 1,
                   0, 0, 1, 1, 1, 1, 1,
                   0, 0, 0, 1, 1, 1, 1,
                   0, 0, 0, 0, 1, 1, 1,
                   0, 0, 0, 0, 0, 1, 1,
                   0, 0, 0, 0, 0, 0, 1), nrow=7, ncol=7, byrow=TRUE)

if (set.mat==1, TRUE) 
  { 
  par(mfrow=c(floor(1+(ncol(set.mat)/2)), nrow(set.mat)));
  iter<-as.data.frame(which(set.mat==1, TRUE));
  
  for (qq in 1:max(iter$row))
    {
    for (pp in 1:max(iter$col)) 
      {hist(trial[qq,pp,],breaks=20,main='ctmat')}
  }
}


#lapply(trial[iter.plot.set$V1,iter.plot.set$V2,],1,FUN=hist)
#trial[iter.plot.set$V1[1],iter.plot.set$V2[1],]
