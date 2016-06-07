genocount <-function(X)
{
loci.dat <- X$tab
vec.geno.num<-duplicated(loci.dat) # a est un vecteur
vec.geno.num
genotype.vec<-vector(length=length(duplicated(loci.dat)))
comp<-1
for (i in 1:(length(vec.geno.num))){
  if(vec.geno.num[i]==FALSE) {genotype.vec[i] <- comp ; comp<-comp+1}}
genotype.vec
for (i in 1:(length(vec.geno.num))){
  if(vec.geno.num[i]==TRUE) {
  id<-FALSE
        for (j in 1:(dim(loci.dat)[1])){
           comp<-0
           for (k in 1:(dim(loci.dat)[2])){
              if(loci.dat[i,k]==loci.dat[j,k]) {comp<-comp+1} else comp<-comp
              }
           if(comp==dim(loci.dat)[2]) {genotype.vec[i]<-genotype.vec[j];break}
           }}}
l<-list(vec=genotype.vec, n=length(unique(genotype.vec)))
return(l)}

