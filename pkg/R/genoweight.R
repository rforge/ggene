genoweight <-function(X, genotypes)
{
n<-dim(X$tab)[1] ; genotype.vec <- genotypes
#- weighting for recurrent genotypes
genotypes.vect<-c(0,sort(unique(as.vector(genotype.vec)))) #sert simplement au calcul suivant
genotypes.counts<-hist(genotype.vec, plot=FALSE, breaks=genotypes.vect)$counts
genotypes.weight<- 1/genotypes.counts
genotype.vect.short <- sort(unique(as.vector(genotype.vec)))
mat.weight.genotype <- matrix(NA, nrow=n, ncol=n) # matrice des pondrations des diff gnotypes
  for (i in 1:(n)){
    for (j in 1:(n)){
      if(j<i){
	    if(genotype.vec[i]==genotype.vec[j]) mat.weight.genotype[i,j] <- 0
	    else
	      mat.weight.genotype[i,j]<- genotypes.weight[genotype.vec[i]]*genotypes.weight[genotype.vec[j]]}
      }
}
return(as.dist(mat.weight.genotype))}