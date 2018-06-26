
connectance=.5
n_plants=15
n_gallers=20
links=c(rep(1,connectance*n_plants*n_gallers),rep(0,(1-connectance)*n_plants*n_gallers))

# Some matrices to store 
A_props=matrix(nrow=100,ncol=3)
colnames(A_props)=c("L_per_plant","L_per_galler","NODF")
B_props=matrix(nrow=10000,ncol=6)
colnames(B_props)=c("N_plants","N_gallers","C","L_per_plant","L_per_galler","NODF")
C_props=matrix(nrow=1000000,ncol=6)
colnames(C_props)=c("N_plants","N_gallers","C","L_per_plant","L_per_galler","NODF")


# Create a matrix with links randomly assigned to give our desired connectance
A=matrix(nrow=n_plants,ncol=n_gallers,data=sample(links,replace=FALSE))

# Assume that we have 100 sampling days.
# The occurrence of each interaction during sampling depends on the process uncertainty of the interaction
# Assume that process uncertainties are uniformly distributed.
# The occurrance of an interaction during any sampling day is drawn from 
# a binomial distribution with n=100, size=



# Assume that each interaction has a process certainty of 0.8 and we have 100 sampling occasions
# I.e., the interaction truly occurs during 80% of sampling occasions.
# Occurrence of the interaction during any sampling day is a binom
process=0.5
B=matrix(nrow=n_plants,ncol=n_gallers)
for(plant in 1:nrow(A)){
	for(galler in 1:ncol(A)){
		B[plant,galler]=A[plant,galler]*rbinom(1,1,process)
	}
}


# Assume that all interactions are detectable in 60% of the samples in which they occur
# Assume that each pair is observed co-occurring a random number of times between 0 and 100
# The number of detected interactions is a Bernoulli process depending on the number of 
# observed co-occurrences
detection=0.5
C=matrix(nrow=n_plants,ncol=n_gallers)
for(plant in 1:nrow(B)){
	for(galler in 1:ncol(B)){
		cooccurrences=round(runif(n=1,min=0,max=1))
		detections_if_occurring=rbinom(1,cooccurrences,detection)
		if(B[plant,galler]==1 & detections_if_occurring>0){
			C[plant,galler]=1
		} else {
			C[plant,galler]=0
		}
	}
}

# Compare A, B, and C.
# Links, plants and gallers with at least 1 link, C, ...


