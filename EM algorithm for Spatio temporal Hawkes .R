## EM algorithm 
## parameters for mu(.) - N(m,S)
## parameters in g(.) - gamma and Sigma 
occurances <- occurances[-1,]  #t_i, eta_i values

#Kernel (Gaussian with Sigma - Maholanonbish)
library(mvtnorm)
K <- function(eta, Sigma, bandwidth){
  return(dmvnorm(eta/bandwidth, sigma = Sigma)/bandwidth^length(eta))
}


#given a P 
#sample offspring and background events
sampled_background <- function(P){
  return(which(runif(n) < diag(P)))
}

#unnormalised background density
#u(eta) which is mu(.) = alpha. u(.), i.e. gives the base denisty  
# u <- function(eta, back_event){
#   
#   vals = numeric(length(back_event))
#   ## put a criteria of length(vals >0 later)
#   for( i in 1:length(vals)){
#     ##d_i =choose from 15th nearest neighbour criteria in the paper
#     d_i = sort(D[i,])[15]
#     vals[i] = K(eta-occurances[back_event[i],-1],Sigma,d_i)
#   }
#   return(mean(vals))  
# }

#base spatial density
u <-function(eta, back_event){
  
  events = occurances[back_event,-1]
  l = length(back_event)
  if(l==0) print("NO background events!!!")
  
  if (l>1){
  m = colMeans(events)
  s = cov(events)
  return(dmvnorm(eta,mean=m, sigma = s))
  }
  
}

## base intensity time exponential density
v <- function(t, back_event){
  times= occurances[back_event,1]
  l = length(back_event)
  if(l==0) print("NO background events!!!")
  if(l >1){
  return(dexp(t, rate = 1/mean(times)))
  }
}

#background density
mu<-function(t,eta, alpha, back_event){
  return(alpha*u(eta,back_event)*v(t,back_event))
}


#load your dataset
## compute ||eta_i-eta_j|| values by distance matrix from occurances[,-1]
D = dist(occurances[,-1])
D = as.matrix(D)
#View(D) #use this to find nearest neighbours by rank 

## n = number of events 
#occurance is the dataset
# #initialise P0
# n = dim(occurances)[1]
# n
# init = matrix(0,n,n)
# init[1,1] =1
# for(i in 2:n){
#   init[i,1:i] = rep(1/i,i)
#  # init[i,1:(i-1)] = rep(1/(2*(i-1)), i-1)
#  # init[i,i] = 0.5
#  }

## start iterations from here ------------------------------
#########################################################
#update P^(k+1) from P^k 
## P = P^k as of now

n = dim(occurances)[1]
## add some randomisation, might help reaching better performance
m0 = rnorm(d,0,1)/10 ##initialise m
S0 = diag(runif(d,0,2))##initialise S
phi_0= 1/mean(occurances[,1])
alpha_0 = 2
beta_0 = 2
Sigma_0 = cov(occurances[,-1]) ##initialise Sigma
gamma_0 = 0.3

#initialise P_0
P0 = matrix(0,n,n)
for(i in 1:n){
  for(j in 1:i){
    if(j!=i) P0[i,j] = beta_0*gaussian_density(occurances[i,2:(d+1)]-occurances[j,2:(d+1)], sigma= Sigma_0)*dexp(occurances[i,1]-occurances[j,1], rate = gamma_0) #wrong in the gaps term #trigger terms
    if(i==j)  P0[i,i] = alpha_0*gaussian_density(occurances[i,2:(d+1)], mean=m0, sigma = S0)*dexp(occurances[i,1], rate = phi_0)  ## mu
  }
  
  P0[i,] = P0[i,]/sum(P0[i,])    
  if(i%%50==0) print(paste0("row ", i, " done"))
}

###########################################
View(P0)
hist(diag(P0)) ## sanity check 
##shouldn't have too many big entries near 1 in the diagonal 


#sample for this iteration
t1=Sys.time()
P=P0
S= S0
m = m0
Sigma = Sigma_0
phi = phi_0
gamma = gamma_0

## how many EM updates ?
iter_max = 20
change = numeric(iter_max)

for(r in 1:iter_max){
back_event = sampled_background(P)
# #resample if not enough samples
# while (length(back_event)<2 || length(back_event)> n-2){
#   back_event = sampled_background(P)
#   print("resampling")
# }

events = c(1:n)
offsping_event = events[-back_event] 

## update m and S
base_locs = occurances[back_event,-1]
l = length(back_event)
if(l==0){ 
  print("No background events!!!")
  m_new = as.numeric(base_locs) #keep m
  S_new = Sigma ## keep S as it is 
}
if(l ==1){
  m_new = as.numeric(base_locs) #new m for mean of base spatial component 
  S_new = Sigma ## keep S as it is 
}

if (l>1){
  m_new = colMeans(base_locs) #new m for mean of base spatial component 
  S_new = cov(base_locs) ## new S for covaraince matrix for base spatial compoennet 
}


#update alpha, beta, gamma, phi 
alpha_new = length(back_event)
w_gaps = matrix(0,n,n)
s=0 ## \sum_{i \in offspring}{ j <i} P_{ij}
for(i in offsping_event)for(j in 1:i){
  if(i>j) {
    w_gaps[i,j] =(occurances[i,1]-occurances[j,1])*P[i,j] 
    s = s+P[i,j]
  }
}

new_gamma = s/sum(w_gaps)
#beta_new = length(offsping_event)/(n-sum(exp(-new_gamma*(T-occurances[,1])))) ##unnecessarily complex 
beta_new = length(offsping_event)/n
phi_new = 1/mean(occurances[back_event,1])

#update sigma
Sigma_new = matrix(0,d,d)
for(i in offsping_event)for(j in 1:i){
  term=(occurances[i,-1]-occurances[j,-1])%*%t(occurances[i,-1]-occurances[j,-1])  #check teh dimensions
  Sigma_new = Sigma_new + term*P[i,j] 
}

## trying to update P
dummy = matrix(0,n,n)
for(i in 1:n){
  for(j in 1:i){
     if(j!=i) dummy[i,j] = beta_new*gaussian_density(occurances[i,-1]-occurances[j,-1], sigma= Sigma_new)*dexp(occurances[i,1]-occurances[j,1], rate = new_gamma) #g(t_i-t_j, eta_i-eta_j)
     if(i==j)  dummy[i,i] = alpha_new*gaussian_density(occurances[i,-1], mean = m_new,sigma = S_new)*dexp(occurances[i,1], rate = phi_new)  ## mu(t_i, eta_i)
  }
  
 dummy[i,] = dummy[i,]/sum(dummy[i,])    
 if(i%%50==0) print(paste0("row ", i, " done"))
 
 ##parameter printing 
 ## current updates
 S= S_new
 m = m_new
 Sigma = Sigma_new
 phi = phi_new
 
}


change[r] = norm(dummy-P, type="F") ## calculate ||P^(k+1)-P^k||
P = dummy ## new P updated
View(P)
beepr::beep()
print(paste0("Iteration ",r, " in EM algo done"))
print(paste0("time taken till iteration ", r, " is ", Sys.time()-t1))
print(paste0("sum of  the diag P, aka sum P_ii's is ", sum(diag(P))))



}


## Most P_ij updates seem to be extremely small and one of them is turning about to be 1? am I doing a mistake? check!!!
## usual suspects - check how the MBN generation is working, do we need to do Cholesky decomposition by hand? 
## How about changing the kernel bandwidth? it seems too many probabilities are too small, so maybe let us use a bigger bandfwidth? 
#how large are the P_ii's?

hist(diag(P))
hist(P)
plot(change, ty="l", col="red", lwd=2)
