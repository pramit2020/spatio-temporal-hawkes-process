#mu(s) \prop e^{-||s||^2/2}, base intensity
# beta = time decay rate
# Sigma - dxd matrix
# s ( spatial component) in R^d 
# element l in the list computes l^th generation of offsprings
# list[[l]] is a matrix of size N_lx(1+d) 

## parameters to generate 
d=10
T = 15
gamma = 0.5
B = mvtnorm::rmvnorm(d,mean=rep(0,d))
Sigma = B%*%t(B)/5
alpha = 3
beta = 1.5
phi = 4

#generate the base intensity events
ST_hawkes = list()
N0 = rpois(1, alpha*(1-exp(-phi*T))) 
G0 <- matrix(0, nrow=N0, ncol = 1+d) #first column is times (s,t)
G0[,1]= rexp(N0,rate = phi)
for( i in 1:N0){
 G0[i,-1] = mvtnorm::rmvnorm(1, mean=rep(0,d), sigma = Sigma) 
}

accepted <- which(G0[,1]<=T)
ST_hawkes[[1]] <- G0[accepted,]

#generate offspring levels
m = beta*(1-exp(-gamma*T)) # mean number of offspring 
print(m)
# level l to l+1

# a function which generates offspringfs given a parent
offspring <- function(v){ #give a vector of form (t,s) and get a matrix of its offsprings 
  
  N = rpois(1,m)
  if(N==0) return("no babies")
  
  if(N>0){
  O <- matrix(0,N,1+d)
  O[,1] <- v[1]+ rexp(N, rate = gamma)
  
  for(i in 1:N){
  O[i,2:(d+1)] = v[2:(d+1)]+ mvtnorm::rmvnorm(1, mean=rep(0,d), sigma = Sigma) 
  }
  
  accepted <- which(O[,1]<=T)
  return(O[accepted,])
  }
}

#Number of generations M
l=1
p=1
c=0 # counter for events 
while(p>0 & c<=20000){
  
  parents = ST_hawkes[[l]] #matrix
  
  num_parents = dim(parents)[1]
  c=c+num_parents
  
  if(c>20000) break
  
  if(num_parents==0) break 
  
  kids <- matrix(0,1,1+d) #dummy 0 row
  for(i in 1:num_parents){
    
    vector = as.numeric(parents[i,])
    temp = offspring(vector)
    if(typeof(temp)=="double") kids <- rbind(kids, temp)
  }
  
  hello <- kids[-1,]
  
  if(sum(dim(hello))==0) p=0
  
  ST_hawkes[[l+1]] <- hello
  print(paste0("We are at generation ", l))
  print(paste0("Added ", dim(hello)[1], " offsprings"))
  l=l+1
}

## straighten things out
g = length(ST_hawkes)
occurances = matrix(0,1,d+1)
for( i in 1:g){
  Iwantthese= ST_hawkes[[i]]
  occurances <- rbind(occurances, Iwantthese)
}

Myoutoput <- occurances[-1,] #get rid of the dummy first row 
## sort based on time 
Myoutoput <- Myoutoput[order(Myoutoput[,1],decreasing=FALSE),]

rownames(Myoutoput) <-c()
colnames(Myoutoput) <- c("time", paste0("eta", c(1:d)))
print(paste0("Total number of events is  ",c))
View(Myoutoput)
