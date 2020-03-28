

#initialize the variables 
N <- 25000         # N = sample size
S <- rep(0,N)      # S = susceptible 
I <- rep(0,N)      # I = infectious
R <- rep(0,N)      # R = removed/recovered
time <- rep(0,N)
nwcas<-rep(0,N)   #nwcas = new cases of the disease 
S[1]<-20000
I[1]<-5
R[1]<-80000
beta<-0.000023     # infection speed 
alpha<-0.337       # recovery time  => 1/alpha = 3 days to recover  

myfct<-function() {
  for(i in 1:(N-1))  
  {
    if(I[i]!=0) {    # if at least 1 person is infectious
      X <-runif(1)   # X~U[0,1] - X follows a uniform distribution 
      lambda = beta*S[i]*I[i]+alpha*I[i]
      # Y = Markov proces / time it takes for the process of ‘S’ → ‘I’ or ‘I’ → ‘R’
      Y <- -log(X)/lambda
      time[i+1]<-time[i] + Y  # time_new = time_old + Y 
      probaSI<-(beta*S[i]*I[i])/(beta*S[i]*I[i]+alpha*I[i])
      #probaIR<-(alpha*I[i])/(beta*S[i]*I[i]+alpha*I[i])
      
      #If X < P[‘S’ → ‘I’] then (S,I,R) → (S-1,I+1,R)
      if (X<=probaSI) {
        S[i+1]<-S[i]-1    #S_new = S_old - 1 
        I[i+1]<-I[i]+1    #I_new = I_old + 1 
        R[i+1]<-R[i]      #R_new = R_old 
        nwcas[i+1]<-1     # 1 new case 
      }
      #Else (S,I,R) → (S,I-1,R+1) 
      else {
        S[i+1]<-S[i]    #S_new = S_old 
        I[i+1]<-I[i]-1  #I_new = I_old - 1 
        R[i+1]<-R[i]+1  #R_new = R_old + 1
        nwcas[i+1]<-0   #no new cases 
      }
    }
  }
  
  #round the numbers of the time value - 1 unit of time = 1 day
  timeclass<-as.numeric(cut(time,breaks=seq(0,max(time),1),include.lowest=TRUE))-1
  #create a contingency table with the time values and the cases (0 = no new case and 1 = a new case)
  contingency<-table(timeclass,nwcas)    
  epidemic<-contingency[,2]     # used to graph the nwcas vs time incidence curve 
  
  items_returned <- list()
  items_returned[["epidemic"]] <- epidemic
  items_returned[["nwcas"]] <- nwcas
  return (items_returned)
}


# Run simulation 60 times and assign nwcas to a data.frame column each time

output <- matrix(nrow=N, ncol=60)

for(i in 1:60) {
  
  items <- myfct()
  
  output[ , i] <- items$nwcas
  
  
}



