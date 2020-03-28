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

nbDays<-200
tnwcas<-rep(0,nbDays)

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
  
  for (i in 1:nbDays) {
    tnwcas[i]<-sum(nwcas[time > i-1 & time < i]) # 0 < time < 1 , time = i, tnwcas[i] = nwcas in terms of i
  }
  
  return (tnwcas)
}



nbSimul<-60
mat <- matrix(NA, nrow  = nbDays, ncol = nbSimul)





plot(NA,type = "n", xlab = "Time (days)", ylab = "Incidence (per 100 000 inhabitants)", 
     xlim = c(0,200), ylim = c(0,400), main = "Incidence curves of 60 epidemics")

for (i in 1:nbSimul) {
  mat[,i]<-myfct()               #mat[,epidemicnumber]
  lines(mat[,i], type = "l", col = i)
}



