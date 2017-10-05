library(deSolve)
library(ggplot2)

phage_dede<- function(time, state, parameters){
  ### Parameters### 
  mu<- parameters["mu"] 
  
  carrying_capacity<- parameters["carrying_capacity"]
  
  D<- parameters["D"] ### dilution rate 
  K_i<- parameters["K_i"]
  C <- parameters["C"] ###carrying capacity
  tau<- parameters["tau"] ###delay time tau
  b<- parameters["b"] ###burst size
  d_p<- parameters["d_p"] ###phage decay rate
  
  
  ### States ###
  X_s<- state["X_s"] ##concentration of bacteria
  P<- state["P"] ### Phage Concentration 
  
  if (time<tau){
    X_lag<-0#X_s
    P_lag<-0#P
  } else{
    lagStates<- lagvalue(time-tau)
    X_lag<- X_s#lagStates[1]
    P_lag<- P#lagStates[2]
  }
  #print(X_lag)
  
  dX_s<- mu * X_s * (1- X_s/C) - (D *X_s)  - (b*K_i * X_s * P)
  dP <- -D* P - (K_i* X_s * P) + (b*K_i* X_lag * P_lag) - d_p* P 
  
  
  return(list(c(dX_s, dP)))
}

theta<- c(D= 0.4, K_i= 6*10^-4,mu= 1.04, C= 100, tau= 1, b= 90, d_p= 0)
initState<- c(X_s= 25, P= 2)
time<- seq(from= 0, to= 150, by= 0.1)


trajectory.df <- data.frame(dede(y= initState, times= time, parms= theta, func= phage_dede))
Conc<- trajectory.df$X_s
phageConc<- trajectory.df$P

ggplot(data= trajectory.df, aes(x= time)) + 
  geom_line(aes(y= Conc, colour= "Yeast Concentration"))+ 
  geom_line(aes(y= phageConc, colour= "YeaBoi Concentration"))

