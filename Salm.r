#-----------------------data-----------------------

Ndoses = 6
Nplates = 3
y = matrix(c(15, 16, 16, 27, 33, 20, 21, 18, 26, 41, 38, 27, 29, 
             21, 33, 60, 41, 42), ncol = 3, byrow = TRUE)
x = c(0, 10, 33, 100, 333, 1000)

#--------------------------MH------------------------

mu_fn=function(alpha, beta, gamma, x, lambda){
  log.mu = alpha + beta * log(x+10) + gamma*x + lambda
  return(exp(log.mu))
}

salm <- function(nchain, init.abg,init.lambda, init.tau,
                 prop.sd.abg, prop.sd.lambda){
  
  chain.abg = matrix(NA, nchain+1, Nplates)
  colnames(chain.abg) = c('alpha', 'beta', 'gamma')
  
  chain.lambda = array(NA, dim = c(nchain+1, Ndoses, Nplates))  
  ##attention : 1 : nb de lignes
  ## 2 : nb de colonne , 3 : nb de matrices
  
  chain.tau = matrix(NA, nchain+1, 1)
  colnames(chain.tau) = "tau"
  
  #init 
  chain.abg[1,] = init.abg
  chain.lambda[1,,1] = init.lambda[,1]
  chain.lambda[1,,2] = init.lambda[,2]
  chain.lambda[1,,3] = init.lambda[,3]
  #init.lambda doit 
  #etre une matrice de taille Ndoses*Nplates
  chain.tau[1] = init.tau
  
  tau_0 = 10^-6
  
  acc_rates.abg = rep(0, 3)
  acc_rates.lambda = matrix(0, nrow = 6, ncol = 3)

  for (iter in 1:nchain){
    current.abg = chain.abg[iter, ]
    current.lambda = chain.lambda[iter,,]
    current.tau = chain.tau[iter, ]
    
    mu = mu_fn(current.abg[1], current.abg[2], current.abg[3], x, current.lambda)
    
    # MAJ alpha beta gamma
    for(p in 1:3){
      prop = current.abg
      prop[p] = rnorm(1, current.abg[p], prop.sd.abg[p])
      
      mu_prop = mu_fn(prop[1], prop[2], prop[3], x, current.lambda)
      
      top = dnorm(prop[p], 0, sqrt(1/tau_0), log=TRUE)+
        sum(dpois( y, mu_prop, log=TRUE))
      
      bottom = dnorm(current.abg[p], 0, sqrt(1/tau_0), log=TRUE) + 
        sum( dpois( y, mu, log=TRUE))
      
      acc_prob.abg <- exp(top - bottom)
      
      if (runif(1) < acc_prob.abg){
        current.abg[p] <- prop[p]
        mu = mu_prop
        acc_rates.abg[p] = acc_rates.abg[p] + 1
      }
    }
    
    #MAJ tau
    update_shape = 10^-3 + Ndoses*Nplates/2 
    update_rate = sum(colSums(current.lambda^2))/2 + 10^-3
    current.tau = rgamma(1, shape = update_shape, rate = update_rate)
    
    #MAJ lambda
    for (i in 1:Ndoses){
      for(j in 1:Nplates){
        prop = current.lambda
        
        prop[i,j] = rnorm(1, current.lambda[i,j], prop.sd.lambda)
        
        mu_prop =  mu_fn(current.abg[1], current.abg[2], 
                         current.abg[3],x[i] ,prop[i,j])
        
        top = dnorm(prop[i,j],0,sqrt(1/current.tau),log=TRUE)+
          dpois(y[i,j],mu_prop,log=TRUE)
        bottom = dnorm(current.lambda[i,j],0,sqrt(1/current.tau),log=TRUE)+
          dpois(y[i,j],mu[i,j],log=TRUE)
        
        acc.prop.lambda = exp(top - bottom)
        
        if(runif(1) < acc.prop.lambda){
          current.lambda = prop
          mu[i,j] = mu_prop
          acc_rates.lambda[i, j] = acc_rates.lambda[i, j] + 1
        }
      }
    }
    
    chain.tau[iter+1, ] = current.tau
    chain.abg[iter+1, ] = current.abg
    chain.lambda[iter+1,,] = current.lambda
  }

  
  return(list(chain.tau = chain.tau, chain.abg = chain.abg,
              chain.lambda = chain.lambda, acc_rates.abg = acc_rates.abg/nchain,
              acc_rates.lambda = acc_rates.lambda/nchain))
  
}


nchain = 11000
init.abg = c(2, 0, 0.001) # On ne touche plus a ça quoi qu'il arrive !!!!!!!!!
init.lambda = matrix(0, 6, 3)
init.tau = 4

prop.sd.abg = c(0.09, 0.01, 0.0001)
prop.sd.lambda = 0.2 # ça c'est pas mal on laisse !

test = salm(nchain, init.abg, init.lambda, init.tau, prop.sd.abg, prop.sd.lambda)
par(mfrow = c(1,3))
for (j in 1:3){
  plot(test$chain.abg[,j], type = 'l')
}

test$acc_rates.abg

plot(test$chain.lambda[,2,3], type = 'l')

test$acc_rates.lambda

#---------------------------Stats---------------------

## Burning
out <- list(chain.tau = test$chain.tau[-(1:1000),],
            chain.lambda = test$chain.lambda[-(1:1000),,],
            chain.abg = test$chain.abg[-(1:1000),])

## Chaines
# alpha, beta et gamma
par(mfrow = c(1, 3))
for (i in (1:3)){
  plot(out$chain.abg[,i], type = 'l', main = "")
}

# lambda (4)
par(mfrow = c(2, 2))
for (i in (1:2)){
  for (j in (1:2)){
    plot(out$chain.lambda[,i,j], type = 'l', main = "") 
  }
}

## densité
# alpha, beta et gamma
par(mfrow = c(1, 3))
for (i in (1:3)){
  plot(density(out$chain.abg[,i]), main = "")
}

# lambda
par(mfrow = c(2, 2))
for (i in (1:2)){
  for (j in 1:2){
    plot(density(out$chain.lambda[,i,j]), main = "") 
  }
}

# tau
par(mfrow = c(1, 1))
plot(density(out$chain.tau), main = "")

## Adéquation du modèle
# simulation
simu <- matrix(NA, nrow(out$chain.abg), 6)
proba <- matrix(NA, 6, 3)
for (t in 1:nrow(out$chain.abg)){
  estim.abg <- out$chain.abg[t,]
  estim.lambda <- out$chain.lambda[t,,]
  proba <- mu_fn(estim.abg[1], estim.abg[2], estim.abg[3], x, estim.lambda)
  
  simu[t,] <- rpois(6, proba)
}

par(mfrow = c(2, 3))
for (j in 1:6){
  plot(table(simu[,j]) / nrow(simu))
  abline(v = median(y[j,]), col = 'red')
}
# A l'air pas mal

## Evolution
dosage <- log(seq(0.01, 1000, length = 1000))
proba_pl1 <- matrix(NA, nrow(out$chain.abg), length(dosage))
proba_pl2 <- matrix(NA, nrow(out$chain.abg), length(dosage))
proba_pl3 <- matrix(NA, nrow(out$chain.abg), length(dosage))

for (t in 1:nrow(out$chain.abg)){
  estim.abg <- out$chain.abg[t,]
  estim.lambda_pl1 <- out$chain.lambda[t,,1]
  estim.lambda_pl2 <- out$chain.lambda[t,,2]
  estim.lambda_pl3 <- out$chain.lambda[t,,3]

  proba_pl1[t,] <- mu_fn(estim.abg[1], estim.abg[2], estim.abg[3], dosage,
                     estim.lambda_pl1)
  proba_pl2[t,] <- mu_fn(estim.abg[1], estim.abg[2], estim.abg[3], dosage,
                     estim.lambda_pl2)
  proba_pl3[t,] <- mu_fn(estim.abg[1], estim.abg[2], estim.abg[3], dosage,
                     estim.lambda_pl3)
}

# Graphiques
par(mfrow = c(1, 3))

# Plaque 1
plot(dosage, apply(proba_pl1, 2, median), col = c(1:6),
     xlab = "Dosage", ylab = "Proba de mutation")
lines(dosage, apply(proba_pl1, 2, quantile, prob = 0.025), lty = 2,
      col = "grey")
lines(dosage, apply(proba_pl1, 2, quantile, prob = 0.975), lty = 2,
      col = "grey")

# Plaque 2
plot(dosage, apply(proba_pl2, 2, median), type = "l", col = "orange",
     xlab = "Dosage", ylab = "Proba de mutation")
lines(dosage, apply(proba_pl2, 2, quantile, prob = 0.025), lty = 2,
      col = "grey")
lines(dosage, apply(proba_pl2, 2, quantile, prob = 0.975), lty = 2,
      col = "grey")

# Plaque 3
plot(dosage, apply(proba_pl3, 2, median), type = "l", col = "orange",
     xlab = "Dosage", ylab = "Proba de mutation")
lines(dosage, apply(proba_pl3, 2, quantile, prob = 0.025), lty = 2,
      col = "grey")
lines(dosage, apply(proba_pl3, 2, quantile, prob = 0.975), lty = 2,
      col = "grey")
# Très moche




