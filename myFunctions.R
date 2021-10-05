###
#Função para amostrar lambdas#
##

#Funçao para amostrar valores de lambda das minhas simulações
sampler <- function(data1,n,nParques){            #data1=matrix, n de samples
  bigSample <- matrix(nrow=n,ncol=nParques)
  for(i in 1:nParques) {
    bigSample[,i] = sample(data1[,i],n, replace= TRUE)
  }            
  return(bigSample)
}

#Função para gerar o S através de Poiss(Lambda)
s.gen <-function(sample){
  s.sample <- matrix(0,nrow=nrow(sample),ncol=ncol(sample))
  for(i in 1: ncol(s.sample)){
    for(j in 1: nrow(s.sample)){
      s.sample[j,i] <- rpois(1,sample[j,i])
    }
  }
  return (s.sample)
}

#Calculador de WAIC
waic <- function(stanfit){
  log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
              pointwise=pointwise, total=total, se=se))
}
  
colVars <- function(a){
  n <- dim(a)[[1]]
  c <- dim(a)[[2]]
  return(.colMeans(((a - matrix(.colMeans(a, n, c), 
                                nrow = n, ncol =c, byrow=TRUE))^2),n,c)* n/(n-1))
}

####################################################################################






