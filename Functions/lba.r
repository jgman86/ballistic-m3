####################################################
# Code to go with the paper by Brown & Heathcote   #
# on the linear ballistic accumulator.             #
####################################################

# Function to calculate the threshold crossing time 
# CDF for one node (node #1). Arguments are:
#  - t = time vector at which CDF is evaluated.
#  - Parameters of the LBA model (A,b,v,s). See paper.
# We call parameter "s" as "sdv" to avoid confusion
# with some "integrate" function parameters.
# The drift rates (v) should be a vector of length
# equal to the number of response accumulators (choices).
# Many functions do not include an argument for the
# non-decision time parameter (t0). To allow for this,
# just pass time vectors (t) in as (t-t0), or add
# t0 to output means.
fptcdf=function(t,A,b,v,sdv) {
  zs=t*sdv ; zu=t*v ; chiminuszu=b-zu ; xx=chiminuszu-A
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  1+(tmp1+tmp2)/A
}

# Function to calculate the threshold crossing time 
# PDF for one node (node #1). Arguments are as for
# fptcdf above.
fptpdf=function(t,A,b,v,sdv) {
  zs=t*sdv ; zu=t*v ; chiminuszu=b-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-A)/zs
  (v*(pnorm(chizu)-pnorm(chizumax)) + 
    sdv*(dnorm(chizumax)-dnorm(chizu)))/A
}

# Function to calculate the CDF for overall
# RT, irrespective of response. Arguments are as
# for the fptcdf function above
allrtCDF=function(t,A,b,v,sdv) {
  # Generates CDF for all RTs irrespective of response.
  N=length(v) # Number of responses.
  tmp=array(dim=c(length(t),N))
  for (i in 1:N) tmp[,i]=fptcdf(t=t,A=A,b=b,v=v[i],sdv=sdv)
  1-apply(1-tmp,1,prod)
}

# Generates defective PDF for responses on node #1. 
# Arguments are as for the fptcdf function above.
# Takes a different calculation approach for binary
# (N=2) choices than multiple (N>2) choices, for 
# computational efficiency.
n1PDF=function(t,A,b,v,sdv) {
  N=length(v) # Number of responses.
  if (N>2) {
    tmp=array(dim=c(length(t),N-1))
    for (i in 2:N) tmp[,i-1]=fptcdf(t=t,A=A,b=b,v=v[i],sdv=sdv)
    G=apply(1-tmp,1,prod)    
  } else {
    G=1-fptcdf(t=t,A=A,b=b,v=v[2],sdv=sdv)
  }
  G*fptpdf(t=t,A=A,b=b,v=v[1],sdv=sdv)
}

# Generates defective CDF for responses on node #1. 
# Arguments are as for the fptcdf function above.
n1CDF=function(t,A,b,v,sdv) {
  # Generates defective CDF for responses on node #1. 
  outs=numeric(length(t)) ; bounds=c(0,t)
  for (i in 1:length(t)) {
    tmp="error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {outs[i]=0;break}
      tmp=try(integrate(f=n1PDF,lower=bounds[i],upper=bounds[i+1],
        A=A,b=b,v=v,sdv=sdv)$value,silent=T)
      if (is.numeric(tmp)) {outs[i]=tmp;break}
      # If we get to here, the integration call failed somehow.
      # Try smart lower bound.
      if (bounds[i]<=0) {
	bounds[i]=max(c((b-0.98*A)/(max(mean(v),v[1])+2*sdv),0))
	next
      }
      # Try smart upper bound.
      if (bounds[i+1]==Inf) {
	bounds[i+1]=0.02*b/(mean(v)-2*sdv)
	next
      }
      stop("Error in n1CDF that I could not catch.")
    }
  }
  cumsum(outs)
}

# Generates mean RT for responses on node #1. 
# Arguments are as for the fptcdf function above.
# Returns a list with mean RT for responses on
# node #1, and the probability of those responses.
# Does not inlcude t0.
n1mean=function(A,b,v,sdv) {
   pc=n1CDF(Inf,A,b,v,sdv)
   fn=function(t,A,b,v,sdv,pc) t*n1PDF(t,A,b,v,sdv)/pc
   tmp=integrate(f=fn,lower=0,upper=100*b,A=A,b=b,pc=pc,
     v=v,sdv=sdv)$value
   list(mean=tmp,p=pc)
}

# Calculates CDF for activation (z) at time (t) 
# for node #1. Arguments are as for the fptcdf
# function above. Does not inlcude t0.
actCDF=function(z,t,A,b,v,sdv) {
  zat=(z-A)/t ; zt=z/t ; s2=2*(sdv^2)
  exp1=exp(-((zt-drift)^2)/s2)
  exp2=exp(-((zat-drift)^2)/s2)
  tmp1=t*sdv*(exp1-exp2)/sqrt(2*pi)
  tmp2=pnorm(zat,mean=drift,sd=sdv)
  tmp3=pnorm(zt,mean=drift,sd=sdv)
  (tmp1+(A-z+drift*t)*tmp2+(z-drift*t)*tmp3)/A
}

# Calculates CDF for activation (z) at time (t) 
# for node #1. Arguments are as for the fptcdf
# function above. Does not inlcude t0.
actPDF=function(z,t,A,b,v,sdv) {
  tmp1=pnorm((z-A)/t,mean=v,sd=sdv)
  tmp2=pnorm(z/t,mean=v,sd=sdv)
  (-tmp1+tmp2)/A
}

# Calculates the probability of a response
# on node #1 at time t in a deadline experiment
# where there are no response thresholds.
# Arguments are as for the ftpcdf function 
# above. Note that t does not include t0.
deadlineaccuracy=function(t,A,b,v,sdv) {
  N=length(v)
  tmpf=function(x,t,A,b,v,sdv) {
    if (N>2) {
      tmp=array(dim=c(length(x),N-1))
      for (i in 2:N) tmp[,i-1]=actCDF(x,t,A,b,v[i],sdv)
      G=apply(tmp,1,prod)
    } else {
      G=actCDF(x,t,A,b,v[2],sdv)
    }
    G*actPDF(x,t,A,b,v[1],sdv)
  }
  outs=numeric(length(t))
  for (i in 1:length(t)) {
    if (t[i]<=0) {
      outs[i]=.5
    } else {
      outs[i]=integrate(f=tmpf,lower=-Inf,upper=Inf,t=t[i],
        A=A,v=v,sdv=sdv)$value
    }
  }
  outs
}
  
# Calculates the predicted quantile RT for probabilities
# given by "qps", on node #1. Arguments are as for the
# ftpcdf function above. Output does not include t0.
n1q=function(A,b,v,sdv,qps=seq(.1,.9,.2)) {
  out=list(predq=numeric(length(qps)))
  out$p=n1CDF(t=Inf,A=A,b=b,v=v,sdv=sdv)
  tmpf=function(t,v,sdv,A,b,p) n1CDF(t=t,A=A,b=b,v=v,sdv=sdv)-p
  for (j in 1:length(qps)) out$predq[j]=uniroot(f=tmpf,lower=0,
    upper=A*20,A=A,b=b,v=v,sdv=sdv,p=qps[j]*out$p)$root
  out
}



# Function for computing probability distribution over possible responses
LBAmulti <- function(t, a, b, drift, sdv) {
  Presponse <- matrix(0, 1, length(drift)) 
  for (r in 1:length(drift)) {
    v <- c(drift[r], drift[setdiff(1:length(drift), r)])   # puts the drift rate of the response currently considered in front
    Presponse[r] <- n1CDF(t=t, A=a, b=b, v=v, sdv=sdv)
  }
  return(Presponse)
}

LBA <- function(A, ch, a=400, b=420, s=0.1) {
  Pr <- matrix(NA, dim(A)[1], sum(ch))
  nsubj <- dim(A)[1]
  for (subj in 1:nsubj) {
    V <- c(A[subj,1], rep(A[subj,2], ch[2]), rep(A[subj,3], ch[3]), rep(A[subj,4], ch[4]), rep(A[subj,5], ch[5]))  
    V <- V/sum(V)  # normalize drift rates, as recommended by Brown & Heathcote
    Pr[subj,] <- LBAmulti(t=Inf, a=a, b=b, drift=V, sdv=s)
  }
  uIdx <- cumsum(ch)
  P <- matrix(Pr[, 1:uIdx[1]], nrow=nsubj)
  for (j in 1:(length(ch)-1)) P <- cbind( P, rowSums(matrix(Pr[, (uIdx[j]+1):uIdx[j+1]], nrow=nsubj)) )
  return(P)
}


LBA_v2 <- function(A, ch, a=400, b=420, s=0.1) {
  Pr <- matrix(NA, dim(A)[1], sum(ch))
  V <- matrix(NA, dim(A)[1], sum(ch))
  nsubj <- dim(A)[1]
  for (subj in 1:nsubj) {
    V[subj,] <- c(A[subj,1], rep(A[subj,2], ch[2]), rep(A[subj,3], ch[3]), rep(A[subj,4], ch[4]), rep(A[subj,5], ch[5]))  
    V[subj,] <- V[subj,]/sum(V[subj,])  # normalize drift rates, as recommended by Brown & Heathcote
    Pr[subj,] <- LBAmulti(t=Inf, a=a, b=b, drift=V[subj,], sdv=s)
  }
  uIdx <- cumsum(ch)
  P <- matrix(Pr[, 1:uIdx[1]], nrow=nsubj)
  
  V_mat <- matrix(V[, 1:uIdx[1]], nrow=nsubj)
  
  
  for (j in 1:(length(ch)-1)) {
    
    P <- cbind( P, rowSums(matrix(Pr[, (uIdx[j]+1):uIdx[j+1]], nrow=nsubj)))
    V_mat <- cbind(V_mat, rowMeans(matrix(V[, (uIdx[j]+1):uIdx[j+1]], nrow=nsubj)))
    
  }
  
  out <- list(P,V_mat)
  return(out)
}


# Function for computing probability distribution over possible responses
LBAmulti <- function(t, a, b, drift, sdv) {
  Presponse <- matrix(0, 1, length(drift)) 
  for (r in 1:length(drift)) {
    v <- c(drift[r], drift[setdiff(1:length(drift), r)])   # puts the drift rate of the response currently considered in front
    Presponse[r] <- n1CDF(t=t, A=a, b=b, v=v, sdv=sdv)
  }
  return(Presponse)
}

rlba=function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
  n.with.extras=ceiling(n*(1+3*prod(pnorm(-vs))))
  drifts=matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE)
  if (truncdrifts) {
    repeat {
      drifts=rbind(drifts,matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE))
      tmp=apply(drifts,1,function(x) any(x>0))
      drifts=drifts[tmp,]
      if (nrow(drifts)>=n) break
    }
  }
  drifts=drifts[1:n,]
  drifts[drifts<0]=0
  starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
  ttf=t((b-t(starts)))/drifts
  rt=apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
  resp=apply(ttf,1,which.min)
  list(rt=rt,resp=resp)
}
