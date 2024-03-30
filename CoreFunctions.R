require("genlasso")

rtaCFR.EST<-function(ct, dt, F_mean = 15.43, F_shape = 2.03, maxsteps = 10000){
  if(length(ct)!=length(dt)){print("ERROR: The length of ct is not equal to that of dt.");break}
  N<-length(ct)
  fs  <-diff(pgamma(q=c(0:N), shape = F_shape,rate = F_shape/F_mean))
  fmat<-matrix(0,nrow = N,ncol = N)
  for(i in c(1:N)){fmat[i,c(1:i)]<-rev(fs[1:i])}
  Xmat    <-fmat%*%diag(ct) #Xi
  D_mat   <-cbind(diag(1,length(dt)-1),0)-cbind(0,diag(1,length(dt)-1))
  a0      <-fusedlasso(y = dt,X = Xmat, D = D_mat, gamma = 0,maxsteps = maxsteps,minlam = 0)
  Sum_a0  <-as.data.frame(summary(a0),row.names = F)
  steps   <-nrow(Sum_a0)
  cd      <-which(as.numeric(apply((a0$beta>=0)*(a0$beta<=1),MARGIN = 2,prod))==1)
  cd_star <-which(Sum_a0$rss==min(Sum_a0$rss[cd]))
  p_hat   <-a0$beta[,cd_star]
  lam_star<-Sum_a0$lambda[cd_star]
  return(list(p_hat=p_hat, lambda_star=lam_star, steps=steps, p_matrix=a0$beta, sol_path=Sum_a0))
}

srtaCFR.SIM<-function(ct, ct_prop_mat, pt_mat, seed = NA, F_mean = 15.43, F_shape = 2.03){
if(is.matrix(ct_prop_mat)==FALSE){print("ERROR: ct_prop_mat is not a matrix.");break}
if(is.matrix(pt_mat)==FALSE){print("ERROR: pt_mat is not a matrix.");break}
if(length(ct)!=nrow(ct_prop_mat)){print("ERROR: The length of ct should match with the number of row of ct_prop_mat.");break}
if(length(ct)!=nrow(pt_mat)){print("ERROR: The length of ct should match with the number of row of pt_mat.");break}
if(ncol(ct_prop_mat)!=ncol(pt_mat)){print("ERROR: The number of column of ct_prop_mat should match with the number of column of pt_mat.");break}
  
N     <-length(ct)
P     <-ncol(ct_prop_mat)
ct_mat<-sapply(c(1:N),function(o){rmultinom(n = 1, size = ct[o],prob = ct_prop_mat[o,])})
fs    <-diff(pgamma(q=c(0:N), shape = F_shape,rate = F_shape/F_mean))
dt_age<-matrix(0,nrow = P,ncol = N)
for(j in (1:P)){
  dt_age[j,]<-rbinom(n = N,size = ct_mat[j,], prob = pt_mat[,j])
}
dt_mat <-c()
for(k in c(1:P)){dt_mat<-rbind(dt_mat, sapply(c(1:N),function(o){sum(dt_age[k,1:o]*rev(fs[1:o]))}))}
dt_sum<-apply(dt_mat,2,sum)

return(list(ct_mat=as.matrix(t(ct_mat)),dt_mat=as.matrix(t(dt_mat))))
}

srtaCFR<-function(ct_mat, dt_mat, q_mat=NA, F_mean = 15.43, F_shape = 2.03, maxsteps = 10000){
  if(is.matrix(ct_mat)==FALSE){print("ERROR: ct_mat is not a matrix.");break}
  if(is.matrix(dt_mat)==FALSE){print("ERROR: dt_mat is not a matrix.");break}
  J<-ncol(ct_mat)
  if(J<2){print("ERROR: At least two categories (J>=2) are needed to stratify the data.");break}
  
  N<-nrow(ct_mat)
  if(is.matrix(q_mat)==FALSE){
    print("Using the default values (even distribution) for q_mat.")
    q_mat<-matrix(1/J,nrow = N,ncol = J)
  }else{
    print("Using the user-defined values for q_mat.")
    if(floor(sum(q_mat))!=N){print("ERROR: Probabilities in a row of q_mat do not sum to 1.");break}
  }
  
  if(!((J==ncol(dt_mat))&(J==ncol(q_mat)))){
    print("ERROR: The #columns of ct_mat, dt_mat, q_mat do not match each other.");break}
  if(!((N==nrow(dt_mat))&(N==nrow(q_mat)))){
    print("ERROR: The #rows of ct_mat, dt_mat, q_mat do not match each other.");break}
  
  fs  <-diff(pgamma(q=c(0:N), shape = F_shape,rate = F_shape/F_mean))
  fmat<-matrix(0,nrow = N,ncol = N)
  for(i in c(1:N)){fmat[i,c(1:i)]<-rev(fs[1:i])}
  
  lam_op<-lam_star<-steps<-numeric(J)
  p_hat_mat<-matrix(0,nrow=N, ncol=J)
  p_hat_tensor<-sol_path_tensor<-vector(mode = "list", length = J)
  for(j in c(1:J)){
    Xmat_j  <-fmat%*%diag(ct_mat[,j]) #Xi
    D_mat_j <-cbind(diag(1,length(dt_mat[,j])-1),0)-cbind(0,diag(1,length(dt_mat[,j])-1))
    a0      <-fusedlasso(y = dt_mat[,j],X = Xmat_j, D = D_mat_j, gamma = 0,maxsteps = maxsteps,minlam = 0)
    Sum_a0  <-as.data.frame(summary(a0),row.names = F)
    steps[j]<-nrow(Sum_a0)
    cd      <-which(as.numeric(apply((a0$beta>=0)*(a0$beta<=1),MARGIN = 2,prod))==1)
    cd_star <-which(Sum_a0$rss==min(Sum_a0$rss[cd]))
    p_hat_mat[,j]<-a0$beta[,cd_star]
    lam_star[j]<-Sum_a0$lambda[cd_star]
    p_hat_tensor[[j]]<-a0$beta
    sol_path_tensor[[j]]<-Sum_a0
  }
  
  p_hat_std<-apply(p_hat_mat*q_mat,1,sum)
  
  return(list(p_hat_std=p_hat_std, p_gp_spec=p_hat_mat, lambda_star_vec=lam_star,
              steps_vec=steps,p_tensor=p_hat_tensor, sol_path_tensor=sol_path_tensor))
}

# #Example
# fit1<-srtaCFR(ct_mat = t(ct_mat), dt_mat = t(dt_mat), lambda_bound = 5000)
# plot(fit1$p_gp_spec[,1], ylab="group-specific fatality rates", xlab="Time", ylim=c(0,0.10), type="l", lwd=2)
# lines(fit1$p_gp_spec[,2], lty=2, type="l", lwd=2, col="red")
# lines(fit1$p_gp_spec[,3], lty=3, type="l", lwd=2, col="blue") 
# plot(fit1$p_hat_std, ylab="srtaCFR", xlab="Time",ylim=c(0,0.06), type="l", lwd=2)

#Example 
#Scenario III
Data <- srtaCFR.SIM(seed = 1,
          ct=10000-50*abs(100-c(1:200)),
          ct_prop_mat=cbind(seq(0.2,0.6,length.out=200),0.2,seq(0.6,0.2,length.out=200)),
          pt_mat=cbind(rep(0.01,200),rep(0.02,200),rep(0.06,200))*replicate(3,exp(c(1:200)*0.004)) 
        )

rt      <-rtaCFR.EST(ct = rowSums(Data$ct_mat), dt = rowSums(Data$dt_mat))
#srt_fit <-srtaCFR(ct_mat = Data$ct_mat, dt_mat = Data$dt_mat, q_mat = matrix(1/3,nrow=nrow(Data$ct_mat), ncol=ncol(Data$ct_mat)))
srt_fit <-srtaCFR(ct_mat = Data$ct_mat, dt_mat = Data$dt_mat)

# plot(srt_fit$p_gp_spec[,1], ylab="Age-group-specific fatality rates", xlab="Time", type="l", #xaxt="n", 
#      lwd=2, col="darkgreen", xlim=c(0,200), ylim=c(0,0.15))
# lines(srt_fit$p_gp_spec[,2], lty=2, type="l", lwd=2, col="red")
# lines(srt_fit$p_gp_spec[,3], lty=3, type="l", lwd=2, col="blue")
# legend("topleft",col=c("darkgreen","red","blue"),lty=c(1,2,3),
#        lwd = c(1.5,1.5,1.5),bty = "n",x.intersp=0.5,text.width=35,
#        legend = c(expression(paste("rtaCFR"^"(1)")),
#                   expression(paste("rtaCFR"^"(2)")),
#                   expression(paste("rtaCFR"^"(3)"))),ncol = 3,cex=1)

plot(srt_fit$p_hat_std, lty=1, type="l", lwd=2, col="darkgreen", xlab="Time", ylab="Fatality rates",
     xlim=c(0,200), ylim=c(0,0.1))
lines(rt$p_hat, lty=2, lwd=2, col="red")
legend("topleft",col=c("darkgreen","red"),lty=c(1,2),
       lwd = c(1.5,1.5),bty = "n",x.intersp=0.5,text.width=40,
       legend = c(expression(paste("srtaCFR")),
                  expression(paste("rtaCFR (non-stratified)"))),ncol = 3,cex=1)



