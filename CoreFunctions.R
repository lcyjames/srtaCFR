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
set.seed(seed)
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
