
nullspace_sol = function(P,Q,b){
  if(dim(Q)[1]==1){
    sv <- svd(Q,nv=dim(Q)[2])
    D = matrix(c(1/sv$d,rep(0,dim(Q)[2]-1)),ncol=1,nrow=dim(Q)[2])
    xp <- sv$v %*% D %*% t(sv$u) %*% b
  }else{
    sv <- svd(Q)
    D = diag(1/sv$d)
    xp <- sv$v %*% D %*% t(sv$u) %*% b
  }
  res = qr(Q)
  rankQ = res$rank
  if(rankQ==dim(Q)[2]){
    c_5 = xp # Null space
  }else{
    Ns = sv$v[,(rankQ+1):dim(Q)[2]]
    za = -solve(t(Ns)%*%t(P)%*%P%*%Ns,t(Ns)%*%t(P)%*%P%*%xp,tol=0)
    c_5 = xp + Ns%*%za # Null space
  }
  return(c_5)
}

