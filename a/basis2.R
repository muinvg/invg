# phi0       : polynomial       / 1,x*x
# phiA       : polynomial       / 1,x,y, ...
# phiB       : polynomial-1     / x,y, ...
# phiC       : bernstein2       / x^i*(1-x)^j*y^k*(1-y)^l
# phiD       : legendre2-1 normal / 
# bnst[[N]]  : bernsteinN
# ebnst[[N]] : bernsteinN[-1,1]
# lgdr[[N]]  : legendreN
# nlgdr[[N]] : legendreN normal
# nmlgdr[[N]]: legendreN-1 normal

phi0=list(function(x,y){  x*0 + 1},function(x,y){x*x})
phi0x=list(function(x,y){  x*0 + 0},function(x,y){2*x})
phi0y=list( function(x,y){  x*0 + 0},function(x,y){  x*0 + 0})
phi0d = list(phi0x,phi0y)
phi0_name = c("1","x*x")


phiA=list(
  function(x,y){  x*0 + 1},
  function(x,y){  x},function(x,y){  y},
  function(x,y){x*x},function(x,y){x*y},function(x,y){y*y},
  function(x,y){x*x*x},function(x,y){x*x*y},function(x,y){x*y*y},function(x,y){y*y*y},
  function(x,y){x*x*x*x},function(x,y){x*x*x*y},function(x,y){x*x*y*y},function(x,y){x*y*y*y},function(x,y){y*y*y*y},
  function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y},function(x,y){x*x*x*y*y},
  function(x,y){x*x*y*y*y},function(x,y){x*y*y*y*y},function(x,y){y*y*y*y*y},
  function(x,y){x*x*x*x*x*x},function(x,y){x*x*x*x*x*y},function(x,y){x*x*x*x*y*y},
  function(x,y){x*x*x*y*y*y},function(x,y){x*x*y*y*y*y},function(x,y){x*y*y*y*y*y},function(x,y){y*y*y*y*y*y})

phiAx=list(
  function(x,y){  x*0 + 0},
  function(x,y){  x*0 + 1},function(x,y){  y*0},
  function(x,y){2*x},function(x,y){y},function(x,y){0*y},
  function(x,y){3*x*x},function(x,y){2*x*y},function(x,y){y*y},function(x,y){0*y},
  function(x,y){4*x*x*x},function(x,y){3*x*x*y},function(x,y){2*x*y*y},function(x,y){y*y*y},function(x,y){0*y},
  function(x,y){5*x*x*x*x},function(x,y){4*x*x*x*y},function(x,y){3*x*x*y*y},
  function(x,y){2*x*y*y*y},function(x,y){y*y*y*y},function(x,y){0*y},
  function(x,y){6*x*x*x*x*x},function(x,y){5*x*x*x*x*y},function(x,y){4*x*x*x*y*y},
  function(x,y){3*x*x*y*y*y},function(x,y){2*x*y*y*y*y},function(x,y){y*y*y*y*y},function(x,y){0*y})

phiAy=list(
  function(x,y){  x*0 + 0},
  function(x,y){  x*0},function(x,y){  y*0 + 1},
  function(x,y){0*x},function(x,y){x},function(x,y){2*y},
  function(x,y){0*x},function(x,y){x*x},function(x,y){x*y*2},function(x,y){3*y*y},
  function(x,y){0*x},function(x,y){x*x*x},function(x,y){x*x*y*2},function(x,y){x*y*y*3},function(x,y){y*y*y*4},
  function(x,y){x*0},function(x,y){x*x*x*x},function(x,y){x*x*x*y*2},
  function(x,y){x*x*y*y*3},function(x,y){x*y*y*y*4},function(x,y){y*y*y*y*5},
  function(x,y){x*0},function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y*2},
  function(x,y){x*x*x*y*y*3},function(x,y){x*x*y*y*y*4},function(x,y){x*y*y*y*y*5},function(x,y){y*y*y*y*y*6})

phiAd = list(phiAx,phiAy)

phiA_name = c("1","x","y","x*x","x*y","y*y",
             "x*x*x","x*x*y","x*y*y","y*y*y","x*x*x*x","x*x*x*y","x*x*y*y","x*y*y*y","y*y*y*y",
             "x*x*x*x*x","x*x*x*x*y","x*x*x*y*y","x*x*y*y*y","x*y*y*y*y","y*y*y*y*y",
             "x*x*x*x*x*x","x*x*x*x*x*y","x*x*x*x*y*y","x*x*x*y*y*y","x*x*y*y*y*y","x*y*y*y*y*y","y*y*y*y*y*y")



phiB=list(
  function(x,y){  x},function(x,y){  y},
  function(x,y){x*x},function(x,y){x*y},function(x,y){y*y},
  function(x,y){x*x*x},function(x,y){x*x*y},function(x,y){x*y*y},function(x,y){y*y*y},
  function(x,y){x*x*x*x},function(x,y){x*x*x*y},function(x,y){x*x*y*y},function(x,y){x*y*y*y},function(x,y){y*y*y*y},
  function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y},function(x,y){x*x*x*y*y},
  function(x,y){x*x*y*y*y},function(x,y){x*y*y*y*y},function(x,y){y*y*y*y*y},
  function(x,y){x*x*x*x*x*x},function(x,y){x*x*x*x*x*y},function(x,y){x*x*x*x*y*y},
  function(x,y){x*x*x*y*y*y},function(x,y){x*x*y*y*y*y},function(x,y){x*y*y*y*y*y},function(x,y){y*y*y*y*y*y})

phiBx=list(
  function(x,y){  x*0 + 1},function(x,y){  y*0},
  function(x,y){2*x},function(x,y){y},function(x,y){0*y},
  function(x,y){3*x*x},function(x,y){2*x*y},function(x,y){y*y},function(x,y){0*y},
  function(x,y){4*x*x*x},function(x,y){3*x*x*y},function(x,y){2*x*y*y},function(x,y){y*y*y},function(x,y){0*y},
  function(x,y){5*x*x*x*x},function(x,y){4*x*x*x*y},function(x,y){3*x*x*y*y},
  function(x,y){2*x*y*y*y},function(x,y){y*y*y*y},function(x,y){0*y},
  function(x,y){6*x*x*x*x*x},function(x,y){5*x*x*x*x*y},function(x,y){4*x*x*x*y*y},
  function(x,y){3*x*x*y*y*y},function(x,y){2*x*y*y*y*y},function(x,y){y*y*y*y*y},function(x,y){0*y})

phiBy=list(
  function(x,y){  x*0},function(x,y){  y*0 + 1},
  function(x,y){0*x},function(x,y){x},function(x,y){2*y},
  function(x,y){0*x},function(x,y){x*x},function(x,y){x*y*2},function(x,y){3*y*y},
  function(x,y){0*x},function(x,y){x*x*x},function(x,y){x*x*y*2},function(x,y){x*y*y*3},function(x,y){y*y*y*4},
  function(x,y){x*0},function(x,y){x*x*x*x},function(x,y){x*x*x*y*2},
  function(x,y){x*x*y*y*3},function(x,y){x*y*y*y*4},function(x,y){y*y*y*y*5},
  function(x,y){x*0},function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y*2},
  function(x,y){x*x*x*y*y*3},function(x,y){x*x*y*y*y*4},function(x,y){x*y*y*y*y*5},function(x,y){y*y*y*y*y*6})

phiBd = list(phiBx,phiBy)

phiB_name = c("x","y","x*x","x*y","y*y",
              "x*x*x","x*x*y","x*y*y","y*y*y","x*x*x*x","x*x*x*y","x*x*y*y","x*y*y*y","y*y*y*y",
              "x*x*x*x*x","x*x*x*x*y","x*x*x*y*y","x*x*y*y*y","x*y*y*y*y","y*y*y*y*y",
              "x*x*x*x*x*x","x*x*x*x*x*y","x*x*x*x*y*y","x*x*x*y*y*y","x*x*y*y*y*y","x*y*y*y*y*y","y*y*y*y*y*y")





phiC=list(
  function(x,y){(1-x)*(1-x)*(1-y)*(1-y)},
  function(x,y){(1-x)*(1-x)*2*y*(1-y)},
  function(x,y){(1-x)*(1-x)*y*y},
  function(x,y){2*x*(1-x)*(1-y)*(1-y)},
  function(x,y){2*x*(1-x)*2*y*(1-y)},
  function(x,y){2*x*(1-x)*y*y},
  function(x,y){x*x*(1-y)*(1-y)},
  function(x,y){x*x*2*y*(1-y)},
  function(x,y){x*x*y*y}
)

phiCx=list(
  function(x,y){-2*(1-x)*(1-y)*(1-y)},
  function(x,y){-2*(1-x)*2*y*(1-y)},
  function(x,y){-2*(1-x)*y*y},
  function(x,y){2*(1-2*x)*(1-y)*(1-y)},
  function(x,y){2*(1-2*x)*2*y*(1-y)},
  function(x,y){2*(1-2*x)*y*y},
  function(x,y){2*x*(1-y)*(1-y)},
  function(x,y){2*x*2*y*(1-y)},
  function(x,y){2*x*y*y}
)
phiCy=list(
  function(x,y){(1-x)*(1-x)*(1-y)*(-2)},
  function(x,y){(1-x)*(1-x)*2*(1-2*y)},
  function(x,y){(1-x)*(1-x)*y*2},
  function(x,y){2*x*(1-x)*(1-y)*(-2)},
  function(x,y){2*x*(1-x)*2*(1-2*y)},
  function(x,y){2*x*(1-x)*y*2},
  function(x,y){x*x*(1-y)*(-2)},
  function(x,y){x*x*2*(1-2*y)},
  function(x,y){x*x*y*2}
)
phiCd = list(phiCx,phiCy)

phiC_name = c("brs2-0-0","brs2-0-1","brs2-0-2",
              "brs2-1-0","brs2-1-1","brs2-1-2",
              "brs2-2-0","brs2-2-1","brs2-2-2")




bern1d <- function(k, n, u) { choose(n, k) * u^k * (1-u)^(n-k) }
dbern1d <- function(k, n, u) {
  if (n == 0) return(0*u)
  # derivative: n * [B_{k-1,n-1}(u) - B_{k,n-1}(u)] with boundary handling
  term1 <- if (k-1 >= 0) bern1d(k-1, n-1, u) else 0*u
  term2 <- if (k <= n-1) bern1d(k, n-1, u) else 0*u
  n * (term1 - term2)
}
# dbern1d <- function(k, n, u) {
#   if (n == 0) return(0*u)
#   if (k == 0) return(-n*(1-u)^(n-1))
#   if (k == n) return(n*u^(n-1))
#   return(choose(n, k) * (k - n*u) * u^(k-1) * (1-u)^(n-k-1) )
# }


#bnst [0,1]
maxdeg = 5
bnst = as.list(1:maxdeg)
bnstd = as.list(1:maxdeg)
bnst_name = as.list(1:maxdeg)
for(dd in 1:5){
  bnst[[dd]] = as.list(1:((dd+1)^2))
  bnstd[[dd]] = as.list(1:2)
  bnstd[[dd]][[1]] = as.list(1:((dd+1)^2))
  bnstd[[dd]][[2]] = as.list(1:((dd+1)^2))
  bnst_name[[dd]]= rep("",((dd+1)^2))
  
  ran = (1:((dd+1)^2)) -1 
  bnst[[dd]] = lapply(ran, function(kk,ee=dd) {
    ii = kk %/% (ee+1)
    jj = kk %% (ee+1)
    return( function(x, y) bern1d(ii,ee,x)*bern1d(jj,ee,y) )
  }) #クロージャをつかった書き方
  bnstd[[dd]][[1]] = lapply(ran, function(kk,ee=dd) {
    ii = kk %/% (ee+1)
    jj = kk %% (ee+1)
    return( function(x, y) dbern1d(ii,ee,x)*bern1d(jj,ee,y) )
  }) #クロージャをつかった書き方
  bnstd[[dd]][[2]] = lapply(ran, function(kk,ee=dd) {
    ii = kk %/% (ee+1)
    jj = kk %% (ee+1)
    return( function(x, y) bern1d(ii,ee,x)*dbern1d(jj,ee,y) )
  }) #クロージャをつかった書き方
  
  for(kk in ran){
    ii = kk %/% (dd+1)
    jj = kk %% (dd+1)
    bnst_name[[dd]][kk+1] = paste0("bnst",dd,"-",ii,"-",jj)
  }
}

#bnst [-1,1]
maxdeg = 5
ebnst = as.list(1:maxdeg)
ebnstd = as.list(1:maxdeg)
ebnst_name = as.list(1:maxdeg)
for(dd in 1:5){
  ebnst[[dd]] = as.list(1:((dd+1)^2))
  ebnstd[[dd]] = as.list(1:2)
  ebnstd[[dd]][[1]] = as.list(1:((dd+1)^2))
  ebnstd[[dd]][[2]] = as.list(1:((dd+1)^2))
  ebnst_name[[dd]]= rep("",((dd+1)^2))
  
  ran = (1:((dd+1)^2)) -1 
  ebnst[[dd]] = lapply(ran, function(kk,ee=dd) {
    ii = kk %/% (ee+1)
    jj = kk %% (ee+1)
    return( function(x, y) bern1d(ii,ee,(x+1)/2)*bern1d(jj,ee,(y+1)/2) )
  }) #クロージャをつかった書き方
  ebnstd[[dd]][[1]] = lapply(ran, function(kk,ee=dd) {
    ii = kk %/% (ee+1)
    jj = kk %% (ee+1)
    return( function(x, y) dbern1d(ii,ee,(x+1)/2)*bern1d(jj,ee,(y+1)/2) )
  }) #クロージャをつかった書き方
  ebnstd[[dd]][[2]] = lapply(ran, function(kk,ee=dd) {
    ii = kk %/% (ee+1)
    jj = kk %% (ee+1)
    return( function(x, y) bern1d(ii,ee,(x+1)/2)*dbern1d(jj,ee,(y+1)/2) )
  }) #クロージャをつかった書き方
  
  for(kk in ran){
    ii = kk %/% (dd+1)
    jj = kk %% (dd+1)
    ebnst_name[[dd]][kk+1] = paste0("ebnst",dd,"-",ii,"-",jj)
  }
}



# Legendre P_n(x) を計算する関数
legendreP1d <- function(n, x) {
  if (n == 0) return(rep(1, length(x)))
  if (n == 1) return(x)
  Pnm1 <- rep(1, length(x))  # P_0
  Pn   <- x
  
  if (n == 1) return(Pn)
  for (k in 1:(n-1)) {
    Pnp1 <- ((2*k + 1) * x * Pn - k * Pnm1) / (k + 1)
    Pnm1 <- Pn
    Pn   <- Pnp1
  }
  return(Pn)
}
# Legendreの導関数
dlegendreP1d <- function(n, x) {
  if (n == 0) return(rep(0, length(x)))
  if (n == 1) return(rep(1, length(x)))
  
  dPnm1 <- rep(0, length(x))  # P_0'
  dPn   <- rep(1, length(x))  # P_1'
  
  for (k in 1:(n-1)) {
    dPnp1 <- (2*k + 1) * legendreP1d(k, x) + dPnm1
    dPnm1 <- dPn
    dPn   <- dPnp1
  }
  return(dPn)
}

# Legendre（正規化） を計算する関数
nlegendreP1d <- function(n, x) {
  sc = sqrt(2/(2*n+1))
  if (n == 0) return(rep(1, length(x))/sc)
  if (n == 1) return(x/sc)
  Pnm1 <- rep(1, length(x))  # P_0
  Pn   <- x
  
  if (n == 1) return(Pn/sc)
  for (k in 1:(n-1)) {
    Pnp1 <- ((2*k + 1) * x * Pn - k * Pnm1) / (k + 1)
    Pnm1 <- Pn
    Pn   <- Pnp1
  }
  return(Pn/sc)
}
# Legendre導関数（正規化）
dnlegendreP1d <- function(n, x) {
  sc = sqrt(2/(2*n+1))
  if (n == 0) return(rep(0, length(x)))
  if (n == 1) return(rep(1, length(x))/sc)
  
  dPnm1 <- rep(0, length(x))  # P_0'
  dPn   <- rep(1, length(x))  # P_1'
  
  for (k in 1:(n-1)) {
    dPnp1 <- (2*k + 1) * legendreP1d(k, x) + dPnm1
    dPnm1 <- dPn
    dPn   <- dPnp1
  }
  return(dPn/sc)
}

enumerate_simple <- function(D) {
  out <- matrix(nrow=(D+1)^2, ncol=2)
  k <- 1
  for (m in 0:D) {
    for (i in 0:m) { out[k,] <- c(i, m); k <- k + 1 }   # (0,m),(1,m),...,(m,m)
    if (m >= 1) for (j in (m-1):0) { out[k,] <- c(m, j); k <- k + 1 } # (m,m-1),...,(m,0)
  }
  out
}


maxdeg = 5
lgdr = as.list(1:((maxdeg+1)^2))
lgdrd = as.list(1:2)
lgdrd[[1]] = as.list(1:((maxdeg+1)^2))
lgdrd[[2]] = as.list(1:((maxdeg+1)^2))
lgdr_name= rep("",((maxdeg+1)^2))

iijj = enumerate_simple(maxdeg)

ran = 1:((maxdeg+1)^2) 
lgdr = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) legendreP1d(ii,x)*legendreP1d(jj,y) )
}) #クロージャをつかった書き方
lgdrd[[1]] = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) dlegendreP1d(ii,x)*legendreP1d(jj,y) )
}) #クロージャをつかった書き方
lgdrd[[2]] = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) legendreP1d(ii,x)*dlegendreP1d(jj,y) )
}) #クロージャをつかった書き方

for(kk in ran){
  ii = iijj[kk,1]
  jj = iijj[kk,2]
  lgdr_name[kk] = paste0("lgdr","-",ii,"-",jj)
}



maxdeg = 5
nlgdr = as.list(1:((maxdeg+1)^2))
nlgdrd = as.list(1:2)
nlgdrd[[1]] = as.list(1:((maxdeg+1)^2))
nlgdrd[[2]] = as.list(1:((maxdeg+1)^2))
nlgdr_name= rep("",((maxdeg+1)^2))

iijj = enumerate_simple(maxdeg)

ran = 1:((maxdeg+1)^2) 
nlgdr = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) nlegendreP1d(ii,x)*nlegendreP1d(jj,y) )
}) #クロージャをつかった書き方
nlgdrd[[1]] = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) dnlegendreP1d(ii,x)*nlegendreP1d(jj,y) )
}) #クロージャをつかった書き方
nlgdrd[[2]] = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) nlegendreP1d(ii,x)*dnlegendreP1d(jj,y) )
}) #クロージャをつかった書き方

for(kk in ran){
  ii = iijj[kk,1]
  jj = iijj[kk,2]
  nlgdr_name[kk] = paste0("nlgdr","-",ii,"-",jj)
}



maxdeg = 5
nmlgdr = as.list(1:(-1+(maxdeg+1)^2))
nmlgdrd = as.list(1:2)
nmlgdrd[[1]] = as.list(1:(-1+(maxdeg+1)^2))
nmlgdrd[[2]] = as.list(1:(-1+(maxdeg+1)^2))
nmlgdr_name= rep("",(-1+(maxdeg+1)^2))

iijj = enumerate_simple(maxdeg)

ran = 2:((maxdeg+1)^2) #2から始める
nmlgdr = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) nlegendreP1d(ii,x)*nlegendreP1d(jj,y) )
}) #クロージャをつかった書き方
nmlgdrd[[1]] = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) dnlegendreP1d(ii,x)*nlegendreP1d(jj,y) )
}) #クロージャをつかった書き方
nmlgdrd[[2]] = lapply(ran, function(kk,ij=iijj) {
  ii = ij[kk,1]
  jj = ij[kk,2]
  return( function(x, y) nlegendreP1d(ii,x)*dnlegendreP1d(jj,y) )
}) #クロージャをつかった書き方

for(kk in ran){
  ii = iijj[kk,1]
  jj = iijj[kk,2]
  nmlgdr_name[kk-1] = paste0("nlgdr","-",ii,"-",jj)
}




nP0_f <- function(x) 1/sqrt(2/(2*0+1))
nP1_f <- function(x) x/sqrt(2/(2*1+1))
nP2_f <- function(x) 0.5 * (3*x^2 - 1)/sqrt(2/(2*2+1))
ndP0_f <- function(x) 0
ndP1_f <- function(x) 1/sqrt(2/(2*1+1))
ndP2_f <- function(x) 3*x/sqrt(2/(2*2+1))

phiD=list(
  function(x,y){  nP0_f(x)*nP1_f(y)},
  function(x,y){  nP1_f(x)*nP1_f(y)},
  function(x,y){  nP1_f(x)*nP0_f(y)},
  function(x,y){  nP0_f(x)*nP2_f(y)},
  function(x,y){  nP1_f(x)*nP2_f(y)},
  function(x,y){  nP2_f(x)*nP2_f(y)},
  function(x,y){  nP2_f(x)*nP1_f(y)},
  function(x,y){  nP2_f(x)*nP0_f(y)}
)
  
phiDx=list(
  function(x,y){  ndP0_f(x)*nP1_f(y)},
  function(x,y){  ndP1_f(x)*nP1_f(y)},
  function(x,y){  ndP1_f(x)*nP0_f(y)},
  function(x,y){  ndP0_f(x)*nP2_f(y)},
  function(x,y){  ndP1_f(x)*nP2_f(y)},
  function(x,y){  ndP2_f(x)*nP2_f(y)},
  function(x,y){  ndP2_f(x)*nP1_f(y)},
  function(x,y){  ndP2_f(x)*nP0_f(y)}
)

phiDy=list(
  function(x,y){  nP0_f(x)*ndP1_f(y)},
  function(x,y){  nP1_f(x)*ndP1_f(y)},
  function(x,y){  nP1_f(x)*ndP0_f(y)},
  function(x,y){  nP0_f(x)*ndP2_f(y)},
  function(x,y){  nP1_f(x)*ndP2_f(y)},
  function(x,y){  nP2_f(x)*ndP2_f(y)},
  function(x,y){  nP2_f(x)*ndP1_f(y)},
  function(x,y){  nP2_f(x)*ndP0_f(y)}
)

phiDd = list(phiDx,phiDy)

phiD_name = c("nlgdr-0-1",
              "nlgdr-1-1",
              "nlgdr-1-0",
              "nlgdr-0-2",
              "nlgdr-1-2",
              "nlgdr-2-2",
              "nlgdr-2-1",
              "nlgdr-2-0")




chk_lgdr=function(){
  P0_fun <- function(x) 1
  P1_fun <- function(x) x
  P2_fun <- function(x) 0.5 * (3*x^2 - 1)
  dP0_fun <- function(x) 0
  dP1_fun <- function(x) 1
  dP2_fun <- function(x) 3*x
  
  X = rbind(
    matrix(c(0,0,0,-1,-1,-1,1,1,1  ,0,-1,1,0,-1,1,0,-1,1),nrow=9,ncol=2),
    matrix(runif(n=9991*2)*2 -1,nrow=9991,ncol=2)
  )
  Tx = dim(X)[1]
  E = matrix(0,nrow=Tx,ncol=3*3)
  for(tt in 1:Tx){
    E[tt,] = c(
      lgdr[[1]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*P0_fun(X[tt,2]),
      lgdr[[2]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*P1_fun(X[tt,2]),
      lgdr[[3]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*P1_fun(X[tt,2]),
      lgdr[[4]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*P0_fun(X[tt,2]),
      lgdr[[5]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*P2_fun(X[tt,2]),
      lgdr[[6]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*P2_fun(X[tt,2]),
      lgdr[[7]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*P2_fun(X[tt,2]),
      lgdr[[8]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*P1_fun(X[tt,2]),
      lgdr[[9]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*P0_fun(X[tt,2])
    )
  }
  Ed = matrix(0,nrow=Tx,ncol=3*3*2)
  for(tt in 1:Tx){
    Ed[tt,] = c(
      lgdrd[[1]][[1]](X[tt,1],X[tt,2]) - dP0_fun(X[tt,1])*P0_fun(X[tt,2]),
      lgdrd[[1]][[2]](X[tt,1],X[tt,2]) - dP0_fun(X[tt,1])*P1_fun(X[tt,2]),
      lgdrd[[1]][[3]](X[tt,1],X[tt,2]) - dP1_fun(X[tt,1])*P1_fun(X[tt,2]),
      lgdrd[[1]][[4]](X[tt,1],X[tt,2]) - dP1_fun(X[tt,1])*P0_fun(X[tt,2]),
      lgdrd[[1]][[5]](X[tt,1],X[tt,2]) - dP0_fun(X[tt,1])*P2_fun(X[tt,2]),
      lgdrd[[1]][[6]](X[tt,1],X[tt,2]) - dP1_fun(X[tt,1])*P2_fun(X[tt,2]),
      lgdrd[[1]][[7]](X[tt,1],X[tt,2]) - dP2_fun(X[tt,1])*P2_fun(X[tt,2]),
      lgdrd[[1]][[8]](X[tt,1],X[tt,2]) - dP2_fun(X[tt,1])*P1_fun(X[tt,2]),
      lgdrd[[1]][[9]](X[tt,1],X[tt,2]) - dP2_fun(X[tt,1])*P0_fun(X[tt,2]),
      lgdrd[[2]][[1]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*dP0_fun(X[tt,2]),
      lgdrd[[2]][[2]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*dP1_fun(X[tt,2]),
      lgdrd[[2]][[3]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*dP1_fun(X[tt,2]),
      lgdrd[[2]][[4]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*dP0_fun(X[tt,2]),
      lgdrd[[2]][[5]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*dP2_fun(X[tt,2]),
      lgdrd[[2]][[6]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*dP2_fun(X[tt,2]),
      lgdrd[[2]][[7]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*dP2_fun(X[tt,2]),
      lgdrd[[2]][[8]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*dP1_fun(X[tt,2]),
      lgdrd[[2]][[9]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*dP0_fun(X[tt,2])
    )
  }
  
  print(apply(abs(E),2,max))
  print(apply(abs(Ed),2,max))
  print(max(abs(E)))
  range(abs(E))
  range(abs(Ed))
  dn();hist(abs(E),breaks = 100)
  dn();hist(abs(Ed),breaks = 100)
  dn();matplot(X[,1],log10(abs(E))[,7:9])
  dn();matplot(X[,1],log10(abs(Ed))[,16:18])
}

chk_nlgdr=function(){
  P0_fun <- function(x) 1/sqrt(2/(2*0+1))
  P1_fun <- function(x) x/sqrt(2/(2*1+1))
  P2_fun <- function(x) 0.5 * (3*x^2 - 1)/sqrt(2/(2*2+1))
  dP0_fun <- function(x) 0
  dP1_fun <- function(x) 1/sqrt(2/(2*1+1))
  dP2_fun <- function(x) 3*x/sqrt(2/(2*2+1))
  
  X = rbind(
    matrix(c(0,0,0,-1,-1,-1,1,1,1  ,0,-1,1,0,-1,1,0,-1,1),nrow=9,ncol=2),
    matrix(runif(n=9991*2)*2 -1,nrow=9991,ncol=2)
  )
  Tx = dim(X)[1]
  E = matrix(0,nrow=Tx,ncol=3*3)
  for(tt in 1:Tx){
    E[tt,] = c(
      nlgdr[[1]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*P0_fun(X[tt,2]),
      nlgdr[[2]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*P1_fun(X[tt,2]),
      nlgdr[[3]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*P1_fun(X[tt,2]),
      nlgdr[[4]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*P0_fun(X[tt,2]),
      nlgdr[[5]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*P2_fun(X[tt,2]),
      nlgdr[[6]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*P2_fun(X[tt,2]),
      nlgdr[[7]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*P2_fun(X[tt,2]),
      nlgdr[[8]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*P1_fun(X[tt,2]),
      nlgdr[[9]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*P0_fun(X[tt,2])
    )
  }
  Ed = matrix(0,nrow=Tx,ncol=3*3*2)
  for(tt in 1:Tx){
    Ed[tt,] = c(
      nlgdrd[[1]][[1]](X[tt,1],X[tt,2]) - dP0_fun(X[tt,1])*P0_fun(X[tt,2]),
      nlgdrd[[1]][[2]](X[tt,1],X[tt,2]) - dP0_fun(X[tt,1])*P1_fun(X[tt,2]),
      nlgdrd[[1]][[3]](X[tt,1],X[tt,2]) - dP1_fun(X[tt,1])*P1_fun(X[tt,2]),
      nlgdrd[[1]][[4]](X[tt,1],X[tt,2]) - dP1_fun(X[tt,1])*P0_fun(X[tt,2]),
      nlgdrd[[1]][[5]](X[tt,1],X[tt,2]) - dP0_fun(X[tt,1])*P2_fun(X[tt,2]),
      nlgdrd[[1]][[6]](X[tt,1],X[tt,2]) - dP1_fun(X[tt,1])*P2_fun(X[tt,2]),
      nlgdrd[[1]][[7]](X[tt,1],X[tt,2]) - dP2_fun(X[tt,1])*P2_fun(X[tt,2]),
      nlgdrd[[1]][[8]](X[tt,1],X[tt,2]) - dP2_fun(X[tt,1])*P1_fun(X[tt,2]),
      nlgdrd[[1]][[9]](X[tt,1],X[tt,2]) - dP2_fun(X[tt,1])*P0_fun(X[tt,2]),
      nlgdrd[[2]][[1]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*dP0_fun(X[tt,2]),
      nlgdrd[[2]][[2]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*dP1_fun(X[tt,2]),
      nlgdrd[[2]][[3]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*dP1_fun(X[tt,2]),
      nlgdrd[[2]][[4]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*dP0_fun(X[tt,2]),
      nlgdrd[[2]][[5]](X[tt,1],X[tt,2]) - P0_fun(X[tt,1])*dP2_fun(X[tt,2]),
      nlgdrd[[2]][[6]](X[tt,1],X[tt,2]) - P1_fun(X[tt,1])*dP2_fun(X[tt,2]),
      nlgdrd[[2]][[7]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*dP2_fun(X[tt,2]),
      nlgdrd[[2]][[8]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*dP1_fun(X[tt,2]),
      nlgdrd[[2]][[9]](X[tt,1],X[tt,2]) - P2_fun(X[tt,1])*dP0_fun(X[tt,2])
    )
  }
  
  print(apply(abs(E),2,max))
  print(apply(abs(Ed),2,max))
  print(max(abs(E)))
  range(abs(E))
  range(abs(Ed))
  dn();hist(abs(E),breaks = 100)
  dn();hist(abs(Ed),breaks = 100)
  dn();matplot(X[,1],log10(abs(E))[,7:9])
  dn();matplot(X[,1],log10(abs(Ed))[,16:18])
}


chk_nmlgdr = function(){
  
  X = rbind(
    matrix(c(0,0,0,-1,-1,-1,1,1,1  ,0,-1,1,0,-1,1,0,-1,1),nrow=9,ncol=2),
    matrix(runif(n=9991*2)*2 -1 ,nrow=9991,ncol=2)
  )
  Tx = dim(X)[1]
  Mm = 8
  n = 2
  
  phi = phiD
  phid = phiDd
  phi_name = phiD_name

  # 軌道に沿った基底関数の計算
  gphi = array(0,dim=c(Mm,Tx))
  gphid = array(0,dim=c(n,Mm,Tx))
  for(ll in 1:Mm){
    for(tt in 1:Tx){
      gphi[ll,tt] = phi[[ll]](X[tt,1],X[tt,2])
      for(kk in 1:n){
        gphid[kk,ll,tt] = phid[[kk]][[ll]](X[tt,1],X[tt,2])
      }
    }
  }
  
  gphi_ = gphi
  gphid_ = gphid
  
  phi = nmlgdr
  phid = nmlgdrd
  phi_name = nmlgdr
  
  # 軌道に沿った基底関数の計算
  gphi = array(0,dim=c(Mm,Tx))
  gphid = array(0,dim=c(n,Mm,Tx))
  for(ll in 1:Mm){
    for(tt in 1:Tx){
      gphi[ll,tt] = phi[[ll]](X[tt,1],X[tt,2])
      for(kk in 1:n){
        gphid[kk,ll,tt] = phid[[kk]][[ll]](X[tt,1],X[tt,2])
      }
    }
  }
  
  
  range(abs(gphi_ - gphi))
  range(abs(gphid_ - gphid))
  
  for(ii in 1:Mm){
    print(c(max(abs(gphi_ - gphi)[ii,]),max(abs(gphid_ - gphid)[1,ii,]),max(abs(gphid_ - gphid)[2,ii,])))
  }
  dn();hist(abs(gphi_ - gphi),breaks = 100)
  dn();hist(abs(gphid_ - gphid),breaks = 100)
  dn();matplot(X[,1],t(log10(abs(gphi_ - gphi))[5:8,]))
  
}






chk_bnst = function(){
  
  X = rbind(
    matrix(c(0,0,0,0.5,0.5,0.5,1,1,1  ,0,0.5,1,0,0.5,1,0,0.5,1),nrow=9,ncol=2),
    matrix(runif(n=9991*2),nrow=9991,ncol=2)
  )
  Tx = dim(X)[1]
  Mm = 9
  n = 2
  
  phi = phiC
  phid = phiCd
  phi_name = phiC_name
  # phi = bnst[[2]]
  # phid = bnstd[[2]]
  # phi_name = bnst_name[[2]]
  
  # 軌道に沿った基底関数の計算
  gphi = array(0,dim=c(Mm,Tx))
  gphid = array(0,dim=c(n,Mm,Tx))
  for(ll in 1:Mm){
    for(tt in 1:Tx){
      gphi[ll,tt] = phi[[ll]](X[tt,1],X[tt,2])
      for(kk in 1:n){
        gphid[kk,ll,tt] = phid[[kk]][[ll]](X[tt,1],X[tt,2])
      }
    }
  }
  
  gphi_ = gphi
  gphid_ = gphid
  
  phi = bnst[[2]]
  phid = bnstd[[2]]
  phi_name = bnst_name[[2]]
  
  # 軌道に沿った基底関数の計算
  gphi = array(0,dim=c(Mm,Tx))
  gphid = array(0,dim=c(n,Mm,Tx))
  for(ll in 1:Mm){
    for(tt in 1:Tx){
      gphi[ll,tt] = phi[[ll]](X[tt,1],X[tt,2])
      for(kk in 1:n){
        gphid[kk,ll,tt] = phid[[kk]][[ll]](X[tt,1],X[tt,2])
      }
    }
  }
  
  
  range(abs(gphi_ - gphi))
  range(abs(gphid_ - gphid))
  dn();hist(abs(gphi_ - gphi),breaks = 100)
  dn();hist(abs(gphid_ - gphid),breaks = 100)
  dn();plot(X[,1],log10(abs(gphi_ - gphi))[1,])
  dn();matplot(X[,1],t(log10(abs(gphi_ - gphi))[1:3,]))
  
  for(ii in 1:Mm){
    print(c(max(abs(gphi_ - gphi)[ii,]),max(abs(gphid_ - gphid)[1,ii,]),max(abs(gphid_ - gphid)[2,ii,])))
  }
}



