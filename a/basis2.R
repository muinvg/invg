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




phiD=list(
  function(x,y){  x},function(x,y){  y},
  function(x,y){x*x},function(x,y){x*y},function(x,y){y*y},
  function(x,y){x*x*x},function(x,y){x*x*y},function(x,y){x*y*y},function(x,y){y*y*y},
  function(x,y){x*x*x*x},function(x,y){x*x*x*y},function(x,y){x*x*y*y},function(x,y){x*y*y*y},function(x,y){y*y*y*y},
  function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y},function(x,y){x*x*x*y*y},
  function(x,y){x*x*y*y*y},function(x,y){x*y*y*y*y},function(x,y){y*y*y*y*y},
  function(x,y){x*x*x*x*x*x},function(x,y){x*x*x*x*x*y},function(x,y){x*x*x*x*y*y},
  function(x,y){x*x*x*y*y*y},function(x,y){x*x*y*y*y*y},function(x,y){x*y*y*y*y*y},function(x,y){y*y*y*y*y*y})

phiDx=list(
  function(x,y){  x*0 + 1},function(x,y){  y*0},
  function(x,y){2*x},function(x,y){y},function(x,y){0*y},
  function(x,y){3*x*x},function(x,y){2*x*y},function(x,y){y*y},function(x,y){0*y},
  function(x,y){4*x*x*x},function(x,y){3*x*x*y},function(x,y){2*x*y*y},function(x,y){y*y*y},function(x,y){0*y},
  function(x,y){5*x*x*x*x},function(x,y){4*x*x*x*y},function(x,y){3*x*x*y*y},
  function(x,y){2*x*y*y*y},function(x,y){y*y*y*y},function(x,y){0*y},
  function(x,y){6*x*x*x*x*x},function(x,y){5*x*x*x*x*y},function(x,y){4*x*x*x*y*y},
  function(x,y){3*x*x*y*y*y},function(x,y){2*x*y*y*y*y},function(x,y){y*y*y*y*y},function(x,y){0*y})

phiDy=list(
  function(x,y){  x*0},function(x,y){  y*0 + 1},
  function(x,y){0*x},function(x,y){x},function(x,y){2*y},
  function(x,y){0*x},function(x,y){x*x},function(x,y){x*y*2},function(x,y){3*y*y},
  function(x,y){0*x},function(x,y){x*x*x},function(x,y){x*x*y*2},function(x,y){x*y*y*3},function(x,y){y*y*y*4},
  function(x,y){x*0},function(x,y){x*x*x*x},function(x,y){x*x*x*y*2},
  function(x,y){x*x*y*y*3},function(x,y){x*y*y*y*4},function(x,y){y*y*y*y*5},
  function(x,y){x*0},function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y*2},
  function(x,y){x*x*x*y*y*3},function(x,y){x*x*y*y*y*4},function(x,y){x*y*y*y*y*5},function(x,y){y*y*y*y*y*6})

phiDd = list(phiDx,phiDy)

phiD_name = c("x","y","x*x","x*y","y*y",
              "x*x*x","x*x*y","x*y*y","y*y*y","x*x*x*x","x*x*x*y","x*x*y*y","x*y*y*y","y*y*y*y",
              "x*x*x*x*x","x*x*x*x*y","x*x*x*y*y","x*x*y*y*y","x*y*y*y*y","y*y*y*y*y",
              "x*x*x*x*x*x","x*x*x*x*x*y","x*x*x*x*y*y","x*x*x*y*y*y","x*x*y*y*y*y","x*y*y*y*y*y","y*y*y*y*y*y")





