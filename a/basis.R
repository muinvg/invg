phi=list(
  function(x,y){  x*0 + 1},
  function(x,y){  x},function(x,y){  y},
  function(x,y){x*x},function(x,y){x*y},function(x,y){y*y},
  function(x,y){x*x*x},function(x,y){x*x*y},function(x,y){x*y*y},function(x,y){y*y*y},
  function(x,y){x*x*x*x},function(x,y){x*x*x*y},function(x,y){x*x*y*y},function(x,y){x*y*y*y},function(x,y){y*y*y*y},
  function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y},function(x,y){x*x*x*y*y},
  function(x,y){x*x*y*y*y},function(x,y){x*y*y*y*y},function(x,y){y*y*y*y*y},
  function(x,y){x*x*x*x*x*x},function(x,y){x*x*x*x*x*y},function(x,y){x*x*x*x*y*y},
  function(x,y){x*x*x*y*y*y},function(x,y){x*x*y*y*y*y},function(x,y){x*y*y*y*y*y},function(x,y){y*y*y*y*y*y})

phix=list(
  function(x,y){  x*0 + 0},
  function(x,y){  x*0 + 1},function(x,y){  y*0},
  function(x,y){2*x},function(x,y){y},function(x,y){0*y},
  function(x,y){3*x*x},function(x,y){2*x*y},function(x,y){y*y},function(x,y){0*y},
  function(x,y){4*x*x*x},function(x,y){3*x*x*y},function(x,y){2*x*y*y},function(x,y){y*y*y},function(x,y){0*y},
  function(x,y){5*x*x*x*x},function(x,y){4*x*x*x*y},function(x,y){3*x*x*y*y},
  function(x,y){2*x*y*y*y},function(x,y){y*y*y*y},function(x,y){0*y},
  function(x,y){6*x*x*x*x*x},function(x,y){5*x*x*x*x*y},function(x,y){4*x*x*x*y*y},
  function(x,y){3*x*x*y*y*y},function(x,y){2*x*y*y*y*y},function(x,y){y*y*y*y*y},function(x,y){0*y})

phiy=list(
  function(x,y){  x*0 + 0},
  function(x,y){  x*0},function(x,y){  y*0 + 1},
  function(x,y){0*x},function(x,y){x},function(x,y){2*y},
  function(x,y){0*x},function(x,y){x*x},function(x,y){x*y*2},function(x,y){3*y*y},
  function(x,y){0*x},function(x,y){x*x*x},function(x,y){x*x*y*2},function(x,y){x*y*y*3},function(x,y){y*y*y*4},
  function(x,y){x*0},function(x,y){x*x*x*x},function(x,y){x*x*x*y*2},
  function(x,y){x*x*y*y*3},function(x,y){x*y*y*y*4},function(x,y){y*y*y*y*5},
  function(x,y){x*0},function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y*2},
  function(x,y){x*x*x*y*y*3},function(x,y){x*x*y*y*y*4},function(x,y){x*y*y*y*y*5},function(x,y){y*y*y*y*y*6})

phid = list(phix,phiy)



phi_name = c("1","x","y","x*x","x*y","y*y",
             "x*x*x","x*x*y","x*y*y","y*y*y","x*x*x*x","x*x*x*y","x*x*y*y","x*y*y*y","y*y*y*y",
             "x*x*x*x*x","x*x*x*x*y","x*x*x*y*y","x*x*y*y*y","x*y*y*y*y","y*y*y*y*y",
             "x*x*x*x*x*x","x*x*x*x*x*y","x*x*x*x*y*y","x*x*x*y*y*y","x*x*y*y*y*y","x*y*y*y*y*y","y*y*y*y*y*y")


bern1d <- function(k, n, x) { choose(n, k) * x^k * (1-x)^(n-k) }
dbern1d <- function(k, n, x) { # derivative: n * [B_{k-1,n-1}(u) - B_{k,n-1}(u)] with boundary handling
  if (n == 0) return(0*x)
    term1 <- if (k-1 >= 0) {bern1d(k-1, n-1, x)} else {0*x}
  term2 <- if (k <= n-1) {bern1d(k, n-1, x)} else {0*x}
  return( n * (term1 - term2))
}


bernstein=list(
  function(x,y){  x*0 + 1},
  function(x,y){  x},function(x,y){  y},
  function(x,y){x*x},function(x,y){x*y},function(x,y){y*y},
  function(x,y){x*x*x},function(x,y){x*x*y},function(x,y){x*y*y},function(x,y){y*y*y},
  function(x,y){x*x*x*x},function(x,y){x*x*x*y},function(x,y){x*x*y*y},function(x,y){x*y*y*y},function(x,y){y*y*y*y},
  function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y},function(x,y){x*x*x*y*y},
  function(x,y){x*x*y*y*y},function(x,y){x*y*y*y*y},function(x,y){y*y*y*y*y},
  function(x,y){x*x*x*x*x*x},function(x,y){x*x*x*x*x*y},function(x,y){x*x*x*x*y*y},
  function(x,y){x*x*x*y*y*y},function(x,y){x*x*y*y*y*y},function(x,y){x*y*y*y*y*y},function(x,y){y*y*y*y*y*y})

# phix=list(
#   function(x,y){  x*0 + 0},
#   function(x,y){  x*0 + 1},function(x,y){  y*0},
#   function(x,y){2*x},function(x,y){y},function(x,y){0*y},
#   function(x,y){3*x*x},function(x,y){2*x*y},function(x,y){y*y},function(x,y){0*y},
#   function(x,y){4*x*x*x},function(x,y){3*x*x*y},function(x,y){2*x*y*y},function(x,y){y*y*y},function(x,y){0*y},
#   function(x,y){5*x*x*x*x},function(x,y){4*x*x*x*y},function(x,y){3*x*x*y*y},
#   function(x,y){2*x*y*y*y},function(x,y){y*y*y*y},function(x,y){0*y},
#   function(x,y){6*x*x*x*x*x},function(x,y){5*x*x*x*x*y},function(x,y){4*x*x*x*y*y},
#   function(x,y){3*x*x*y*y*y},function(x,y){2*x*y*y*y*y},function(x,y){y*y*y*y*y},function(x,y){0*y})
# 
# phiy=list(
#   function(x,y){  x*0 + 0},
#   function(x,y){  x*0},function(x,y){  y*0 + 1},
#   function(x,y){0*x},function(x,y){x},function(x,y){2*y},
#   function(x,y){0*x},function(x,y){x*x},function(x,y){x*y*2},function(x,y){3*y*y},
#   function(x,y){0*x},function(x,y){x*x*x},function(x,y){x*x*y*2},function(x,y){x*y*y*3},function(x,y){y*y*y*4},
#   function(x,y){x*0},function(x,y){x*x*x*x},function(x,y){x*x*x*y*2},
#   function(x,y){x*x*y*y*3},function(x,y){x*y*y*y*4},function(x,y){y*y*y*y*5},
#   function(x,y){x*0},function(x,y){x*x*x*x*x},function(x,y){x*x*x*x*y*2},
#   function(x,y){x*x*x*y*y*3},function(x,y){x*x*y*y*y*4},function(x,y){x*y*y*y*y*5},function(x,y){y*y*y*y*y*6})
# 
# phid = list(phix,phiy)
# 
# 
# 
# phi_name = c("1","x","y","x*x","x*y","y*y",
#              "x*x*x","x*x*y","x*y*y","y*y*y","x*x*x*x","x*x*x*y","x*x*y*y","x*y*y*y","y*y*y*y",
#              "x*x*x*x*x","x*x*x*x*y","x*x*x*y*y","x*x*y*y*y","x*y*y*y*y","y*y*y*y*y",
#              "x*x*x*x*x*x","x*x*x*x*x*y","x*x*x*x*y*y","x*x*x*y*y*y","x*x*y*y*y*y","x*y*y*y*y*y","y*y*y*y*y*y")
