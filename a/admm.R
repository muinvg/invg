
admm01 = function(A,Q,r,eps,rho0=1,max_iter = 200){  
  #------------------------------
  # Problem dimensions
  #------------------------------
  n <- dim(A)[2]    # dimension of x
  m <- dim(A)[1]    # rows of A (tall matrix)
  k <- dim(Q)[1]      # rows of Q (here 1x n)
  
  r_vec <- r             # r: k-dimensional
  
  # eps <- 1.0                    # epsilon in ||Ax||_2 <= eps
  rho <- rho0                    # ADMM penalty parameter
  
  #------------------------------
  # Helper: soft-thresholding
  #------------------------------
  soft <- function(x, tau) {
    sign(x) * pmax(abs(x) - tau, 0)
  }
  
  #------------------------------
  # ADMM variables initialization
  #------------------------------
  x <- rep(0, n)     # primal variable
  z <- rep(0, n)     # for L1
  y <- rep(0, m)     # = Ax

  u <- rep(0, n)     # dual for z = x
  w <- rep(0, m)     # dual for y = Ax
  
  # # 別解法で得た解
  # x0 <- x_initial   # ここにあなたの別解法の解を入れる
  # 
  # # ADMM の整合的な初期化
  # x <- x0
  # z <- x0           # z = x を最初から満たす
  # y <- A %*% x0     # y = Ax を最初から満たす
  # 
  # u <- rep(0, n)    # 双対変数は 0 で OK
  # w <- rep(0, m)
  # 
  #------------------------------
  # Precompute KKT matrix
  #------------------------------
  # KKT system for x-update:
  #   minimize  (rho/2)||x - (z - u)||^2 + (rho/2)||A x - (y - w)||^2
  #   subject to Qx = r
  #
  # Normal equations part: (I + A^T A) x + Q^T lambda = (z - u) + A^T(y - w)
  # Constraint:            Q x                     = r
  
  I <- diag(n)
  AtA <- crossprod(A)   # A^T A  (n x n)
  M11 <- I + AtA        # (n x n)
  M12 <- t(Q)          # (n x k)
  M21 <- Q              # (k x n)
  M22 <- matrix(0, k, k)
  
  KKT <- rbind(
    cbind(M11, M12),
    cbind(M21, M22)
  )
  # KKT is (n+k) x (n+k)
  
  # For moderate n,k, we can factor once;ここでは簡単に solve を毎回使う
  # 実際には chol 分解などで高速化してもよい
  
  #------------------------------
  # ADMM iterations
  #------------------------------
  
  history <- list()
  
  for (iter in 1:max_iter) {
    
    #--------------------------
    # 1. x-update: solve KKT
    #--------------------------
    rhs1 <- (z - u) + crossprod(A, (y - w))  # n-dim
    rhs2 <- r_vec                            # k-dim
    rhs <- c(rhs1, rhs2)                     # (n+k)-dim
    
    sol <- solve(KKT, rhs)
    x_new <- sol[1:n]
    lambda_new <- sol[(n+1):(n+k)]            # not really used explicitly
    
    x <- x_new
    
    #--------------------------
    # 2. y-update: projection onto L2 ball
    #    y <- argmin_y (rho/2)||Ax - y + w||^2  s.t. ||y||_2 <= eps
    #    => projection of v = Ax + w onto {y: ||y||<=eps}
    #--------------------------
    v <- as.vector(A %*% x + w)
    v_norm <- sqrt(sum(v^2))
    
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    #--------------------------
    # 3. z-update: soft thresholding
    #    z <- argmin_z ||z||_1 + (rho/2)||x - z + u||^2
    #--------------------------
    z <- soft(x + u, 1 / rho)
    
    #--------------------------
    # 4. Dual updates
    #--------------------------
    u <- u + (x - z)
    w <- w + (A %*% x - y)
    
    #--------------------------
    # 5. Record history (optional)
    #--------------------------
    obj <- sum(abs(z))         # objective ~ ||x||_1 (z ~ x)
    feas_Ax <- sqrt(sum((A %*% x)^2))
    feas_Q <- sqrt(sum((Q %*% x - r_vec)^2))
    
    history[[iter]] <- list(
      obj = obj,
      feas_Ax = feas_Ax,
      feas_Q = feas_Q
    )
    
    if (iter %% 20 == 0) {
      cat(sprintf("iter %d: obj=%.4f, ||Ax||=%.4f, ||Qx-r||=%.4e\n",
                  iter, obj, feas_Ax, feas_Q))
    }
  }
  
  cat("ADMM finished.\n")

  return(list(x=x,h=history))
}  



admm01_fast_m_large <- function(A, Q, r, eps, rho0 = 1, max_iter = 200, reg_eps = 1e-12) {
  # A: m x n (m large ~ 40000, n small < 100)
  # Q: 1 x n (k = 1)
  n <- dim(A)[2]
  m <- dim(A)[1]
  k <- dim(Q)[1]
  if (k != 1) stop("This variant assumes k == 1 (Q has one row).")
  
  r_vec <- r
  rho <- rho0
  
  # soft-threshold
  soft <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)
  
  # ADMM variables
  x <- rep(0, n)
  z <- rep(0, n)
  y <- rep(0, m)
  u <- rep(0, n)
  w <- rep(0, m)
  
  # Precompute M11 = I + A^T A (n x n) and its Cholesky factor
  AtA <- crossprod(A)            # n x n
  M11 <- AtA
  diag(M11) <- diag(M11) + 1.0   # M11 = I + A^T A
  # ensure symmetry
  M11 <- (M11 + t(M11)) / 2
  
  cholM11 <- chol(M11)           # upper triangular R such that t(R) %*% R = M11
  
  # helper: solve M11 * z = b using cholM11 (two backsolve)
  solve_M11 <- function(b) {
    # solve R^T y = b  (transpose = TRUE)
    ytmp <- backsolve(cholM11, b, transpose = TRUE)
    # solve R z = ytmp
    z <- backsolve(cholM11, ytmp, transpose = FALSE)
    return(as.vector(z))
  }
  
  # q = Q^T (n-vector)
  q <- as.vector(t(Q))
  # precompute M11^{-1} q
  M11inv_q <- solve_M11(q)
  S <- - as.numeric(Q %*% M11inv_q)   # scalar Schur complement
  if (abs(S) < 1e-14) {
    # tiny regularization to avoid division by zero
    S <- S - sign(S) * reg_eps
    if (S == 0) S <- -reg_eps
  }
  
  history <- vector("list", max_iter)
  
  for (iter in 1:max_iter) {
    # 1) x-update via Schur complement using M11 factor
    rhs1 <- (z - u) + crossprod(A, (y - w))   # n-vector
    
    M11inv_rhs1 <- solve_M11(rhs1)
    rhs_lambda <- as.numeric(r_vec - Q %*% M11inv_rhs1)
    lambda <- rhs_lambda / S
    
    rhs_x <- rhs1 - q * lambda
    x_new <- solve_M11(rhs_x)
    x <- x_new
    
    # 2) y-update: projection onto L2 ball
    Ax = A %*% x
    v <- as.vector(Ax + w)
    v_norm <- sqrt(sum(v^2))
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    # 3) z-update: soft thresholding
    z <- soft(x + u, 1 / rho)
    
    # 4) dual updates
    u <- u + (x - z)
    w <- w + (Ax - y)
    
    # 5) history
    obj <- sum(abs(z))
    feas_Ax <- sqrt(sum((Ax)^2))
    feas_Q <- sqrt(sum((Q %*% x - r_vec)^2))
    history[[iter]] <- list(obj = obj, feas_Ax = feas_Ax, feas_Q = feas_Q)
    
    if (iter %% 20 == 0) {
      cat(sprintf("iter %d: obj=%.4f, ||Ax||=%.4f, ||Qx-r||=%.4e\n",
                  iter, obj, feas_Ax, feas_Q))
    }
  }
  
  cat("ADMM finished.\n")
  return(list(x = x, h = history))
}


admm02 = function(A,Q,r,eps,rho0=1,max_iter = 200){  
  #------------------------------
  # Problem dimensions
  #------------------------------
  n <- dim(A)[2]    # dimension of x
  m <- dim(A)[1]    # rows of A (tall matrix)
  k <- dim(Q)[1]      # rows of Q (here 1x n)
  
  r_vec <- r             # r: k-dimensional
  
  # eps <- 1.0                    # epsilon in ||Ax||_2 <= eps
  rho <- rho0                    # ADMM penalty parameter
  
  #------------------------------
  # Helper: soft-thresholding
  #------------------------------
  soft <- function(x, tau) {
    sign(x) * pmax(abs(x) - tau, 0)
  }
  
  #------------------------------
  # ADMM variables initialization
  #------------------------------
  x <- rep(0, n)     # primal variable
  z <- rep(0, n)     # for L1
  y <- rep(0, m)     # = Ax
  
  u <- rep(0, n)     # dual for z = x
  w <- rep(0, m)     # dual for y = Ax
  
  # # 別解法で得た解
  # x0 <- x_initial   # ここにあなたの別解法の解を入れる
  # 
  # # ADMM の整合的な初期化
  # x <- x0
  # z <- x0           # z = x を最初から満たす
  # y <- A %*% x0     # y = Ax を最初から満たす
  # 
  # u <- rep(0, n)    # 双対変数は 0 で OK
  # w <- rep(0, m)
  # 
  #------------------------------
  # Precompute KKT matrix
  #------------------------------
  # KKT system for x-update:
  #   minimize  (rho/2)||x - (z - u)||^2 + (rho/2)||A x - (y - w)||^2
  #   subject to Qx = r
  #
  # Normal equations part: (I + A^T A) x + Q^T lambda = (z - u) + A^T(y - w)
  # Constraint:            Q x                     = r
  
  I <- diag(n)
  AtA <- crossprod(A)   # A^T A  (n x n)
  M11 <- I + AtA        # (n x n)
  M12 <- t(Q)/rho          # (n x k)
  M21 <- Q              # (k x n)
  M22 <- matrix(0, k, k)
  
  KKT <- rbind(
    cbind(M11, M12),
    cbind(M21, M22)
  )
  # KKT is (n+k) x (n+k)
  
  # For moderate n,k, we can factor once;ここでは簡単に solve を毎回使う
  # 実際には chol 分解などで高速化してもよい
  
  #------------------------------
  # ADMM iterations
  #------------------------------
  
  history <- list()
  
  for (iter in 1:max_iter) {
    
    #--------------------------
    # 1. x-update: solve KKT
    #--------------------------
    rhs1 <- (z - u) + crossprod(A, (y - w))  # n-dim
    rhs2 <- r_vec                            # k-dim
    rhs <- c(rhs1, rhs2)                     # (n+k)-dim
    
    sol <- solve(KKT, rhs,tol=0)
    x_new <- sol[1:n]
    lambda_new <- sol[(n+1):(n+k)]            # not really used explicitly
    
    x <- x_new
    
    #--------------------------
    # 2. y-update: projection onto L2 ball
    #    y <- argmin_y (rho/2)||Ax - y + w||^2  s.t. ||y||_2 <= eps
    #    => projection of v = Ax + w onto {y: ||y||<=eps}
    #--------------------------
    Ax = A %*% x
    v <- as.vector(Ax + w)
    v_norm <- sqrt(sum(v^2))
    
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    #--------------------------
    # 3. z-update: soft thresholding
    #    z <- argmin_z ||z||_1 + (rho/2)||x - z + u||^2
    #--------------------------
    z <- soft(x + u, 1 / rho)
    
    #--------------------------
    # 4. Dual updates
    #--------------------------
    u <- u + (x - z)
    w <- w + (Ax - y)
    
    #--------------------------
    # 5. Record history (optional)
    #--------------------------
    obj <- sum(abs(z))         # objective ~ ||x||_1 (z ~ x)
    feas_Ax <- sqrt(sum((Ax)^2))
    feas_Q <- sqrt(sum((Q %*% x - r_vec)^2))
    
    history[[iter]] <- list(
      obj = obj,
      feas_Ax = feas_Ax,
      feas_Q = feas_Q
    )
    
    if (iter %% 20 == 0) {
      cat(sprintf("iter %d: obj=%.4f, ||Ax||=%.4f, ||Qx-r||=%.4e\n",
                  iter, obj, feas_Ax, feas_Q))
    }
  }
  
  cat("ADMM finished.\n")
  
  return(list(x=x,h=history))
}  



admm_with_psd <- function(A, Q, r, eps, B_list,
                          rho = 1.0, rho_s = 1.0,
                          max_iter = 500, eps_abs = 1e-4, eps_rel = 1e-3,
                          verbose = TRUE) {
  # A: m x n
  # Q: k x n
  # r: k-vector
  # eps: scalar for ||Ax||_2 <= eps
  # B_list: list of length K, each is 3 x n matrix with order [s11, s12, s22]
  # rho: ADMM penalty for x-z and Ax-y parts
  # rho_s: ADMM penalty for PSD constraints (can be = rho)
  # Returns list(x, history)
  
  n <- ncol(A)
  m <- nrow(A)
  kq <- nrow(Q)
  K <- length(B_list)
  r_vec <- as.vector(r)
  
  # soft-threshold
  soft <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)
  
  # initialize primal/dual variables
  x <- rep(0, n)
  z <- rep(0, n)
  y <- rep(0, m)        # for Ax
  u <- rep(0, n)        # dual for z = x
  w <- rep(0, m)        # dual for y = Ax
  
  # PSD auxiliary variables and duals
  s_list <- vector("list", K)
  y_s_list <- vector("list", K)
  for (kk in 1:K) {
    s_list[[kk]] <- rep(0, 3)
    y_s_list[[kk]] <- rep(0, 3)
  }
  
  # Precompute some matrices
  AtA <- crossprod(A)   # n x n
  M_base <- diag(n) + AtA   # I + A^T A
  
  
  
  history <- list()
  
  for (iter in 1:max_iter) {
    # -----------------------
    # x-update: assemble LHS and RHS for KKT
    # We use scaled form: multiply (I + A^T A) by rho to avoid 1/rho fractions
    # LHS_x = rho * (I + A^T A) + sum_k rho_s * Bk^T Bk
    # RHS_x = rho * ((z - u) + A^T (y - w)) + sum_k rho_s * Bk^T (s_k - y_k)
    # KKT: [LHS_x  Q^T; Q  0] [x; lambda] = [RHS_x; r]
    # -----------------------
    PSD_term <- matrix(0, n, n)
    rhs_extra <- rep(0, n)
    Bkx = as.list(1:K)
    for (kk in 1:K) {
      Bk <- B_list[[kk]]            # 3 x n
      PSD_term <- PSD_term + rho_s * crossprod(Bk)   # n x n
      rhs_extra <- rhs_extra + rho_s * as.vector(crossprod(Bk, (s_list[[kk]] - y_s_list[[kk]])))
      Bkx[[kk]] = Bk %*% x
    }
    
    LHS_x <- rho * M_base + PSD_term
    RHS_x <- rho * ((z - u) + crossprod(A, (y - w))) + rhs_extra
    
    # Assemble KKT
    M11 <- LHS_x
    M12 <- t(Q)   # scaled form uses Q^T (no 1/rho)
    M21 <- Q
    M22 <- matrix(0, nrow = kq, ncol = kq)
    
    KKT <- rbind(cbind(M11, M12), cbind(M21, M22))
    rhs_kkt <- c(RHS_x, r_vec)
    
    sol <- tryCatch({
      solve(KKT, rhs_kkt,tol=0)
    }, error = function(e) {
      stop("KKT solve failed: ", e$message)
    })
    x <- as.vector(sol[1:n])
    # lambda <- sol[(n+1):(n+kq)]  # not used further here
    
    # -----------------------
    # y-update: projection of v = A x + w onto L2 ball ||y|| <= eps
    # -----------------------
    Ax = A %*% x
    v <- as.vector(Ax + w)
    v_norm <- sqrt(sum(v^2))
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    # -----------------------
    # z-update: soft thresholding
    # -----------------------
    z_old <- z
    z <- soft(x + u, 1 / rho)
    
    # -----------------------
    # s_k update: PSD projection for each k
    # order of s vector is [s11, s12, s22]
    # reconstruct 2x2 symmetric matrix M = [s11 s12; s12 s22]
    # project M onto PSD cone via eigen decomposition
    # -----------------------
    s_old_list <- s_list
    for (kk in 1:K) {
      # Bk <- B_list[[kk]] #Bkx[[kk]]
      # b_tilde <- as.vector(Bk %*% x + y_s_list[[kk]])   # length 3
      b_tilde <- as.vector(Bkx[[kk]] + y_s_list[[kk]])   # length 3
      Mmat <- matrix(c(b_tilde[1], b_tilde[2], b_tilde[2], b_tilde[3]), nrow = 2, byrow = TRUE)
      ev <- eigen(Mmat, symmetric = TRUE)
      Dpos <- pmax(ev$values, 0)
      Mpos <- ev$vectors %*% diag(Dpos) %*% t(ev$vectors)
      s_list[[kk]] <- c(Mpos[1,1], Mpos[1,2], Mpos[2,2])  # [s11, s12, s22]
    }
    
    # -----------------------
    # dual updates for PSD constraints
    # -----------------------
    for (kk in 1:K) {
      # y_s_list[[kk]] <- y_s_list[[kk]] + as.vector(B_list[[kk]] %*% x - s_list[[kk]])
      y_s_list[[kk]] <- y_s_list[[kk]] + as.vector(Bkx[[kk]] - s_list[[kk]])
    }
    
    # -----------------------
    # dual updates for z and w
    # -----------------------
    u <- u + (x - z)
    w <- w + (Ax - y)
    
    # -----------------------
    # compute residuals and history
    # primal residuals: r1 = x - z, r2 = A x - y, r3_k = Bk x - s_k
    # dual residuals: s_z = rho * (z - z_old), s_s = rho_s * (s - s_old)
    # -----------------------
    r_prim_xz <- sqrt(sum((x - z)^2))
    r_prim_Ay <- sqrt(sum((Ax - y)^2))
    r_prim_Bs <- 0
    for (kk in 1:K) {
      # r_prim_Bs <- max(r_prim_Bs, sqrt(sum((B_list[[kk]] %*% x - s_list[[kk]])^2)))
      r_prim_Bs <- max(r_prim_Bs, sqrt(sum((Bkx[[kk]] - s_list[[kk]])^2)))
    }
    s_dual_z <- rho * sqrt(sum((z - z_old)^2))
    s_dual_s <- 0
    for (kk in 1:K) {
      s_dual_s <- max(s_dual_s, rho_s * sqrt(sum((s_list[[kk]] - s_old_list[[kk]])^2)))
    }
    
    obj_val <- sum(abs(z))   # approximate objective
    
    history[[iter]] <- list(obj = obj_val,
                            r_xz = r_prim_xz,
                            r_Ay = r_prim_Ay,
                            r_Bs = r_prim_Bs,
                            s_z = s_dual_z,
                            s_s = s_dual_s)
    
    # stopping criteria (use max of primal residuals and max of dual residuals)
    r_prim_max <- max(r_prim_xz, r_prim_Ay, r_prim_Bs)
    s_dual_max <- max(s_dual_z, s_dual_s)
    
    eps_primal <- sqrt(n + m + 3*K) * eps_abs + eps_rel * max(sqrt(sum(x^2)), sqrt(sum(z^2)), sqrt(sum(y^2)))
    eps_dual <- sqrt(n + 3*K) * eps_abs + eps_rel * sqrt(sum(u^2) + sum(unlist(y_s_list)^2))
    
    if (verbose && (iter %% 50 == 0 || iter == 1)) {
      cat(sprintf("iter %d: obj=%.6f, r_prim=%.4e, s_dual=%.4e\n", iter, obj_val, r_prim_max, s_dual_max))
    }
    
    if (r_prim_max <= eps_primal && s_dual_max <= eps_dual) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      break
    }
  } # end for iter
  
  return(list(x = x, history = history, iter = iter))
}

admm_with_psd_mode <- function(A, Q, r, eps, B_list,
                          rho = 1.0, rho_s = 1.0,
                          max_iter = 500, eps_abs = 1e-4, eps_rel = 1e-3,
                          verbose = TRUE,
                          scale_mode = c("transform", "normalize", "none"),
                          col_eps = 1e-12) {
  # A: m x n
  # Q: k x n
  # r: k-vector
  # eps: scalar for ||Ax||_2 <= eps
  # B_list: list of length K, each is 3 x n matrix with order [s11, s12, s22]
  # rho: ADMM penalty for x-z and Ax-y parts
  # rho_s: ADMM penalty for PSD constraints (can be = rho)
  # scale_mode: "none", "normalize", or "transform" (variable transform with weighted L1)
  # col_eps: small floor to avoid division by zero
  #
  # Returns list(x = solution in original scale, history, iter, D (diagonal scaling vector))
  
  scale_mode <- match.arg(scale_mode)
  
  n <- ncol(A)
  m <- nrow(A)
  kq <- nrow(Q)
  K <- length(B_list)
  r_vec <- as.vector(r)
  
  # compute column norms (squared sums provided by user earlier)
  col_norm2 <- sqrt(pmax(colSums(A^2), col_eps))
  # diagonal scaling vector d_vec such that we will multiply columns by d_vec (i.e., D = diag(d_vec))
  # For transform mode we want D = diag(1/col_norm) to make columns ~ unit norm
  d_vec <- 1 / col_norm2
  
  # Build scaled matrices depending on mode
  if (scale_mode == "none") {
    A_s <- A
    Q_s <- Q
    B_list_s <- B_list
    D_vec <- rep(1, n)   # identity scaling
    transform_mode <- FALSE
  } else if (scale_mode == "normalize") {
    # scale columns of A (and B) but keep variables x in same meaning (no weighted L1)
    A_s <- sweep(A, 2, col_norm2, FUN = "/")   # A %*% diag(1/col_norm)
    Q_s <- sweep(Q, 2, col_norm2, FUN = "/")
    B_list_s <- lapply(B_list, function(B) B %*% diag(1/col_norm2))
    D_vec <- rep(1, n)   # no variable transform; final x is directly solution
    transform_mode <- FALSE
  } else { # "transform"
    # variable transform: x = D x', with D = diag(1/col_norm2)
    # Solve for x' with weighted L1 prox; scale A, Q, B accordingly: A_s = A %*% D
    A_s <- A %*% diag(1/col_norm2)
    Q_s <- Q %*% diag(1/col_norm2)
    B_list_s <- lapply(B_list, function(B) B %*% diag(1/col_norm2))
    D_vec <- 1/col_norm2   # so x = D_vec * x_prime (elementwise)
    transform_mode <- TRUE
  }
  
  # soft-threshold (scalar threshold or vector threshold)
  soft <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)
  soft_weighted <- function(x, thresh_vec) {
    sign(x) * pmax(abs(x) - thresh_vec, 0)
  }
  
  # initialize primal/dual variables (work in x' space if transform_mode)
  x <- rep(0, n)        # this is x' when transform_mode TRUE, else actual x
  z <- rep(0, n)
  y <- rep(0, m)        # for A_s x
  u <- rep(0, n)        # dual for z = x
  w <- rep(0, m)        # dual for y = A_s x
  
  # PSD auxiliary variables and duals (operate with B_list_s)
  s_list <- vector("list", K)
  y_s_list <- vector("list", K)
  for (kk in 1:K) {
    s_list[[kk]] <- rep(0, 3)
    y_s_list[[kk]] <- rep(0, 3)
  }
  
  # Precompute some matrices for scaled system
  AtA <- crossprod(A_s)   # n x n
  M_base <- diag(n) + AtA   # I + A_s^T A_s
  
  history <- list()
  
  for (iter in 1:max_iter) {
    # -----------------------
    # x-update: assemble LHS and RHS for KKT
    # LHS_x = rho * (I + A_s^T A_s) + sum_k rho_s * Bk_s^T Bk_s
    # RHS_x = rho * ((z - u) + A_s^T (y - w)) + sum_k rho_s * Bk_s^T (s_k - y_k)
    # KKT: [LHS_x  Q_s^T; Q_s  0] [x; lambda] = [RHS_x; r]
    # -----------------------
    PSD_term <- matrix(0, n, n)
    rhs_extra <- rep(0, n)
    Bkx <- vector("list", K)
    for (kk in 1:K) {
      Bk <- B_list_s[[kk]]            # 3 x n
      PSD_term <- PSD_term + rho_s * crossprod(Bk)   # n x n
      rhs_extra <- rhs_extra + rho_s * as.vector(crossprod(Bk, (s_list[[kk]] - y_s_list[[kk]])))
      Bkx[[kk]] <- as.vector(Bk %*% x)
    }
    
    LHS_x <- rho * M_base + PSD_term
    RHS_x <- rho * ((z - u) + crossprod(A_s, (y - w))) + rhs_extra
    
    # Assemble KKT
    M11 <- LHS_x
    M12 <- t(Q_s)   # Q_s is scaled Q
    M21 <- Q_s
    M22 <- matrix(0, nrow = kq, ncol = kq)
    
    KKT <- rbind(cbind(M11, M12), cbind(M21, M22))
    rhs_kkt <- c(RHS_x, r_vec)
    
    sol <- tryCatch({
      solve(KKT, rhs_kkt)
    }, error = function(e) {
      stop("KKT solve failed: ", e$message)
    })
    x <- as.vector(sol[1:n])
    # lambda <- sol[(n+1):(n+kq)]  # not used further here
    
    # -----------------------
    # y-update: projection of v = A_s x + w onto L2 ball ||y|| <= eps
    # -----------------------
    Ax <- A_s %*% x
    v <- as.vector(Ax + w)
    v_norm <- sqrt(sum(v^2))
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    # -----------------------
    # z-update: soft thresholding
    # If transform_mode, objective is ||D x'||_1 = sum_j D_vec[j] * |x'_j|
    # so threshold_j = (1/rho) * D_vec[j]
    # If not transform_mode and scale_mode == "normalize" or "none", use standard soft
    # -----------------------
    z_old <- z
    if (transform_mode) {
      thresh_vec <- (1 / rho) * D_vec
      z <- soft_weighted(x + u, thresh_vec)
    } else {
      z <- soft(x + u, 1 / rho)
    }
    
    # -----------------------
    # s_k update: PSD projection for each k
    # order of s vector is [s11, s12, s22]
    # reconstruct 2x2 symmetric matrix M = [s11 s12; s12 s22]
    # project M onto PSD cone via eigen decomposition
    # -----------------------
    s_old_list <- s_list
    for (kk in 1:K) {
      b_tilde <- as.vector(Bkx[[kk]] + y_s_list[[kk]])   # length 3
      Mmat <- matrix(c(b_tilde[1], b_tilde[2], b_tilde[2], b_tilde[3]), nrow = 2, byrow = TRUE)
      ev <- eigen(Mmat, symmetric = TRUE)
      Dpos <- pmax(ev$values, 0)
      Mpos <- ev$vectors %*% diag(Dpos) %*% t(ev$vectors)
      s_list[[kk]] <- c(Mpos[1,1], Mpos[1,2], Mpos[2,2])  # [s11, s12, s22]
    }
    
    # -----------------------
    # dual updates for PSD constraints
    # -----------------------
    for (kk in 1:K) {
      y_s_list[[kk]] <- y_s_list[[kk]] + as.vector(Bkx[[kk]] - s_list[[kk]])
    }
    
    # -----------------------
    # dual updates for z and w
    # -----------------------
    u <- u + (x - z)
    w <- w + (Ax - y)
    
    # -----------------------
    # compute residuals and history
    # -----------------------
    r_prim_xz <- sqrt(sum((x - z)^2))
    r_prim_Ay <- sqrt(sum((Ax - y)^2))
    r_prim_Bs <- 0
    for (kk in 1:K) {
      r_prim_Bs <- max(r_prim_Bs, sqrt(sum((Bkx[[kk]] - s_list[[kk]])^2)))
    }
    s_dual_z <- rho * sqrt(sum((z - z_old)^2))
    s_dual_s <- 0
    for (kk in 1:K) {
      s_dual_s <- max(s_dual_s, rho_s * sqrt(sum((s_list[[kk]] - s_old_list[[kk]])^2)))
    }
    
    # objective: if transform_mode objective is sum_j D_vec[j] * |x_j|
    if (transform_mode) {
      obj_val <- sum(D_vec * abs(z))
    } else {
      obj_val <- sum(abs(z))
    }
    
    history[[iter]] <- list(obj = obj_val,
                            r_xz = r_prim_xz,
                            r_Ay = r_prim_Ay,
                            r_Bs = r_prim_Bs,
                            s_z = s_dual_z,
                            s_s = s_dual_s)
    
    # stopping criteria
    r_prim_max <- max(r_prim_xz, r_prim_Ay, r_prim_Bs)
    s_dual_max <- max(s_dual_z, s_dual_s)
    
    eps_primal <- sqrt(n + m + 3*K) * eps_abs + eps_rel * max(sqrt(sum(x^2)), sqrt(sum(z^2)), sqrt(sum(y^2)))
    eps_dual <- sqrt(n + 3*K) * eps_abs + eps_rel * sqrt(sum(u^2) + sum(unlist(y_s_list)^2))
    
    if (verbose && (iter %% 50 == 0 || iter == 1)) {
      cat(sprintf("iter %d: obj=%.6f, r_prim=%.4e, s_dual=%.4e\n", iter, obj_val, r_prim_max, s_dual_max))
    }
    
    if (r_prim_max <= eps_primal && s_dual_max <= eps_dual) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      break
    }
  } # end for iter
  
  # map back to original x if transform_mode
  if (transform_mode) {
    x_original <- D_vec * x   # elementwise scaling
  } else if (scale_mode == "normalize") {
    # if normalized columns but no variable transform, we must undo column scaling:
    # we solved for x (same meaning) but A was scaled; to get original-sense x we need to divide by column scaling
    # In normalize mode we scaled A <- A / col_norm, so original x = x (no change) if we treated objective as ||x||_1.
    # However if user expects original-scale coefficients relative to original A, multiply:
    x_original <- x / col_norm2   # because A_s = A / col_norm => A x_original = A_s (col_norm * x_original)
    # Note: this branch depends on interpretation; keep for completeness.
  } else {
    x_original <- x
  }
  
  return(list(x = as.vector(x_original), history = history, iter = iter, D_vec = D_vec, scale_mode = scale_mode))
}


admm_with_psd_mode_02 <- function(A, Q, r, eps, B_list,
                               rho = 1.0, rho_s = 1.0,
                               max_iter = 1000, eps_abs = 1e-4, eps_rel = 1e-3,
                               verbose = TRUE,
                               scale_mode = c("transform", "normalize", "none"),
                               col_eps = 1e-12,
                               # adaptive options (追加引数)
                               adapt_method = c("none", "residual", "adaptive"),
                               update_freq = 10,
                               tau = 2.0,
                               mu = 10.0,
                               rho_bounds = c(1e-6, 1e6),
                               rho_s_bounds = c(1e-6, 1e6),
                               eps_update = 1e-12) {
  # A: m x n
  # Q: k x n
  # r: k-vector
  # eps: scalar for ||Ax||_2 <= eps
  # B_list: list of length K, each is 3 x n matrix with order [s11, s12, s22]
  # rho: ADMM penalty for x-z and Ax-y parts
  # rho_s: ADMM penalty for PSD constraints (can be = rho)
  # scale_mode: "none", "normalize", or "transform" (variable transform with weighted L1)
  # col_eps: small floor to avoid division by zero
  # adapt_method: "none" (default), "residual" (residual balancing), "adaptive" (smooth ratio)
  # update_freq: how often (in iterations) to consider updating rho/rho_s
  # tau, mu: parameters for residual balancing
  # rho_bounds, rho_s_bounds: lower/upper bounds for rho and rho_s
  # eps_update: small floor to avoid division by zero in updates
  #
  # Returns list(x = solution in original scale, history, iter, D_vec, final_rho, final_rho_s, scale_mode)
  
  scale_mode <- match.arg(scale_mode)
  adapt_method <- match.arg(adapt_method)
  
  ## basic validation for rho/rho_s
  if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho)) stop("rho must be a finite numeric scalar")
  if (!is.numeric(rho_s) || length(rho_s) != 1 || !is.finite(rho_s)) stop("rho_s must be a finite numeric scalar")
  
  n <- ncol(A)
  m <- nrow(A)
  kq <- nrow(Q)
  K <- length(B_list)
  r_vec <- as.vector(r)
  
  # compute column norms (squared sums provided by user earlier)
  col_norm2 <- sqrt(pmax(colSums(A^2), col_eps))
  # diagonal scaling vector d_vec such that we will multiply columns by d_vec (i.e., D = diag(d_vec))
  # For transform mode we want D = diag(1/col_norm) to make columns ~ unit norm
  d_vec <- 1 / col_norm2
  
  # Build scaled matrices depending on mode
  if (scale_mode == "none") {
    A_s <- A
    Q_s <- Q
    B_list_s <- B_list
    D_vec <- rep(1, n)   # identity scaling
    transform_mode <- FALSE
  } else if (scale_mode == "normalize") {
    # scale columns of A (and B) but keep variables x in same meaning (no weighted L1)
    A_s <- sweep(A, 2, col_norm2, FUN = "/")   # A %*% diag(1/col_norm)
    Q_s <- sweep(Q, 2, col_norm2, FUN = "/")
    B_list_s <- lapply(B_list, function(B) as.matrix(B) %*% diag(1/col_norm2))
    D_vec <- rep(1, n)   # no variable transform; final x is directly solution
    transform_mode <- FALSE
  } else { # "transform"
    # variable transform: x = D x', with D = diag(1/col_norm2)
    # Solve for x' with weighted L1 prox; scale A, Q, B accordingly: A_s = A %*% D
    A_s <- A %*% diag(1/col_norm2)
    Q_s <- Q %*% diag(1/col_norm2)
    B_list_s <- lapply(B_list, function(B) as.matrix(B) %*% diag(1/col_norm2))
    D_vec <- 1/col_norm2   # so x = D_vec * x_prime (elementwise)
    transform_mode <- TRUE
  }
  
  # soft-threshold (scalar threshold or vector threshold)
  soft <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)
  soft_weighted <- function(x, thresh_vec) {
    sign(x) * pmax(abs(x) - thresh_vec, 0)
  }
  
  # initialize primal/dual variables (work in x' space if transform_mode)
  x <- rep(0, n)        # this is x' when transform_mode TRUE, else actual x
  z <- rep(0, n)
  y <- rep(0, m)        # for A_s x
  u <- rep(0, n)        # dual for z = x (scaled dual)
  w <- rep(0, m)        # dual for y = A_s x (scaled dual)
  
  # PSD auxiliary variables and duals (operate with B_list_s)
  s_list <- vector("list", K)
  y_s_list <- vector("list", K)
  for (kk in 1:K) {
    # ensure numeric vectors of correct length (assume each Bk has 3 rows as documented)
    Bk <- as.matrix(B_list_s[[kk]])
    p <- nrow(Bk)
    s_list[[kk]] <- rep(0, p)
    y_s_list[[kk]] <- rep(0, p)
  }
  
  # Precompute some matrices for scaled system
  AtA <- crossprod(A_s)   # n x n
  M_base <- diag(n) + AtA   # I + A_s^T A_s
  
  # helper: scale scaled-dual variables when rho or rho_s changes
  scale_duals_for_rho_change <- function(u, w, y_s_list, rho_old, rho_new) {
    if (rho_old == rho_new) return(list(u = u, w = w, y_s_list = y_s_list))
    factor <- rho_old / rho_new
    u <- u * factor
    w <- w * factor
    y_s_list <- lapply(y_s_list, function(ys) ys * factor)
    list(u = u, w = w, y_s_list = y_s_list)
  }
  scale_duals_for_rho_s_change <- function(y_s_list, rho_s_old, rho_s_new) {
    if (rho_s_old == rho_s_new) return(y_s_list)
    factor <- rho_s_old / rho_s_new
    lapply(y_s_list, function(ys) ys * factor)
  }
  
  history <- list()
  
  for (iter in 1:max_iter) {
    # -----------------------
    # x-update: assemble LHS and RHS for KKT
    # LHS_x = rho * (I + A_s^T A_s) + sum_k rho_s * Bk_s^T Bk_s
    # RHS_x = rho * ((z - u) + A_s^T (y - w)) + sum_k rho_s * Bk_s^T (s_k - y_k)
    # KKT: [LHS_x  Q_s^T; Q_s  0] [x; lambda] = [RHS_x; r]
    # -----------------------
    PSD_term <- matrix(0, n, n)
    rhs_extra <- rep(0, n)
    Bkx <- vector("list", K)
    for (kk in 1:K) {
      Bk <- as.matrix(B_list_s[[kk]])            # p x n
      PSD_term <- PSD_term + rho_s * crossprod(Bk)   # n x n
      # ensure vector subtraction yields numeric vector
      vec_sy <- as.numeric(s_list[[kk]] - y_s_list[[kk]])
      rhs_extra <- rhs_extra + rho_s * as.numeric(crossprod(Bk, vec_sy))
      Bkx[[kk]] <- as.numeric(Bk %*% x)
    }
    
    LHS_x <- rho * M_base + PSD_term
    RHS_x <- rho * ((z - u) + crossprod(A_s, (y - w))) + rhs_extra
    
    # Assemble KKT
    M11 <- LHS_x
    M12 <- t(Q_s)   # Q_s is scaled Q
    M21 <- Q_s
    M22 <- matrix(0, nrow = kq, ncol = kq)
    
    KKT <- rbind(cbind(M11, M12), cbind(M21, M22))
    rhs_kkt <- c(RHS_x, r_vec)
    
    sol <- tryCatch({
      solve(KKT, rhs_kkt)
    }, error = function(e) {
      stop("KKT solve failed: ", e$message)
    })
    x <- as.vector(sol[1:n])
    # lambda <- sol[(n+1):(n+kq)]  # not used further here
    
    # -----------------------
    # y-update: projection of v = A_s x + w onto L2 ball ||y|| <= eps
    # -----------------------
    Ax <- A_s %*% x
    v <- as.vector(Ax + w)
    v_norm <- sqrt(sum(v^2))
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    # -----------------------
    # z-update: soft thresholding
    # If transform_mode, objective is ||D x'||_1 = sum_j D_vec[j] * |x'_j|
    # so threshold_j = (1/rho) * D_vec[j]
    # If not transform_mode and scale_mode == "normalize" or "none", use standard soft
    # -----------------------
    z_old <- z
    if (transform_mode) {
      thresh_vec <- (1 / rho) * D_vec
      z <- soft_weighted(x + u, thresh_vec)
    } else {
      z <- soft(x + u, 1 / rho)
    }
    
    # -----------------------
    # s_k update: PSD projection for each k
    # order of s vector is [s11, s12, s22]
    # reconstruct 2x2 symmetric matrix M = [s11 s12; s12 s22]
    # project M onto PSD cone via eigen decomposition
    # -----------------------
    s_old_list <- s_list
    for (kk in 1:K) {
      b_tilde <- as.vector(Bkx[[kk]] + y_s_list[[kk]])   # length p (expected 3)
      if (length(b_tilde) == 3) {
        Mmat <- matrix(c(b_tilde[1], b_tilde[2], b_tilde[2], b_tilde[3]), nrow = 2, byrow = TRUE)
        ev <- eigen(Mmat, symmetric = TRUE)
        Dpos <- pmax(ev$values, 0)
        Mpos <- ev$vectors %*% diag(Dpos) %*% t(ev$vectors)
        s_list[[kk]] <- c(Mpos[1,1], Mpos[1,2], Mpos[2,2])  # [s11, s12, s22]
      } else {
        # fallback: if Bk produces different packing, keep b_tilde (user must ensure packing)
        s_list[[kk]] <- b_tilde
      }
    }
    
    # -----------------------
    # dual updates for PSD constraints
    # -----------------------
    for (kk in 1:K) {
      y_s_list[[kk]] <- y_s_list[[kk]] + as.vector(Bkx[[kk]] - s_list[[kk]])
    }
    
    # -----------------------
    # dual updates for z and w
    # -----------------------
    u <- u + (x - z)
    w <- w + (Ax - y)
    
    # -----------------------
    # compute residuals and history
    # -----------------------
    r_prim_xz <- sqrt(sum((x - z)^2))
    r_prim_Ay <- sqrt(sum((Ax - y)^2))
    r_prim_Bs <- 0
    for (kk in 1:K) {
      r_prim_Bs <- max(r_prim_Bs, sqrt(sum((Bkx[[kk]] - s_list[[kk]])^2)))
    }
    s_dual_z <- rho * sqrt(sum((z - z_old)^2))
    s_dual_s <- 0
    for (kk in 1:K) {
      s_dual_s <- max(s_dual_s, rho_s * sqrt(sum((s_list[[kk]] - s_old_list[[kk]])^2)))
    }
    
    # objective: if transform_mode objective is sum_j D_vec[j] * |x_j|
    if (transform_mode) {
      obj_val <- sum(D_vec * abs(z))
    } else {
      obj_val <- sum(abs(z))
    }
    
    history[[iter]] <- list(obj = obj_val,
                            r_xz = r_prim_xz,
                            r_Ay = r_prim_Ay,
                            r_Bs = r_prim_Bs,
                            s_z = s_dual_z,
                            s_s = s_dual_s,
                            rho = rho,
                            rho_s = rho_s)
    
    # -----------------------
    # adaptive updates for rho and rho_s (every update_freq iterations)
    # -----------------------
    if (adapt_method != "none" && (iter %% update_freq == 0)) {
      # aggregate primal and dual residuals for rho
      r_primal_rho <- max(r_prim_xz, r_prim_Ay)
      s_dual_rho <- s_dual_z
      
      rho_old <- rho
      rho_s_old <- rho_s
      
      if (adapt_method == "residual") {
        # classic residual balancing (multiplicative)
        if (r_primal_rho > mu * s_dual_rho + eps_update) {
          rho <- min(rho * tau, rho_bounds[2])
        } else if (s_dual_rho > mu * r_primal_rho + eps_update) {
          rho <- max(rho / tau, rho_bounds[1])
        }
        # PSD-specific: use r_prim_Bs and s_dual_s
        if (r_prim_Bs > mu * s_dual_s + eps_update) {
          rho_s <- min(rho_s * tau, rho_s_bounds[2])
        } else if (s_dual_s > mu * r_prim_Bs + eps_update) {
          rho_s <- max(rho_s / tau, rho_s_bounds[1])
        }
      } else if (adapt_method == "adaptive") {
        # smoother scaling: scale by sqrt(r_primal / s_dual)
        factor_rho <- sqrt((r_primal_rho + eps_update) / (s_dual_rho + eps_update))
        factor_rho <- min(max(factor_rho, 1 / tau), tau)
        rho <- min(max(rho * factor_rho, rho_bounds[1]), rho_bounds[2])
        
        factor_rho_s <- sqrt((r_prim_Bs + eps_update) / (s_dual_s + eps_update))
        factor_rho_s <- min(max(factor_rho_s, 1 / tau), tau)
        rho_s <- min(max(rho_s * factor_rho_s, rho_s_bounds[1]), rho_s_bounds[2])
      }
      
      # if changed, scale scaled-dual variables to keep consistency:
      if (abs(rho - rho_old) > 0) {
        scaled <- scale_duals_for_rho_change(u, w, y_s_list, rho_old, rho)
        u <- scaled$u; w <- scaled$w; y_s_list <- scaled$y_s_list
      }
      if (abs(rho_s - rho_s_old) > 0) {
        y_s_list <- scale_duals_for_rho_s_change(y_s_list, rho_s_old, rho_s)
      }
    }
    
    # stopping criteria
    r_prim_max <- max(r_prim_xz, r_prim_Ay, r_prim_Bs)
    s_dual_max <- max(s_dual_z, s_dual_s)
    
    eps_primal <- sqrt(n + m + 3*K) * eps_abs + eps_rel * max(sqrt(sum(x^2)), sqrt(sum(z^2)), sqrt(sum(y^2)))
    eps_dual <- sqrt(n + 3*K) * eps_abs + eps_rel * sqrt(sum(u^2) + sum(unlist(y_s_list)^2))
    
    if (verbose && (iter %% 50 == 0 || iter == 1)) {
      cat(sprintf("iter %d: obj=%.6f, r_prim=%.4e, s_dual=%.4e, rho=%.3e, rho_s=%.3e\n",
                  iter, obj_val, r_prim_max, s_dual_max, rho, rho_s))
    }
    
    if (r_prim_max <= eps_primal && s_dual_max <= eps_dual) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      break
    }
  } # end for iter
  
  # map back to original x if transform_mode
  if (transform_mode) {
    x_original <- D_vec * x   # elementwise scaling
  } else if (scale_mode == "normalize") {
    x_original <- x / col_norm2
  } else {
    x_original <- x
  }
  
  return(list(x = as.vector(x_original),
              history = history,
              iter = iter,
              D_vec = D_vec,
              scale_mode = scale_mode,
              final_rho = rho,
              final_rho_s = rho_s))
}


admm_with_psd_mode_l2 <- function(A, Q, r, eps, B_list,
                                  rho = 1.0, rho_s = 1.0,
                                  max_iter = 500, eps_abs = 1e-4, eps_rel = 1e-3,
                                  verbose = TRUE,
                                  scale_mode = c("transform", "normalize", "none"),
                                  col_eps = 1e-12) {
  # A: m x n
  # Q: k x n
  # r: k-vector
  # eps: scalar for ||Ax||_2 <= eps
  # B_list: list of length K, each is 3 x n matrix with order [s11, s12, s22]
  # rho: ADMM penalty for x-z and Ax-y parts
  # rho_s: ADMM penalty for PSD constraints (can be = rho)
  # scale_mode: "none", "normalize", or "transform" (variable transform)
  # col_eps: small floor to avoid division by zero
  #
  # Objective assumed: (1/2) * ||x||_2^2
  # If transform_mode: objective is (1/2) * ||D x'||_2^2 = (1/2) sum_j (D_vec[j]^2 * x'_j^2)
  #
  # Returns list(x = solution in original scale, history, iter, D_vec)
  
  scale_mode <- match.arg(scale_mode)
  
  n <- ncol(A)
  m <- nrow(A)
  kq <- nrow(Q)
  K <- length(B_list)
  r_vec <- as.vector(r)
  
  # compute column norms
  col_norm2 <- sqrt(pmax(colSums(A^2), col_eps))
  d_vec <- 1 / col_norm2
  
  # Build scaled matrices depending on mode
  if (scale_mode == "none") {
    A_s <- A
    Q_s <- Q
    B_list_s <- B_list
    D_vec <- rep(1, n)
    transform_mode <- FALSE
  } else if (scale_mode == "normalize") {
    A_s <- sweep(A, 2, col_norm2, FUN = "/")
    Q_s <- sweep(Q, 2, col_norm2, FUN = "/")
    B_list_s <- lapply(B_list, function(B) B %*% diag(1/col_norm2))
    D_vec <- rep(1, n)
    transform_mode <- FALSE
  } else { # "transform"
    A_s <- A %*% diag(1/col_norm2)
    Q_s <- Q %*% diag(1/col_norm2)
    B_list_s <- lapply(B_list, function(B) B %*% diag(1/col_norm2))
    D_vec <- 1/col_norm2
    transform_mode <- TRUE
  }
  
  # initialize variables (work in x' if transform_mode)
  x <- rep(0, n)
  z <- rep(0, n)
  y <- rep(0, m)
  u <- rep(0, n)
  w <- rep(0, m)
  
  s_list <- vector("list", K)
  y_s_list <- vector("list", K)
  for (kk in 1:K) {
    s_list[[kk]] <- rep(0, 3)
    y_s_list[[kk]] <- rep(0, 3)
  }
  
  # Precompute
  AtA <- crossprod(A_s)
  M_base <- diag(n) + AtA
  
  history <- list()
  
  # weight vector for L2 objective: w_j = D_vec[j]^2 if transform_mode else 1
  if (transform_mode) {
    w_obj <- D_vec^2
  } else {
    w_obj <- rep(1, n)
  }
  
  for (iter in 1:max_iter) {
    # PSD terms
    PSD_term <- matrix(0, n, n)
    rhs_extra <- rep(0, n)
    Bkx <- vector("list", K)
    for (kk in 1:K) {
      Bk <- B_list_s[[kk]]            # 3 x n
      PSD_term <- PSD_term + rho_s * crossprod(Bk)
      rhs_extra <- rhs_extra + rho_s * as.vector(crossprod(Bk, (s_list[[kk]] - y_s_list[[kk]])))
      Bkx[[kk]] <- as.vector(Bk %*% x)
    }
    
    LHS_x <- rho * M_base + PSD_term
    RHS_x <- rho * ((z - u) + crossprod(A_s, (y - w))) + rhs_extra
    
    # Assemble KKT for equality Q_s x = r
    M11 <- LHS_x
    M12 <- t(Q_s)
    M21 <- Q_s
    M22 <- matrix(0, nrow = kq, ncol = kq)
    
    KKT <- rbind(cbind(M11, M12), cbind(M21, M22))
    rhs_kkt <- c(RHS_x, r_vec)
    
    sol <- tryCatch({
      solve(KKT, rhs_kkt)
    }, error = function(e) {
      stop("KKT solve failed: ", e$message)
    })
    x <- as.vector(sol[1:n])
    
    # y-update: projection onto L2 ball ||y|| <= eps
    Ax <- A_s %*% x
    v <- as.vector(Ax + w)
    v_norm <- sqrt(sum(v^2))
    if (v_norm <= eps) {
      y <- v
    } else {
      y <- eps * v / v_norm
    }
    
    # z-update: prox for (1/2) * sum_j w_obj[j] * z_j^2  + (rho/2) ||z - (x+u)||^2
    # closed form: z_j = (rho / (rho + w_obj[j])) * (x_j + u_j)
    z_old <- z
    v_z <- x + u
    denom <- rho + w_obj
    z <- (rho / denom) * v_z
    
    # s_k update: PSD projection
    s_old_list <- s_list
    for (kk in 1:K) {
      b_tilde <- as.vector(Bkx[[kk]] + y_s_list[[kk]])
      Mmat <- matrix(c(b_tilde[1], b_tilde[2], b_tilde[2], b_tilde[3]), nrow = 2, byrow = TRUE)
      ev <- eigen(Mmat, symmetric = TRUE)
      Dpos <- pmax(ev$values, 0)
      Mpos <- ev$vectors %*% diag(Dpos) %*% t(ev$vectors)
      s_list[[kk]] <- c(Mpos[1,1], Mpos[1,2], Mpos[2,2])
    }
    
    # dual updates for PSD constraints
    for (kk in 1:K) {
      y_s_list[[kk]] <- y_s_list[[kk]] + as.vector(Bkx[[kk]] - s_list[[kk]])
    }
    
    # dual updates for z and w
    u <- u + (x - z)
    w <- w + (Ax - y)
    
    # residuals and history
    r_prim_xz <- sqrt(sum((x - z)^2))
    r_prim_Ay <- sqrt(sum((Ax - y)^2))
    r_prim_Bs <- 0
    for (kk in 1:K) {
      r_prim_Bs <- max(r_prim_Bs, sqrt(sum((Bkx[[kk]] - s_list[[kk]])^2)))
    }
    s_dual_z <- rho * sqrt(sum((z - z_old)^2))
    s_dual_s <- 0
    for (kk in 1:K) {
      s_dual_s <- max(s_dual_s, rho_s * sqrt(sum((s_list[[kk]] - s_old_list[[kk]])^2)))
    }
    
    # objective value: (1/2) * sum_j w_obj[j] * z_j^2
    obj_val <- 0.5 * sum(w_obj * (z^2))
    
    history[[iter]] <- list(obj = obj_val,
                            r_xz = r_prim_xz,
                            r_Ay = r_prim_Ay,
                            r_Bs = r_prim_Bs,
                            s_z = s_dual_z,
                            s_s = s_dual_s)
    
    # stopping criteria
    r_prim_max <- max(r_prim_xz, r_prim_Ay, r_prim_Bs)
    s_dual_max <- max(s_dual_z, s_dual_s)
    
    eps_primal <- sqrt(n + m + 3*K) * eps_abs + eps_rel * max(sqrt(sum(x^2)), sqrt(sum(z^2)), sqrt(sum(y^2)))
    eps_dual <- sqrt(n + 3*K) * eps_abs + eps_rel * sqrt(sum(u^2) + sum(unlist(y_s_list)^2))
    
    if (verbose && (iter %% 50 == 0 || iter == 1)) {
      cat(sprintf("iter %d: obj=%.6f, r_prim=%.4e, s_dual=%.4e\n", iter, obj_val, r_prim_max, s_dual_max))
    }
    
    if (r_prim_max <= eps_primal && s_dual_max <= eps_dual) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      break
    }
  } # end for iter
  
  # map back to original x if transform_mode
  if (transform_mode) {
    x_original <- D_vec * x
  } else if (scale_mode == "normalize") {
    x_original <- x / col_norm2
  } else {
    x_original <- x
  }
  
  return(list(x = as.vector(x_original), history = history, iter = iter, D_vec = D_vec, scale_mode = scale_mode))
}



solve_admm_3dB <- function(P, Q, r, B_list, eps, eta, scale_mode = "l2", rho_global = 1.0,
                           rho_blocks = NULL, rho_update_rule = "residual_balance",
                           rho_update_freq = 10, rho_mu = 10.0, rho_tau = 2.0, rho_clip = c(1e-6, 1e6),
                           relax_alpha = 1.0, max_iter = 1000, tol_primal = 1e-4, tol_dual = 1e-4,
                           warm_start = FALSE, cholesky_recompute = "on_rho_change",
                           parallel_workers = 1, eig_batch_size = 128, verbose = TRUE, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)
  
  # ADMM solver variant that assumes each B_list[[i]] is a 3D array D x D x M
  # Inputs:
  #  P: T x M matrix
  #  Q: K x M matrix
  #  r: K-length vector
  #  B_list: list of length N, each element is a 3D array with dim = c(D, D, M)
  #  eps, eta: scalars
  # Hyperparameters: same style as previous implementation (defaults provided)
  
  # Basic dims
  M <- ncol(P)
  if (ncol(Q) != M) stop("P and Q must have same number of columns (M).")
  N <- length(B_list)
  if (N < 1) stop("B_list must be non-empty list of 3D arrays.")
  
  # Validate B_list elements are 3D arrays and consistent
  dims_first <- dim(B_list[[1]])
  if (length(dims_first) != 3) stop("Each B_list[[i]] must be a 3D array D x D x M.")
  D <- dims_first[1]
  if (dims_first[2] != D) stop("Each B_list[[i]] must have shape D x D x M (square first two dims).")
  if (dims_first[3] != M) stop(sprintf("Third dimension of B_list[[1]] must equal M=%d.", M))
  # check all elements have same dims
  for (i in seq_len(N)) {
    di <- dim(B_list[[i]])
    if (length(di) != 3 || di[1] != D || di[2] != D || di[3] != M) {
      stop(sprintf("B_list[[%d]] has inconsistent dimensions: expected %dx%dx%d.", i, D, D, M))
    }
  }
  
  # Convert each 3D array Bi (D x D x M) -> matrix (D^2) x M where column j = vec(B[,,j]) (column-major)
  B_mat_list <- lapply(B_list, function(Bi) {
    D2 <- D * D
    Mat <- matrix(0, nrow = D2, ncol = M)
    for (j in 1:M) Mat[, j] <- as.vector(Bi[ , , j])
    return(Mat)
  })
  # Now B_mat_list[[i]] is (D^2) x M
  
  # --- scaling utilities (apply same scale to P and Q) ---------------------
  scale_apply <- function(A, mode) {
    if (mode == "none") return(list(A = A, meta = list(mode = "none")))
    if (mode == "l2") {
      norms <- apply(A, 2, function(col) sqrt(sum(col^2)))
      norms[norms == 0] <- 1.0
      S <- diag(1 / norms, nrow = length(norms))
      return(list(A = A %*% S, meta = list(mode = "l2", scale = norms)))
    }
    if (mode == "standard") {
      mu <- colMeans(A)
      s <- apply(A, 2, sd)
      s[s == 0] <- 1.0
      S <- diag(1 / s, nrow = length(s))
      A2 <- sweep(A, 2, mu, "-") %*% S
      return(list(A = A2, meta = list(mode = "standard", mean = mu, scale = s)))
    }
    if (mode == "maxabs") {
      mabs <- apply(abs(A), 2, max)
      mabs[mabs == 0] <- 1.0
      S <- diag(1 / mabs, nrow = length(mabs))
      return(list(A = A %*% S, meta = list(mode = "maxabs", scale = mabs)))
    }
    if (mode == "whiten") {
      C <- crossprod(A) / nrow(A)
      ev <- eigen(C, symmetric = TRUE)
      vals <- ev$values
      vals[vals < 1e-12] <- 1e-12
      W <- ev$vectors %*% diag(1 / sqrt(vals)) %*% t(ev$vectors)
      A2 <- A %*% W
      return(list(A = A2, meta = list(mode = "whiten", W = W)))
    }
    stop("Unknown scale_mode")
  }
  
  sp <- scale_apply(P, scale_mode); P_s <- sp$A; scale_meta_P <- sp$meta
  sq <- scale_apply(Q, scale_mode); Q_s <- sq$A; scale_meta_Q <- sq$meta
  # For simplicity we assume same scale_mode applied to P and Q; if different, user must pre-transform r and B_list.
  
  # If diagonal scaling was used, transform B_mat_list accordingly: B_new = B_old %*% S_inv
  B_list_s <- B_mat_list
  if (!is.null(scale_meta_P$scale)) {
    Sdiag <- scale_meta_P$scale
    Sinv <- diag(1 / Sdiag, nrow = length(Sdiag))
    B_list_s <- lapply(B_mat_list, function(Bi) Bi %*% Sinv)
  } else if (!is.null(scale_meta_P$W)) {
    Winv <- solve(scale_meta_P$W)
    B_list_s <- lapply(B_mat_list, function(Bi) Bi %*% Winv)
  }
  
  # --- initialize variables ------------------------------------------------
  x <- rep(0, M); y <- rep(0, M)
  u <- rep(0, nrow(P_s)); v <- rep(0, nrow(Q_s))
  S_list <- lapply(1:N, function(i) matrix(0, nrow = D, ncol = D))
  lam_y <- rep(0, M); lam_P <- rep(0, nrow(P_s)); lam_Q <- rep(0, nrow(Q_s))
  lam_G <- lapply(1:N, function(i) rep(0, D * D))
  
  if (is.null(rho_blocks)) {
    rho_y <- rho_global; rho_P <- rho_global; rho_Q <- rho_global; rho_G <- rho_global
  } else {
    rho_y <- ifelse(!is.null(rho_blocks$y), rho_blocks$y, rho_global)
    rho_P <- ifelse(!is.null(rho_blocks$P), rho_blocks$P, rho_global)
    rho_Q <- ifelse(!is.null(rho_blocks$Q), rho_blocks$Q, rho_global)
    rho_G <- ifelse(!is.null(rho_blocks$G), rho_blocks$G, rho_global)
  }
  
  # Precompute PtP, QtQ, sum BtB
  PtP <- crossprod(P_s); QtQ <- crossprod(Q_s)
  sumBtB <- matrix(0, nrow = M, ncol = M)
  for (i in 1:N) sumBtB <- sumBtB + crossprod(B_list_s[[i]])
  
  soft_threshold <- function(z, kappa) { pmax(0, z - kappa) - pmax(0, -z - kappa) }
  project_to_PSD <- function(A) {
    A_sym <- (A + t(A)) / 2
    ev <- eigen(A_sym, symmetric = TRUE)
    vals <- pmax(ev$values, 0)
    ev$vectors %*% diag(vals) %*% t(ev$vectors)
  }
  
  build_A <- function(rho_y, rho_P, rho_Q, rho_G) {
    A <- rho_y * diag(M) + rho_P * PtP + rho_Q * QtQ + rho_G * sumBtB
    A <- (A + t(A)) / 2 + 1e-12 * diag(M)
    return(A)
  }
  
  A <- build_A(rho_y, rho_P, rho_Q, rho_G)
  cholA <- tryCatch(chol(A), error = function(e) {
    A <- A + 1e-8 * diag(M); chol(A)
  })
  
  Pt <- t(P_s); Qt <- t(Q_s)
  
  stats <- list(primal_res = numeric(), dual_res = numeric(), obj = numeric(), rho_hist = list())
  rho_history <- list(list(rho_y = rho_y, rho_P = rho_P, rho_Q = rho_Q, rho_G = rho_G))
  
  # Main ADMM loop
  for (k in 1:max_iter) {
    # x-update
    b <- rho_y * (y - lam_y / rho_y) +
      rho_P * (Pt %*% (u - lam_P / rho_P)) +
      rho_Q * (Qt %*% (v - lam_Q / rho_Q))
    tmpG <- rep(0, M)
    for (i in 1:N) {
      Bi <- B_list_s[[i]]
      si_vec <- as.vector(S_list[[i]])
      tmpG <- tmpG + rho_G * (crossprod(Bi, (si_vec - lam_G[[i]] / rho_G)))
    }
    b <- b + tmpG
    x_new <- backsolve(cholA, forwardsolve(t(cholA), b))
    
    # y-update
    y_new <- soft_threshold(x_new + lam_y / rho_y, 1 / rho_y)
    
    # u-update
    z_u <- P_s %*% x_new + lam_P / rho_P
    norm_z_u <- sqrt(sum(z_u^2))
    u_new <- if (norm_z_u <= eps) z_u else (eps / norm_z_u) * z_u
    
    # v-update
    z_v <- Q_s %*% x_new + lam_Q / rho_Q
    diff_v <- z_v - r
    norm_diff_v <- sqrt(sum(diff_v^2))
    v_new <- if (norm_diff_v <= eta) z_v else r + (eta / norm_diff_v) * diff_v
    
    # S_i updates (PSD projection) - process in batches if desired
    S_list_new <- vector("list", N)
    batch_size <- min(eig_batch_size, N)
    for (start in seq(1, N, by = batch_size)) {
      end <- min(start + batch_size - 1, N)
      for (i in start:end) {
        Bi <- B_list_s[[i]]
        vec_target <- Bi %*% x_new + lam_G[[i]] / rho_G
        mat_target <- matrix(vec_target, nrow = D, ncol = D)
        S_list_new[[i]] <- project_to_PSD(mat_target)
      }
    }
    
    # residuals and dual updates
    r_y <- x_new - y_new
    r_P <- P_s %*% x_new - u_new
    r_Q <- Q_s %*% x_new - v_new
    r_G_vecs <- numeric()
    for (i in 1:N) r_G_vecs <- c(r_G_vecs, as.vector(B_list_s[[i]] %*% x_new - as.vector(S_list_new[[i]])))
    
    s_y <- rho_y * (y_new - y)
    s_P <- rho_P * (t(P_s) %*% (u_new - u))
    s_Q <- rho_Q * (t(Q_s) %*% (v_new - v))
    s_G_vec <- numeric()
    for (i in 1:N) s_G_vec <- c(s_G_vec, rho_G * (t(B_list_s[[i]]) %*% (as.vector(S_list_new[[i]]) - as.vector(S_list[[i]]))))
    
    lam_y <- lam_y + rho_y * r_y
    lam_P <- lam_P + rho_P * r_P
    lam_Q <- lam_Q + rho_Q * r_Q
    for (i in 1:N) lam_G[[i]] <- lam_G[[i]] + rho_G * (as.vector(B_list_s[[i]] %*% x_new - as.vector(S_list_new[[i]])))
    
    norm_primal <- sqrt(sum(r_y^2) + sum(r_P^2) + sum(r_Q^2) + sum(r_G_vecs^2))
    norm_dual <- sqrt(sum(s_y^2) + sum(s_P^2) + sum(s_Q^2) + sum(s_G_vec^2))
    obj_val <- sum(abs(x_new))
    
    stats$primal_res <- c(stats$primal_res, norm_primal)
    stats$dual_res <- c(stats$dual_res, norm_dual)
    stats$obj <- c(stats$obj, obj_val)
    stats$rho_hist[[length(stats$rho_hist) + 1]] <- list(rho_y = rho_y, rho_P = rho_P, rho_Q = rho_Q, rho_G = rho_G)
    
    if (verbose && (k <= 5 || k %% max(1, floor(max_iter/10)) == 0)) {
      cat(sprintf("iter %4d: obj=%.6g, primal=%.3e, dual=%.3e\n", k, obj_val, norm_primal, norm_dual))
    }
    
    if (norm_primal < tol_primal && norm_dual < tol_dual) {
      x <- x_new; y <- y_new; u <- u_new; v <- v_new; S_list <- S_list_new
      if (verbose) cat("Converged at iter", k, "\n")
      break
    }
    
    # rho auto-update
    if (rho_update_rule == "residual_balance" && (k %% rho_update_freq == 0)) {
      if (norm_primal > rho_mu * norm_dual) {
        rho_y_old <- rho_y; rho_P_old <- rho_P; rho_Q_old <- rho_Q; rho_G_old <- rho_G
        rho_y <- min(rho_y * rho_tau, rho_clip[2])
        rho_P <- min(rho_P * rho_tau, rho_clip[2])
        rho_Q <- min(rho_Q * rho_tau, rho_clip[2])
        rho_G <- min(rho_G * rho_tau, rho_clip[2])
        lam_y <- lam_y * (rho_y_old / rho_y)
        lam_P <- lam_P * (rho_P_old / rho_P)
        lam_Q <- lam_Q * (rho_Q_old / rho_Q)
        for (i in 1:N) lam_G[[i]] <- lam_G[[i]] * (rho_G_old / rho_G)
        if (cholesky_recompute == "on_rho_change") {
          A <- build_A(rho_y, rho_P, rho_Q, rho_G); cholA <- chol(A)
        }
      } else if (norm_dual > rho_mu * norm_primal) {
        rho_y_old <- rho_y; rho_P_old <- rho_P; rho_Q_old <- rho_Q; rho_G_old <- rho_G
        rho_y <- max(rho_y / rho_tau, rho_clip[1])
        rho_P <- max(rho_P / rho_tau, rho_clip[1])
        rho_Q <- max(rho_Q / rho_tau, rho_clip[1])
        rho_G <- max(rho_G / rho_tau, rho_clip[1])
        lam_y <- lam_y * (rho_y_old / rho_y)
        lam_P <- lam_P * (rho_P_old / rho_P)
        lam_Q <- lam_Q * (rho_Q_old / rho_Q)
        for (i in 1:N) lam_G[[i]] <- lam_G[[i]] * (rho_G_old / rho_G)
        if (cholesky_recompute == "on_rho_change") {
          A <- build_A(rho_y, rho_P, rho_Q, rho_G); cholA <- chol(A)
        }
      }
      rho_history[[length(rho_history) + 1]] <- list(iter = k, rho_y = rho_y, rho_P = rho_P, rho_Q = rho_Q, rho_G = rho_G)
    }
    
    # update for next iter
    x <- x_new; y <- y_new; u <- u_new; v <- v_new; S_list <- S_list_new
  } # end loop
  
  # Recover original-scale x if scaling applied
  x_scaled <- x
  x_orig <- x_scaled
  if (!is.null(scale_meta_P$W)) {
    x_orig <- scale_meta_P$W %*% x_scaled
  } else if (!is.null(scale_meta_P$scale)) {
    x_orig <- diag(scale_meta_P$scale) %*% x_scaled
  }
  
  result <- list(
    x = as.numeric(x_orig),
    x_scaled = as.numeric(x_scaled),
    stats = list(iters = length(stats$obj), primal_res = stats$primal_res, dual_res = stats$dual_res, obj = stats$obj, rho_hist = stats$rho_hist),
    meta = list(scale_meta = scale_meta_P, rho_history = rho_history)
  )
  return(result)
}


solve_admm_01 <- function(P, Q, r, B_list, eps, eta,
                       scale_mode = "none",    # "none","l2","standard","maxabs","operator_norm","whiten"
                       rho_global = 1.0, rho_blocks = NULL,    # list with names "y","P","Q","G" or NULL
                       rho_update_rule = "residual_balance", # "none","residual_balance"
                       rho_update_freq = 10, rho_mu = 10.0, rho_tau = 2.0, rho_clip = c(1e-6, 1e6),
                       relax_alpha = 1.0,    # over-relaxation
                       max_iter = 1000, tol_primal = 1e-4, tol_dual = 1e-4,
                       warm_start = FALSE, cholesky_recompute = "on_rho_change", # "always","on_rho_change","never"
                       parallel_workers = 1, eig_batch_size = 128, verbose = TRUE) {
  
  # ADMM solver for:
  # minimize ||x||_1
  # s.t. ||P x||_2 <= eps, ||Q x - r||_2 <= eta, G^(i)(x) PSD for all i
  #
  # Inputs:
  #  P: T x M matrix
  #  Q: K x M matrix
  #  r: K-length vector
  #  B_list: list of length N, each is (D^2) x M matrix mapping x -> vec(G^(i))
  #  eps, eta: scalars
  #
  # Hyperparameters are passed as named arguments (see defaults).
  #
  
  
  t_start <- proc.time()[3]
  
  # Dimensions
  M <- ncol(P)
  if (ncol(Q) != M) stop("P and Q must have same number of columns (M).")
  K <- nrow(Q)
  N <- length(B_list)
  if (N == 0) stop("B_list must contain at least one B_i.")
  D2 <- nrow(B_list[[1]])
  D <- round(sqrt(D2))
  if (D * D != D2) stop("Each B_i must have D^2 rows (square matrix when reshaped).")
  
  # --- Scaling utilities ---------------------------------------------------
  scale_meta <- list()
  scale_apply <- function(A, mode) {
    if (mode == "none") return(list(A=A, meta=list(mode="none")))
    if (mode == "l2") {
      norms <- apply(A, 2, function(col) sqrt(sum(col^2)))
      norms[norms == 0] <- 1.0
      S <- diag(1 / norms, nrow = length(norms))
      return(list(A = A %*% S, meta = list(mode="l2", scale = norms)))
    }
    if (mode == "standard") {
      mu <- colMeans(A)
      s <- apply(A, 2, sd)
      s[s == 0] <- 1.0
      S <- diag(1 / s, nrow = length(s))
      A2 <- sweep(A, 2, mu, "-") %*% S
      return(list(A = A2, meta = list(mode="standard", mean = mu, scale = s)))
    }
    if (mode == "maxabs") {
      mabs <- apply(abs(A), 2, max)
      mabs[mabs == 0] <- 1.0
      S <- diag(1 / mabs, nrow = length(mabs))
      return(list(A = A %*% S, meta = list(mode="maxabs", scale = mabs)))
    }
    if (mode == "operator_norm") {
      # approximate operator norm by largest singular value via svd (costly for large T)
      svals <- apply(A, 2, function(col) sqrt(sum(col^2))) # fallback to column norms
      svals[svals == 0] <- 1.0
      S <- diag(1 / svals, nrow = length(svals))
      return(list(A = A %*% S, meta = list(mode="operator_norm", scale = svals)))
    }
    if (mode == "whiten") {
      # PCA whitening on columns: compute covariance of columns (M x M)
      C <- crossprod(A) / nrow(A)
      ev <- eigen(C, symmetric = TRUE)
      vals <- ev$values
      vecs <- ev$vectors
      # regularize tiny eigenvalues
      vals[vals < 1e-12] <- 1e-12
      W <- vecs %*% diag(1 / sqrt(vals)) %*% t(vecs)
      A2 <- A %*% W
      return(list(A = A2, meta = list(mode="whiten", W = W)))
    }
    stop("Unknown scale_mode")
  }
  
  # Apply scaling to P and Q (keep meta for inverse transform)
  sp <- scale_apply(P, scale_mode); P_s <- sp$A; scale_meta$P <- sp$meta
  sq <- scale_apply(Q, scale_mode); Q_s <- sq$A; scale_meta$Q <- sq$meta
  
  # If whiten, we must also transform B_list accordingly: x_new = W^{-1} x_old
  # We will keep track of transform on x: x_s = Sx * x_orig + shift (shift only for standard)
  x_transform <- list() # to recover original x
  if (scale_meta$P$mode == "whiten") {
    W <- scale_meta$P$W
    # For whiten, we applied same transform to P and Q; assume same W
    # Transform B_i: B_i_new corresponds to mapping x_new -> vec(G) where x_new = W^{-1} x_orig
    Winv <- solve(W)
    B_list_s <- lapply(B_list, function(Bi) Bi %*% Winv)
    x_transform$W <- W
  } else if (scale_meta$P$mode %in% c("l2","standard","maxabs","operator_norm")) {
    # These are diagonal scalings S such that A_s = A %*% S
    if (!is.null(scale_meta$P$scale)) {
      Sdiag <- diag(scale_meta$P$scale, nrow = length(scale_meta$P$scale))
      # A_s = A %*% diag(1/scale) in our implementation; we stored scale as original norms
      # For B: B_i_new = B_i %*% S^{-1} where S^{-1} = diag(1/scale)
      Sinv <- diag(1 / scale_meta$P$scale, nrow = length(scale_meta$P$scale))
      B_list_s <- lapply(B_list, function(Bi) Bi %*% Sinv)
      x_transform$Sdiag <- scale_meta$P$scale
    } else {
      B_list_s <- B_list
    }
  } else {
    B_list_s <- B_list
  }
  
  # If scale_mode == "standard", note mean shift for Q*r consistency (we used column centering)
  # For simplicity, we assume r is in same scale as Q; if standardization used, we must transform r:
  if (scale_meta$Q$mode == "standard") {
    muQ <- scale_meta$Q$mean
    sQ <- scale_meta$Q$scale
    # Q_s = (Q - muQ) %*% diag(1/sQ)
    # So Q x - r  -> Q_s x_s - r_s where r_s = (r - muQ %*% x?)  -- complicated if x scaling differs.
    # To avoid complexity, we assume standardization is applied identically to P and Q (common scale_mode).
    # If user uses different scaling for P and Q, they must handle r externally.
  }
  
  # --- initialize ADMM variables -------------------------------------------
  # Variables: x (M), y (M), u (T), v (K), S_i (D x D matrices), and duals
  x <- rep(0, M)
  y <- rep(0, M)
  u <- rep(0, nrow(P_s))
  v <- rep(0, nrow(Q_s))
  S_list <- lapply(1:N, function(i) matrix(0, nrow = D, ncol = D))
  
  # duals
  lam_y <- rep(0, M)
  lam_P <- rep(0, nrow(P_s))
  lam_Q <- rep(0, nrow(Q_s))
  lam_G <- lapply(1:N, function(i) rep(0, D2))
  
  # rho handling: either global or block-specific
  if (is.null(rho_blocks)) {
    rho_y <- rho_global
    rho_P <- rho_global
    rho_Q <- rho_global
    rho_G <- rho_global
  } else {
    rho_y <- ifelse(!is.null(rho_blocks$y), rho_blocks$y, rho_global)
    rho_P <- ifelse(!is.null(rho_blocks$P), rho_blocks$P, rho_global)
    rho_Q <- ifelse(!is.null(rho_blocks$Q), rho_blocks$Q, rho_global)
    rho_G <- ifelse(!is.null(rho_blocks$G), rho_blocks$G, rho_global)
  }
  
  rho_history <- list()
  rho_history[[1]] <- list(rho_y = rho_y, rho_P = rho_P, rho_Q = rho_Q, rho_G = rho_G)
  
  # Precompute P^T P, Q^T Q, sum B_i^T B_i
  PtP <- crossprod(P_s)   # M x M
  QtQ <- crossprod(Q_s)   # M x M
  sumBtB <- matrix(0, nrow = M, ncol = M)
  for (i in 1:N) {
    Bi <- B_list_s[[i]]
    sumBtB <- sumBtB + crossprod(Bi)
  }
  
  # Helper: soft-threshold
  soft_threshold <- function(z, kappa) {
    pmax(0, z - kappa) - pmax(0, -z - kappa)
  }
  
  # Helper: PSD projection of D x D matrix
  project_to_PSD <- function(A) {
    # ensure symmetry
    A_sym <- (A + t(A)) / 2
    ev <- eigen(A_sym, symmetric = TRUE)
    vals <- ev$values
    vecs <- ev$vectors
    vals[vals < 0] <- 0
    return(vecs %*% diag(vals) %*% t(vecs))
  }
  
  # Build initial A matrix for x-update: A = rho_y I + rho_P PtP + rho_Q QtQ + rho_G sumBtB
  build_A <- function(rho_y, rho_P, rho_Q, rho_G) {
    A <- rho_y * diag(M) + rho_P * PtP + rho_Q * QtQ + rho_G * sumBtB
    # regularize slightly for numerical stability
    A <- (A + t(A)) / 2 + 1e-12 * diag(M)
    return(A)
  }
  
  A <- build_A(rho_y, rho_P, rho_Q, rho_G)
  cholA <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(cholA)) {
    # fallback to adding small diagonal
    A <- A + 1e-8 * diag(M)
    cholA <- chol(A)
  }
  
  # Precompute some constants
  Pt <- t(P_s); Qt <- t(Q_s)
  
  # Stats storage
  stats <- list(primal_res = numeric(), dual_res = numeric(), obj = numeric(), rho_hist = list())
  
  # ADMM main loop
  for (k in 1:max_iter) {
    # --- x-update: solve quadratic minimization
    # RHS b = rho_y*(y - lam_y/rho_y) + rho_P * P^T*(u - lam_P/rho_P) + rho_Q * Q^T*(v - lam_Q/rho_Q)
    #       + rho_G * sum_i B_i^T*(vec(S_i) - lam_G_i/rho_G)
    b <- rho_y * (y - lam_y / rho_y) +
      rho_P * (Pt %*% (u - lam_P / rho_P)) +
      rho_Q * (Qt %*% (v - lam_Q / rho_Q))
    # add G contributions
    tmpG <- rep(0, M)
    for (i in 1:N) {
      Bi <- B_list_s[[i]]
      si_vec <- as.vector(S_list[[i]])
      tmpG <- tmpG + rho_G * (crossprod(Bi, (si_vec - lam_G[[i]] / rho_G)))
    }
    b <- b + tmpG
    
    # solve A x = b using Cholesky
    # if rho changed and cholesky_recompute requires, recompute
    x_new <- backsolve(cholA, forwardsolve(t(cholA), b))
    
    # --- y-update: soft-thresholding (prox of l1)
    z_y <- x_new + lam_y / rho_y
    y_new <- soft_threshold(z_y, 1 / rho_y)
    
    # --- u-update: projection onto l2 ball ||u||_2 <= eps
    z_u <- P_s %*% x_new + lam_P / rho_P
    norm_z_u <- sqrt(sum(z_u^2))
    if (norm_z_u <= eps) {
      u_new <- z_u
    } else {
      u_new <- (eps / norm_z_u) * z_u
    }
    
    # --- v-update: projection onto ball centered at r with radius eta
    z_v <- Q_s %*% x_new + lam_Q / rho_Q
    diff_v <- z_v - r
    norm_diff_v <- sqrt(sum(diff_v^2))
    if (norm_diff_v <= eta) {
      v_new <- z_v
    } else {
      v_new <- r + (eta / norm_diff_v) * diff_v
    }
    
    # --- S_i updates: PSD projection
    S_list_new <- vector("list", N)
    for (i in 1:N) {
      Bi <- B_list_s[[i]]
      vec_target <- Bi %*% x_new + lam_G[[i]] / rho_G
      mat_target <- matrix(vec_target, nrow = D, ncol = D)
      Sproj <- project_to_PSD(mat_target)
      S_list_new[[i]] <- Sproj
    }
    
    # --- dual updates (with optional over-relaxation)
    # primal residuals
    r_y <- x_new - y_new
    r_P <- P_s %*% x_new - u_new
    r_Q <- Q_s %*% x_new - v_new
    r_G_vecs <- numeric()
    for (i in 1:N) {
      rGi <- as.vector(B_list_s[[i]] %*% x_new - as.vector(S_list_new[[i]]))
      r_G_vecs <- c(r_G_vecs, rGi)
    }
    
    # dual residuals (scaled)
    s_y <- rho_y * (y_new - y)
    s_P <- rho_P * (t(P_s) %*% (u_new - u))
    s_Q <- rho_Q * (t(Q_s) %*% (v_new - v))
    s_G_vec <- numeric()
    for (i in 1:N) {
      sGi <- rho_G * (t(B_list_s[[i]]) %*% (as.vector(S_list_new[[i]]) - as.vector(S_list[[i]])))
      s_G_vec <- c(s_G_vec, sGi)
    }
    
    # update dual variables
    lam_y <- lam_y + rho_y * r_y
    lam_P <- lam_P + rho_P * r_P
    lam_Q <- lam_Q + rho_Q * r_Q
    for (i in 1:N) {
      lam_G[[i]] <- lam_G[[i]] + rho_G * (as.vector(B_list_s[[i]] %*% x_new - as.vector(S_list_new[[i]])))
    }
    
    # compute norms for stopping and rho update
    norm_primal <- sqrt(sum(r_y^2) + sum(r_P^2) + sum(r_Q^2) + sum(r_G_vecs^2))
    norm_dual <- sqrt(sum(s_y^2) + sum(s_P^2) + sum(s_Q^2) + sum(s_G_vec^2))
    
    # objective (l1 norm)
    obj_val <- sum(abs(x_new))
    
    # store stats
    stats$primal_res <- c(stats$primal_res, norm_primal)
    stats$dual_res <- c(stats$dual_res, norm_dual)
    stats$obj <- c(stats$obj, obj_val)
    stats$rho_hist[[length(stats$rho_hist) + 1]] <- list(rho_y = rho_y, rho_P = rho_P, rho_Q = rho_Q, rho_G = rho_G)
    
    if (verbose && (k %% max(1, floor(max_iter/10)) == 0 || k <= 5)) {
      cat(sprintf("iter %4d: obj=%.6g, primal=%.3e, dual=%.3e\n", k, obj_val, norm_primal, norm_dual))
    }
    
    # check stopping
    if (norm_primal < tol_primal && norm_dual < tol_dual) {
      x <- x_new
      y <- y_new
      u <- u_new
      v <- v_new
      S_list <- S_list_new
      if (verbose) cat("Converged at iter", k, "\n")
      break
    }
    
    # rho auto-update (residual balancing) every rho_update_freq iterations
    if (rho_update_rule == "residual_balance" && (k %% rho_update_freq == 0)) {
      # compare primal and dual residuals
      if (norm_primal > rho_mu * norm_dual) {
        # increase rho
        rho_y_old <- rho_y; rho_P_old <- rho_P; rho_Q_old <- rho_Q; rho_G_old <- rho_G
        rho_y <- min(rho_y * rho_tau, rho_clip[2])
        rho_P <- min(rho_P * rho_tau, rho_clip[2])
        rho_Q <- min(rho_Q * rho_tau, rho_clip[2])
        rho_G <- min(rho_G * rho_tau, rho_clip[2])
        # scale duals
        lam_y <- lam_y * (rho_y_old / rho_y)
        lam_P <- lam_P * (rho_P_old / rho_P)
        lam_Q <- lam_Q * (rho_Q_old / rho_Q)
        for (i in 1:N) lam_G[[i]] <- lam_G[[i]] * (rho_G_old / rho_G)
        # rebuild A and chol if needed
        if (cholesky_recompute == "on_rho_change") {
          A <- build_A(rho_y, rho_P, rho_Q, rho_G)
          cholA <- chol(A)
        }
      } else if (norm_dual > rho_mu * norm_primal) {
        # decrease rho
        rho_y_old <- rho_y; rho_P_old <- rho_P; rho_Q_old <- rho_Q; rho_G_old <- rho_G
        rho_y <- max(rho_y / rho_tau, rho_clip[1])
        rho_P <- max(rho_P / rho_tau, rho_clip[1])
        rho_Q <- max(rho_Q / rho_tau, rho_clip[1])
        rho_G <- max(rho_G / rho_tau, rho_clip[1])
        lam_y <- lam_y * (rho_y_old / rho_y)
        lam_P <- lam_P * (rho_P_old / rho_P)
        lam_Q <- lam_Q * (rho_Q_old / rho_Q)
        for (i in 1:N) lam_G[[i]] <- lam_G[[i]] * (rho_G_old / rho_G)
        if (cholesky_recompute == "on_rho_change") {
          A <- build_A(rho_y, rho_P, rho_Q, rho_G)
          cholA <- chol(A)
        }
      }
      rho_history[[length(rho_history) + 1]] <- list(iter = k, rho_y = rho_y, rho_P = rho_P, rho_Q = rho_Q, rho_G = rho_G)
    }
    
    # update variables for next iter
    x <- x_new; y <- y_new; u <- u_new; v <- v_new; S_list <- S_list_new
    
    # optionally recompute chol if policy is always
    if (cholesky_recompute == "always" && (k %% rho_update_freq == 0)) {
      A <- build_A(rho_y, rho_P, rho_Q, rho_G)
      cholA <- chol(A)
    }
  } # end for
  
  t_end <- proc.time()[3]
  elapsed <- t_end - t_start
  
  # Recover x to original scale if scaling applied
  x_orig <- x
  if (!is.null(x_transform$W)) {
    # x_s = W^{-1} x_orig => original x = W^{-1}^{-1} x_s = W x_s
    x_orig <- x_transform$W %*% x
  } else if (!is.null(x_transform$Sdiag)) {
    # we used Sdiag as original column norms; recall we scaled A by S^{-1} so x_s = S^{-1} x_orig
    # thus x_orig = S * x_s
    x_orig <- diag(x_transform$Sdiag) %*% x
  }
  
  result <- list(
    x = as.numeric(x_orig),
    x_scaled = as.numeric(x),
    stats = list(iters = length(stats$obj),
                 primal_res = stats$primal_res,
                 dual_res = stats$dual_res,
                 obj = stats$obj,
                 rho_hist = stats$rho_hist,
                 elapsed = elapsed),
    meta = list(scale_meta = scale_meta, rho_history = rho_history, cholesky_recompute = cholesky_recompute)
  )
  
  return(result)
}


#遅すぎる???
admm_sdp01 <- function(P, Q, r, B_list, rho = 1.0, alpha = 1.0, max_iter = 1000, tol = 1e-4, verbose = TRUE) {

  # ADMM solver for: min x1 s.t. P x = 0, Q x = r, G^i(x) PSD
  # B_list: list of arrays dim = c(D, D, M), B_list[[i]][k,l,] is length-M vector

  # Basic sizes
  M <- ncol(P)
  K <- nrow(Q)
  N <- length(B_list)
  D <- dim(B_list[[1]])[1]
  if (dim(B_list[[1]])[2] != D) stop("B arrays must be D x D x M")
  nArows <- N * D * D

  # Build A matrix: each block A_i is (D^2) x M, stack vertically
  A <- matrix(0, nrow = nArows, ncol = M)
  row_idx <- 1
  for (i in seq_len(N)) {
    Bi <- B_list[[i]]  # D x D x M
    for (k in seq_len(D)) {
      for (l in seq_len(D)) {
        A[row_idx, ] <- Bi[k, l, ]  # length M
        row_idx <- row_idx + 1
      }
    }
  }
  # Helper to reshape vec -> D x D for block i
  vec_to_mat <- function(v, i) {
    start <- (i-1)*D*D + 1
    m <- matrix(v[start:(start + D*D - 1)], nrow = D, ncol = D, byrow = TRUE)
    return(m)
  }
  mat_to_vec <- function(mat, i) {
    as.vector(t(mat))  # row-major to match construction
  }

  # Initializations
  x <- rep(0, M)
  # Make initial x feasible for P x = 0, Q x = r by projection
  project_affine <- function(x0) {
    # Solve min ||x - x0||^2 s.t. P x = 0, Q x = r
    # KKT: [I  P^T  Q^T; P 0 0; Q 0 0] [x; lam; mu] = [x0; 0; r]
    Aeq <- rbind(P, Q)
    meq <- nrow(Aeq)
    KKT <- rbind(
      cbind(diag(M), t(Aeq)),
      cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
    )
    rhs <- c(x0, rep(0, nrow(P)), r)
    sol <- tryCatch(solve(KKT, rhs), error = function(e) NULL)
    if (is.null(sol)) stop("KKT solve failed in projection; check P,Q ranks")
    return(sol[1:M])
  }
  x <- project_affine(x)

  # Duals and local variables
  Y_vec <- rep(0, nArows)      # stacked vec(Y^i)
  U <- rep(0, nArows)          # scaled duals (u = lambda / rho)

  # Precompute A^T A and A^T
  At <- t(A)
  AtA <- At %*% A

  # Precompute KKT matrix structure for x-update:
  # Solve: min c^T x + (rho/2)||A x - b||^2 s.t. P x = 0, Q x = r
  # => (rho AtA) x + P^T lam + Q^T mu = -c + rho At b
  # KKT matrix:
  Aeq <- rbind(P, Q)
  meq <- nrow(Aeq)
  # We'll build KKT each iteration because RHS changes; but LHS (rho AtA, Aeq) can be reused
  KKT_LHS_top <- rho * AtA
  KKT_LHS <- rbind(
    cbind(KKT_LHS_top, t(Aeq)),
    cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
  )

  # objective linear term c (minimize x1)
  cvec <- rep(0, M); cvec[1] <- 1.0

  # History
  hist <- list(primal = numeric(), dual = numeric())

  # ADMM loop
  for (k in seq_len(max_iter)) {
    # x-update: compute b = vec(Y) - U
    b_vec <- Y_vec - U
    Atb <- At %*% b_vec
    rhs_top <- -cvec + rho * Atb
    rhs <- c(rhs_top, rep(0, nrow(P)), r)
    sol <- tryCatch(solve(KKT_LHS, rhs), error = function(e) NULL)
    if (is.null(sol)) stop("KKT solve failed in x-update; consider regularization")
    x_old <- x
    x <- sol[1:M]

    # Over-relaxation for A x
    Ax <- A %*% x
    Ax_hat <- alpha * Ax + (1 - alpha) * (Y_vec)  # note: Y_vec is previous Y

    # Y-update: for each block, form S = reshape(Ax_hat + U) and project to PSD
    primal_res_norm_sq <- 0
    for (i in seq_len(N)) {
      start <- (i-1)*D*D + 1
      end <- start + D*D - 1
      s_vec <- Ax_hat[start:end] + U[start:end]
      S <- matrix(s_vec, nrow = D, ncol = D, byrow = TRUE)
      # Symmetrize
      S_sym <- (S + t(S)) / 2
      # Eigen-decompose and project
      ev <- eigen(S_sym, symmetric = TRUE)
      vals <- ev$values
      vecs <- ev$vectors
      vals_proj <- pmax(vals, 0)
      Ymat <- vecs %*% diag(vals_proj) %*% t(vecs)
      Yv <- as.vector(t(Ymat))
      Y_vec[start:end] <- Yv
      # accumulate primal residual norm
      primal_res_norm_sq <- primal_res_norm_sq + sum((Ax[start:end] - Yv)^2)
    }
    primal_res <- sqrt(primal_res_norm_sq)

    # Dual update: U <- U + Ax_hat - Y
    dual_res_norm_sq <- 0
    U_old <- U
    U <- U + (Ax_hat - Y_vec)
    # dual residual: rho * A^T (Y - Y_prev)
    dual_vec <- rho * (At %*% (Y_vec - (Y_vec - (Ax_hat - U + (Ax_hat - U)))))
    # The above is messy; simpler compute dual residual as rho * A^T (Y - Y_prev)
    # But we didn't store Y_prev separately; compute Y_prev = Y_vec - (Ax_hat - U)
    # Recompute properly:
    # Instead, store Y_prev at start of Y-update. Let's adjust: (we'll compute dual residual approx)
    # For simplicity, compute dual residual as rho * norm(At %*% (Y_vec - (Ax_hat - U)))
    # But to avoid confusion, compute dual residual via change in Y_vec:
    # (we have U_old and U, but we need Y change)
    # We'll approximate dual residual using rho * ||At %*% (Y_vec - Y_prev)||.
    # To do that, we need Y_prev; compute Y_prev = Y_vec - (Y_vec - previous Y) impossible now.
    # Simpler: compute dual residual via rho * ||At %*% (Y_vec - Y_prev)|| using stored Y_prev.
    # For correctness, we store Y_prev before update. (Refactor: move Y_prev earlier)

    # For clarity and correctness, redo Y-update block with Y_prev stored:
    # (To keep code readable, we implement a correct version below and break loop to re-run)
    break
  }

  # --- Re-implement loop correctly with Y_prev tracking ---
  # Reinitialize variables and history
  x <- project_affine(rep(0, M))
  Y_vec <- rep(0, nArows)
  U <- rep(0, nArows)
  hist$primal <- numeric(0); hist$dual <- numeric(0)

  for (k in seq_len(max_iter)) {
    # x-update
    b_vec <- Y_vec - U
    Atb <- At %*% b_vec
    rhs_top <- -cvec + rho * Atb
    rhs <- c(rhs_top, rep(0, nrow(P)), r)
    sol <- tryCatch(solve(KKT_LHS, rhs), error = function(e) NULL)
    if (is.null(sol)) stop("KKT solve failed in x-update; consider regularization")
    x_old <- x
    x <- sol[1:M]

    # compute Ax and over-relaxed version
    Ax <- A %*% x
    Ax_hat <- alpha * Ax + (1 - alpha) * (Y_vec)

    # Y-update with Y_prev
    Y_prev <- Y_vec
    primal_res_norm_sq <- 0
    for (i in seq_len(N)) {
      start <- (i-1)*D*D + 1
      end <- start + D*D - 1
      s_vec <- Ax_hat[start:end] + U[start:end]
      S <- matrix(s_vec, nrow = D, ncol = D, byrow = TRUE)
      S_sym <- (S + t(S)) / 2
      ev <- eigen(S_sym, symmetric = TRUE)
      vals_proj <- pmax(ev$values, 0)
      Ymat <- ev$vectors %*% diag(vals_proj) %*% t(ev$vectors)
      Yv <- as.vector(t(Ymat))
      Y_vec[start:end] <- Yv
      primal_res_norm_sq <- primal_res_norm_sq + sum((Ax[start:end] - Yv)^2)
    }
    primal_res <- sqrt(primal_res_norm_sq)

    # Dual update
    U <- U + (Ax_hat - Y_vec)
    # dual residual: rho * ||A^T (Y - Y_prev)||
    dual_res <- sqrt(sum((rho * (At %*% (Y_vec - Y_prev)))^2))

    hist$primal <- c(hist$primal, primal_res)
    hist$dual <- c(hist$dual, dual_res)

    if (verbose && (k %% 10 == 0 || k == 1)) {
      cat(sprintf("iter %4d: primal = %.3e, dual = %.3e\n", k, primal_res, dual_res))
    }

    # stopping
    if (primal_res < tol && dual_res < tol) {
      if (verbose) cat("Converged at iter", k, "\n")
      return(list(x = x, history = hist, info = list(iter = k, primal = primal_res, dual = dual_res)))
    }

    # adaptive rho heuristic (residual balancing)
    mu <- 10
    tau_incr <- 2
    tau_decr <- 2
    if (primal_res > mu * dual_res) {
      rho_old <- rho
      rho <- rho * tau_incr
      # scale duals
      U <- U / tau_incr
      # update KKT_LHS_top and KKT_LHS
      KKT_LHS_top <- rho * AtA
      KKT_LHS <- rbind(
        cbind(KKT_LHS_top, t(Aeq)),
        cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
      )
    } else if (dual_res > mu * primal_res) {
      rho_old <- rho
      rho <- rho / tau_decr
      U <- U * tau_decr
      KKT_LHS_top <- rho * AtA
      KKT_LHS <- rbind(
        cbind(KKT_LHS_top, t(Aeq)),
        cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
      )
    }
  }

  # If reached here, max_iter exceeded
  if (verbose) cat("Max iter reached\n")
  return(list(x = x, history = hist, info = list(iter = max_iter, primal = tail(hist$primal,1), dual = tail(hist$dual,1))))
}


#遅すぎる???
admm_sdp02 <- function(P, Q, r, B_list,  rho = 1.0, alpha = 1.0, max_iter = 1000, tol = 1e-4, verbose = TRUE,
                     scaling = c("l2","none","standard","max","ruiz","sinkhorn","custom"),
                     ruiz_iters = 10, sinkhorn_iters = 50, custom_s = NULL) {
  
  # ADMM solver with selectable column-scaling methods (base R only)
  
  scaling <- match.arg(scaling)
  M <- ncol(P)
  K <- nrow(Q)
  N <- length(B_list)
  D <- dim(B_list[[1]])[1]
  if (dim(B_list[[1]])[2] != D) stop("B arrays must be D x D x M")

  eps_small <- 1e-12

  # Build A (stacked D^2 x M blocks) as before
  print("Build A (stacked D^2 x M blocks) as before" )
  nArows <- N * D * D
  A <- matrix(0, nrow = nArows, ncol = M)
  row_idx <- 1
  for (i in seq_len(N)) {
    Bi <- B_list[[i]]
    for (k in seq_len(D)) for (l in seq_len(D)) {
      A[row_idx, ] <- Bi[k, l, ]
      row_idx <- row_idx + 1
    }
  }

  # ---------- Scaling utilities ----------
  # compute column norms or stats on the combined matrix Mbig = rbind(P, Q, A)
  Mbig <- rbind(P, Q, A)

  compute_scaling_vector <- function(method) {
    if (method == "none") return(rep(1, M))
    if (method == "custom") {
      if (is.null(custom_s)) stop("custom_s must be provided for scaling='custom'")
      if (length(custom_s) != M) stop("custom_s length mismatch")
      if (any(custom_s <= 0)) stop("custom_s must be positive")
      return(as.numeric(custom_s))
    }
    if (method == "l2") {
      colnorms <- sqrt(colSums(Mbig^2))
      colnorms[colnorms < eps_small] <- 1.0
      return(colnorms)
    }
    if (method == "standard") {
      # center not meaningful for sparse 0/1 B; use std dev
      colsd <- apply(Mbig, 2, sd)
      colsd[colsd < eps_small] <- 1.0
      return(colsd)
    }
    if (method == "max") {
      colmax <- apply(abs(Mbig), 2, max)
      colmax[colmax < eps_small] <- 1.0
      return(colmax)
    }
    if (method == "ruiz") {
      # Ruiz balancing: find diagonal left (r) and right (c) scalings so that row/col norms ~1
      # We'll return column scaling s = c (positive)
      Atemp <- Mbig
      nr <- nrow(Atemp); nc <- ncol(Atemp)
      rvec <- rep(1, nr); cvec <- rep(1, nc)
      for (it in seq_len(ruiz_iters)) {
        # scale rows to unit norm
        rownorms <- sqrt(rowSums((Atemp * (matrix(cvec, nrow = nr, ncol = nc, byrow = TRUE)))^2))
        rownorms[rownorms < eps_small] <- 1.0
        rscale <- 1 / sqrt(rownorms)
        # apply row scaling
        Atemp <- (matrix(rscale, nrow = nr, ncol = 1) * Atemp)
        rvec <- rvec * rscale
        # scale cols to unit norm
        colnorms <- sqrt(colSums(Atemp^2))
        colnorms[colnorms < eps_small] <- 1.0
        cscale <- 1 / sqrt(colnorms)
        Atemp <- Atemp %*% diag(cscale)
        cvec <- cvec * cscale
      }
      # final column scaling factors: we will scale columns by dividing by cvec
      # but to keep positive and not too small, clip
      cvec[cvec < eps_small] <- 1.0
      return(cvec)
    }
    if (method == "sinkhorn") {
      # Sinkhorn-like row/col sum normalization (requires nonnegativity ideally)
      Atemp <- abs(Mbig) + eps_small
      nr <- nrow(Atemp); nc <- ncol(Atemp)
      rvec <- rep(1, nr); cvec <- rep(1, nc)
      for (it in seq_len(sinkhorn_iters)) {
        # normalize rows to sum 1
        rowsum <- rowSums(Atemp)
        rowsum[rowsum < eps_small] <- 1.0
        Atemp <- Atemp / rowsum
        # normalize cols to sum 1
        colsum <- colSums(Atemp)
        colsum[colsum < eps_small] <- 1.0
        Atemp <- t(t(Atemp) / colsum)
      }
      # derive column scaling as original column sums before normalization
      colsum_orig <- colSums(abs(Mbig)) + eps_small
      colsum_orig[colsum_orig < eps_small] <- 1.0
      return(colsum_orig)
    }
    stop("unknown scaling method")
  }

  # compute scale vector s (positive)
  print("compute scale vector s (positive)" )
  s_vec <- compute_scaling_vector(scaling)
  # ensure positivity and clip
  s_vec[is.na(s_vec) | s_vec <= 0] <- 1.0
  s_vec <- pmax(s_vec, eps_small)

  # We will transform variable: x = S^{-1} z. So replace matrices by M %*% diag(1/s)
  invS_diag <- 1.0 / s_vec
  # scale columns of P, Q, and each B_list element
  P_s <- P %*% diag(invS_diag)
  Q_s <- Q %*% diag(invS_diag)
  # scale B_list in-place: divide each length-M vector by s_j
  B_list_s <- vector("list", N)
  for (i in seq_len(N)) {
    Bi <- B_list[[i]]
    Bi_s <- Bi
    for (k in seq_len(D)) for (l in seq_len(D)) {
      Bi_s[k, l, ] <- Bi[k, l, ] * invS_diag
    }
    B_list_s[[i]] <- Bi_s
  }
  # Rebuild A from scaled B_list_s
  A_s <- matrix(0, nrow = nArows, ncol = M)
  row_idx <- 1
  for (i in seq_len(N)) {
    Bi <- B_list_s[[i]]
    for (k in seq_len(D)) for (l in seq_len(D)) {
      A_s[row_idx, ] <- Bi[k, l, ]
      row_idx <- row_idx + 1
    }
  }

  # Now proceed with ADMM on scaled matrices P_s, Q_s, A_s.
  # For clarity, we copy the earlier ADMM implementation but using P_s, Q_s, A_s.

  # Helper: project initial x to satisfy P_s x = 0, Q_s x = r
  project_affine <- function(x0) {
    Aeq <- rbind(P_s, Q_s)
    meq <- nrow(Aeq)
    KKT <- rbind(
      cbind(diag(M), t(Aeq)),
      cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
    )
    rhs <- c(x0, rep(0, nrow(P_s)), r)
    sol <- tryCatch(solve(KKT, rhs), error = function(e) NULL)
    if (is.null(sol)) stop("KKT solve failed in projection; check P,Q ranks")
    return(sol[1:M])
  }

  project_affine_fast <- function(x0, P, Q, r, tol = 1e-12) {
    Aeq <- rbind(P, Q)
    b <- c(rep(0, nrow(P)), r)

    # 1) compute a particular solution x_p to Aeq x = b
    # Use qr.solve which uses QR factorization (efficient for tall matrices)
    # If system is underdetermined, qr.solve returns a least-squares solution.
    print("1) compute a particular solution x_p to Aeq x = b")
    x_p <- tryCatch(qr.solve(Aeq, b), error = function(e) NULL)
    if (is.null(x_p)) {
      # fallback: use lm.fit (more robust in some cases)
      fit <- lm.fit(Aeq, b)
      x_p <- fit$coefficients
      x_p[is.na(x_p)] <- 0
    }

    # 2) compute nullspace basis Z of Aeq using QR of t(Aeq)
    # qr(t(Aeq)) returns Q (M x M) where columns span R^M; rank = qr$rank
    print("2) compute nullspace basis Z of Aeq using QR of t(Aeq)")
    qr_t <- qr(t(Aeq))
    rank <- qr_t$rank
    M <- ncol(Aeq)
    if (rank < M) {
      Qmat <- qr.Q(qr_t)            # M x M orthonormal
      Z <- Qmat[, (rank+1):M, drop = FALSE]  # M x (M-rank)
      # 3) projection onto affine set
      diff <- x0 - x_p
      x_proj <- x_p + Z %*% (t(Z) %*% diff)
      return(as.numeric(x_proj))
    } else {
      # full rank: unique solution x_p
      return(as.numeric(x_p))
    }
  }

  # 高速版 project_affine_fast2
  project_affine_fast2 <- function(x0, P, Q, r, precomp = NULL,
                                   tol_eig = 1e-10, reg = 1e-12) {
    # Aeq and RHS
    Aeq <- rbind(P, Q)
    b <- c(rep(0, nrow(P)), r)
    M <- ncol(Aeq)

    # If precomp not provided, compute Gram and particular solution and eigen-decomp
    if (is.null(precomp)) {
      # particular solution x_p: least squares solve Aeq x = b
      print("particular solution x_p: least squares solve Aeq x = b")
      x_p <- tryCatch({
        # qr.solve is efficient for tall Aeq
        qr.solve(Aeq, b)
      }, error = function(e) {
        # fallback to lm.fit
        fit <- lm.fit(Aeq, b)
        coef <- fit$coefficients
        coef[is.na(coef)] <- 0
        coef
      })

      # Gram matrix G = Aeq^T Aeq (M x M)
      G <- crossprod(Aeq)  # cheap when M small

      # regularize G slightly for numerical stability
      maxdiag <- max(abs(diag(G)))
      eps_reg <- reg * max(1.0, maxdiag)
      G_reg <- G + eps_reg * diag(M)

      # eigen-decomposition of small MxM matrix
      ev <- eigen(G_reg, symmetric = TRUE)
      vals <- ev$values
      vecs <- ev$vectors

      # determine numerical rank by threshold
      thresh <- max(vals) * tol_eig
      rank <- sum(vals > thresh)

      precomp <- list(x_p = x_p, G = G, G_reg = G_reg,
                      eigvals = vals, eigvecs = vecs, rank = rank,
                      tol_eig = tol_eig, reg = eps_reg)
    } else {
      # use provided precomp; ensure fields exist
      if (is.null(precomp$eigvals) || is.null(precomp$eigvecs) || is.null(precomp$x_p)) {
        stop("precomp must contain eigvals, eigvecs, and x_p")
      }
      x_p <- precomp$x_p
    }

    # If full column rank, unique solution x_p is the projection
    if (precomp$rank >= M) {
      return(list(x_proj = as.numeric(precomp$x_p), precomp = precomp))
    }

    # Build nullspace basis Z from eigenvectors corresponding to near-zero eigenvalues
    # eigenvectors are columns of eigvecs; select indices where eigvals <= thresh
    vals <- precomp$eigvals
    vecs <- precomp$eigvecs
    thresh <- max(vals) * precomp$tol_eig
    null_idx <- which(vals <= thresh)
    if (length(null_idx) == 0) {
      # numerical corner: treat as full rank
      return(list(x_proj = as.numeric(precomp$x_p), precomp = precomp))
    }
    Z <- vecs[, null_idx, drop = FALSE]  # M x (M-rank)

    # projection: x = x_p + Z Z^T (x0 - x_p)
    diff <- x0 - precomp$x_p
    # compute Z^T diff (small)
    ztd <- crossprod(Z, diff)   # (M-rank) vector
    x_proj <- precomp$x_p + Z %*% ztd

    return(list(x_proj = as.numeric(x_proj), precomp = precomp))
  }

  # initialize
  print("initialize")
  x_scaled <- project_affine_fast(rep(0, M),P,Q,r)
  Y_vec <- rep(0, nArows)
  U <- rep(0, nArows)
  At <- t(A_s)
  AtA <- At %*% A_s
  Aeq <- rbind(P_s, Q_s)
  meq <- nrow(Aeq)
  KKT_LHS_top <- rho * AtA
  KKT_LHS <- rbind(
    cbind(KKT_LHS_top, t(Aeq)),
    cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
  )
  cvec <- rep(0, M); cvec[1] <- 1.0  # minimize x1 -> in scaled variable minimize z1 (since x = S^{-1} z)
  hist <- list(primal = numeric(), dual = numeric())

  print("iteration start")
  for (k in seq_len(max_iter)) {
    # x (scaled variable z) update
    b_vec <- Y_vec - U
    Atb <- At %*% b_vec
    rhs_top <- -cvec + rho * Atb
    rhs <- c(rhs_top, rep(0, nrow(P_s)), r)
    sol <- tryCatch(solve(KKT_LHS, rhs), error = function(e) NULL)
    if (is.null(sol)) stop("KKT solve failed in x-update; consider regularization")
    z_old <- x_scaled
    x_scaled <- sol[1:M]

    # Ax and over-relaxation
    Ax <- A_s %*% x_scaled
    Ax_hat <- alpha * Ax + (1 - alpha) * (Y_vec)

    # Y-update with PSD projection per block
    Y_prev <- Y_vec
    primal_res_norm_sq <- 0
    for (i in seq_len(N)) {
      start <- (i-1)*D*D + 1
      end <- start + D*D - 1
      s_vec_block <- Ax_hat[start:end] + U[start:end]
      S <- matrix(s_vec_block, nrow = D, ncol = D, byrow = TRUE)
      S_sym <- (S + t(S)) / 2
      ev <- eigen(S_sym, symmetric = TRUE)
      vals_proj <- pmax(ev$values, 0)
      Ymat <- ev$vectors %*% diag(vals_proj) %*% t(ev$vectors)
      Yv <- as.vector(t(Ymat))
      Y_vec[start:end] <- Yv
      primal_res_norm_sq <- primal_res_norm_sq + sum((Ax[start:end] - Yv)^2)
    }
    primal_res <- sqrt(primal_res_norm_sq)

    # Dual update
    U <- U + (Ax_hat - Y_vec)
    dual_res <- sqrt(sum((rho * (At %*% (Y_vec - Y_prev)))^2))

    hist$primal <- c(hist$primal, primal_res)
    hist$dual <- c(hist$dual, dual_res)
    if (verbose && (k %% 10 == 0 || k == 1)) {
      cat(sprintf("iter %4d: primal = %.3e, dual = %.3e, rho = %.3e\n", k, primal_res, dual_res, rho))
    }

    # stopping
    if (primal_res < tol && dual_res < tol) {
      if (verbose) cat("Converged at iter", k, "\n")
      break
    }

    # adaptive rho (residual balancing)
    mu <- 10; tau_incr <- 2; tau_decr <- 2
    if (primal_res > mu * dual_res) {
      rho <- rho * tau_incr
      U <- U / tau_incr
      KKT_LHS_top <- rho * AtA
      KKT_LHS <- rbind(
        cbind(KKT_LHS_top, t(Aeq)),
        cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
      )
    } else if (dual_res > mu * primal_res) {
      rho <- rho / tau_decr
      U <- U * tau_decr
      KKT_LHS_top <- rho * AtA
      KKT_LHS <- rbind(
        cbind(KKT_LHS_top, t(Aeq)),
        cbind(Aeq, matrix(0, nrow = meq, ncol = meq))
      )
    }
  }

  # Recover original-scale x from scaled solution: x = S^{-1} z
  x_sol <- x_scaled * invS_diag  # since x = S^{-1} z and invS_diag = 1/s
  info <- list(iter = k, primal = tail(hist$primal, 1), dual = tail(hist$dual, 1), scaling = scaling)
  return(list(x = x_sol, history = hist, info = info, s = s_vec))
}



# 
# 
# 
# 
# # soft-threshold
# soft_threshold <- function(v, tau) {
#   sign(v) * pmax(abs(v) - tau, 0)
# }
# 
# # ADMM solver without Matrix package
# admm_sparse_residual_baseR <- function(P, Q, b,
#                                        lambda = 1.0, rho = 1.0,
#                                        eps_reg = 1e-8,
#                                        max_iter = 2000,
#                                        eps_abs = 1e-6, eps_rel = 1e-3,
#                                        x0 = NULL, y0 = NULL,
#                                        auto_rho = TRUE,
#                                        rho_update_factor = 2,
#                                        rho_balance_tol = 10,
#                                        verbose = FALSE) {
#   n_rows <- nrow(P)
#   n_col  <- ncol(P)
#   
#   # Column scaling (norm) and store for rescaling x at the end
#   colnorms <- sqrt(colSums(P^2))
#   colnorms[colnorms == 0] <- 1
#   Psc <- sweep(P, 2, colnorms, "/")    # scaled P
#   Pt_sc <- t(Psc)
#   G <- crossprod(Psc)                  # n_col x n_col
#   
#   # Projection helpers for Q x = b (Q is 1 x n_col)
#   Q_vec <- as.numeric(Q)
#   QQt <- as.numeric(Q_vec %*% Q_vec)   # scalar
#   Qt <- as.numeric(Q_vec)              # length n_col
#   
#   # Precompute matrix for x-update: (2 G + rho I + eps_reg I)
#   A_mat <- 2 * G + rho * diag(n_col) + eps_reg * diag(n_col)
#   cholA <- chol(A_mat)
#   
#   # Initialize x, y, u (scaled dual)
#   if (!is.null(x0)) {
#     if (length(x0) != n_col) stop("x0 length mismatch")
#     x <- as.numeric(x0)
#     # project to satisfy Qx=b
#     qx_minus_b <- as.numeric(Q_vec %*% x - b)
#     if (abs(qx_minus_b) > 0) x <- x - (qx_minus_b / QQt) * Qt
#     x_sc <- x * colnorms
#   } else {
#     x_sc <- rep(0, n_col)
#   }
#   
#   if (!is.null(y0)) {
#     if (length(y0) != n_rows) stop("y0 length mismatch")
#     y_dense <- as.numeric(y0)
#   } else {
#     y_dense <- rep(0, n_rows)
#   }
#   
#   # scaled dual u = Px - y (dense)
#   Px_sc <- as.numeric(Psc %*% x_sc)
#   u <- Px_sc - y_dense
#   
#   sqrt_n <- sqrt(n_rows)
#   sqrt_m <- sqrt(n_col)
#   
#   for (k in 1:max_iter) {
#     # x-update: solve (2 G + rho I) x_sc = 2 P' y + rho P' (y + u)
#     rhs <- 2 * (Pt_sc %*% y_dense) + rho * (Pt_sc %*% (y_dense + u))
#     x_sc_new <- backsolve(cholA, forwardsolve(t(cholA), rhs))
#     
#     # project to satisfy Q x = b (x in original scale)
#     x_candidate <- x_sc_new / colnorms
#     qx_minus_b <- as.numeric(Q_vec %*% x_candidate - b)
#     if (abs(qx_minus_b) > 0) {
#       x_candidate <- x_candidate - (qx_minus_b / QQt) * Qt
#       x_sc_new <- x_candidate * colnorms
#     }
#     
#     # y-update: soft-threshold on P x - u
#     Px_sc_new <- as.numeric(Psc %*% x_sc_new)
#     v <- Px_sc_new - u
#     tau <- lambda / rho
#     y_new_dense <- soft_threshold(v, tau)
#     
#     # dual update u <- u + (y_new - Px_new)
#     u_new <- u + (y_new_dense - Px_sc_new)
#     
#     # residuals
#     prim_res <- sqrt(sum((Px_sc_new - y_new_dense)^2))
#     dvec <- Pt_sc %*% (y_new_dense - y_dense)
#     dual_res <- sqrt(sum((rho * dvec)^2))
#     
#     # stopping tolerances (absolute + relative)
#     eps_pri <- sqrt_n * eps_abs + eps_rel * max(sqrt(sum(Px_sc_new^2)), sqrt(sum(y_new_dense^2)))
#     eps_dual <- sqrt_m * eps_abs + eps_rel * sqrt(sum((rho * Pt_sc %*% u_new)^2))
#     
#     if (verbose && (k %% 50 == 0 || k == 1)) {
#       nz_count <- sum(abs(y_new_dense) > 0)
#       cat(sprintf("iter %4d: prim_res=%.3e dual_res=%.3e eps_pri=%.3e eps_dual=%.3e nz=%d rho=%.3e\n",
#                   k, prim_res, dual_res, eps_pri, eps_dual, nz_count, rho))
#     }
#     
#     # numerical checks
#     if (!is.finite(prim_res) || !is.finite(dual_res) || any(!is.finite(x_sc_new)) || any(!is.finite(u_new))) {
#       warning("Numerical instability detected; consider increasing eps_reg, scaling P, or reducing rho")
#       break
#     }
#     
#     # convergence check
#     if (prim_res <= eps_pri && dual_res <= eps_dual) {
#       x_sc <- x_sc_new; y_dense <- y_new_dense; u <- u_new
#       break
#     }
#     
#     # auto rho update (Boyd rule)
#     if (auto_rho) {
#       if (prim_res > rho_balance_tol * dual_res) {
#         rho <- rho * rho_update_factor
#         u <- u / rho_update_factor
#         A_mat <- 2 * G + rho * diag(n_col) + eps_reg * diag(n_col)
#         cholA <- chol(A_mat)
#       } else if (dual_res > rho_balance_tol * prim_res) {
#         rho <- rho / rho_update_factor
#         u <- u * rho_update_factor
#         A_mat <- 2 * G + rho * diag(n_col) + eps_reg * diag(n_col)
#         cholA <- chol(A_mat)
#       }
#     }
#     
#     # update variables
#     x_sc <- x_sc_new
#     y_dense <- y_new_dense
#     u <- u_new
#   }
#   
#   # final x in original scale
#   x_final <- x_sc / colnorms
#   y_final <- y_dense
#   support <- which(abs(y_final) > 0)
#   
#   # support refinement: re-solve x with final y fixed
#   if (length(support) > 0) {
#     rhs2 <- 2 * (Pt_sc %*% y_final)
#     x_sc_ref <- backsolve(cholA, forwardsolve(t(cholA), rhs2))
#     x_ref <- x_sc_ref / colnorms
#     qx_minus_b <- as.numeric(Q_vec %*% x_ref - b)
#     if (abs(qx_minus_b) > 0) x_ref <- x_ref - (qx_minus_b / QQt) * Qt
#     x_final <- x_ref
#   }
#   
#   list(x = as.numeric(x_final),
#        y = as.numeric(y_final),
#        support = support,
#        iter = k,
#        prim_res = prim_res,
#        dual_res = dual_res,
#        rho = rho)
# }


