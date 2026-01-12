
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
    
    sol <- solve(KKT, rhs)
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
      solve(KKT, rhs_kkt)
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



# soft-threshold
soft_threshold <- function(v, tau) {
  sign(v) * pmax(abs(v) - tau, 0)
}

# ADMM solver without Matrix package
admm_sparse_residual_baseR <- function(P, Q, b,
                                       lambda = 1.0, rho = 1.0,
                                       eps_reg = 1e-8,
                                       max_iter = 2000,
                                       eps_abs = 1e-6, eps_rel = 1e-3,
                                       x0 = NULL, y0 = NULL,
                                       auto_rho = TRUE,
                                       rho_update_factor = 2,
                                       rho_balance_tol = 10,
                                       verbose = FALSE) {
  n_rows <- nrow(P)
  n_col  <- ncol(P)
  
  # Column scaling (norm) and store for rescaling x at the end
  colnorms <- sqrt(colSums(P^2))
  colnorms[colnorms == 0] <- 1
  Psc <- sweep(P, 2, colnorms, "/")    # scaled P
  Pt_sc <- t(Psc)
  G <- crossprod(Psc)                  # n_col x n_col
  
  # Projection helpers for Q x = b (Q is 1 x n_col)
  Q_vec <- as.numeric(Q)
  QQt <- as.numeric(Q_vec %*% Q_vec)   # scalar
  Qt <- as.numeric(Q_vec)              # length n_col
  
  # Precompute matrix for x-update: (2 G + rho I + eps_reg I)
  A_mat <- 2 * G + rho * diag(n_col) + eps_reg * diag(n_col)
  cholA <- chol(A_mat)
  
  # Initialize x, y, u (scaled dual)
  if (!is.null(x0)) {
    if (length(x0) != n_col) stop("x0 length mismatch")
    x <- as.numeric(x0)
    # project to satisfy Qx=b
    qx_minus_b <- as.numeric(Q_vec %*% x - b)
    if (abs(qx_minus_b) > 0) x <- x - (qx_minus_b / QQt) * Qt
    x_sc <- x * colnorms
  } else {
    x_sc <- rep(0, n_col)
  }
  
  if (!is.null(y0)) {
    if (length(y0) != n_rows) stop("y0 length mismatch")
    y_dense <- as.numeric(y0)
  } else {
    y_dense <- rep(0, n_rows)
  }
  
  # scaled dual u = Px - y (dense)
  Px_sc <- as.numeric(Psc %*% x_sc)
  u <- Px_sc - y_dense
  
  sqrt_n <- sqrt(n_rows)
  sqrt_m <- sqrt(n_col)
  
  for (k in 1:max_iter) {
    # x-update: solve (2 G + rho I) x_sc = 2 P' y + rho P' (y + u)
    rhs <- 2 * (Pt_sc %*% y_dense) + rho * (Pt_sc %*% (y_dense + u))
    x_sc_new <- backsolve(cholA, forwardsolve(t(cholA), rhs))
    
    # project to satisfy Q x = b (x in original scale)
    x_candidate <- x_sc_new / colnorms
    qx_minus_b <- as.numeric(Q_vec %*% x_candidate - b)
    if (abs(qx_minus_b) > 0) {
      x_candidate <- x_candidate - (qx_minus_b / QQt) * Qt
      x_sc_new <- x_candidate * colnorms
    }
    
    # y-update: soft-threshold on P x - u
    Px_sc_new <- as.numeric(Psc %*% x_sc_new)
    v <- Px_sc_new - u
    tau <- lambda / rho
    y_new_dense <- soft_threshold(v, tau)
    
    # dual update u <- u + (y_new - Px_new)
    u_new <- u + (y_new_dense - Px_sc_new)
    
    # residuals
    prim_res <- sqrt(sum((Px_sc_new - y_new_dense)^2))
    dvec <- Pt_sc %*% (y_new_dense - y_dense)
    dual_res <- sqrt(sum((rho * dvec)^2))
    
    # stopping tolerances (absolute + relative)
    eps_pri <- sqrt_n * eps_abs + eps_rel * max(sqrt(sum(Px_sc_new^2)), sqrt(sum(y_new_dense^2)))
    eps_dual <- sqrt_m * eps_abs + eps_rel * sqrt(sum((rho * Pt_sc %*% u_new)^2))
    
    if (verbose && (k %% 50 == 0 || k == 1)) {
      nz_count <- sum(abs(y_new_dense) > 0)
      cat(sprintf("iter %4d: prim_res=%.3e dual_res=%.3e eps_pri=%.3e eps_dual=%.3e nz=%d rho=%.3e\n",
                  k, prim_res, dual_res, eps_pri, eps_dual, nz_count, rho))
    }
    
    # numerical checks
    if (!is.finite(prim_res) || !is.finite(dual_res) || any(!is.finite(x_sc_new)) || any(!is.finite(u_new))) {
      warning("Numerical instability detected; consider increasing eps_reg, scaling P, or reducing rho")
      break
    }
    
    # convergence check
    if (prim_res <= eps_pri && dual_res <= eps_dual) {
      x_sc <- x_sc_new; y_dense <- y_new_dense; u <- u_new
      break
    }
    
    # auto rho update (Boyd rule)
    if (auto_rho) {
      if (prim_res > rho_balance_tol * dual_res) {
        rho <- rho * rho_update_factor
        u <- u / rho_update_factor
        A_mat <- 2 * G + rho * diag(n_col) + eps_reg * diag(n_col)
        cholA <- chol(A_mat)
      } else if (dual_res > rho_balance_tol * prim_res) {
        rho <- rho / rho_update_factor
        u <- u * rho_update_factor
        A_mat <- 2 * G + rho * diag(n_col) + eps_reg * diag(n_col)
        cholA <- chol(A_mat)
      }
    }
    
    # update variables
    x_sc <- x_sc_new
    y_dense <- y_new_dense
    u <- u_new
  }
  
  # final x in original scale
  x_final <- x_sc / colnorms
  y_final <- y_dense
  support <- which(abs(y_final) > 0)
  
  # support refinement: re-solve x with final y fixed
  if (length(support) > 0) {
    rhs2 <- 2 * (Pt_sc %*% y_final)
    x_sc_ref <- backsolve(cholA, forwardsolve(t(cholA), rhs2))
    x_ref <- x_sc_ref / colnorms
    qx_minus_b <- as.numeric(Q_vec %*% x_ref - b)
    if (abs(qx_minus_b) > 0) x_ref <- x_ref - (qx_minus_b / QQt) * Qt
    x_final <- x_ref
  }
  
  list(x = as.numeric(x_final),
       y = as.numeric(y_final),
       support = support,
       iter = k,
       prim_res = prim_res,
       dual_res = dual_res,
       rho = rho)
}


