source("../../lg.R")
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
# 
# library(gridExtra)
# library(CVXR)

source("./basis2.R")
source("./nss.R")
source("./admm.R")

# データ準備
# df = read.csv("model2_sol.csv",header = TRUE,stringsAsFactors = FALSE)
# df = read.csv("modelX_with_forces_1_0_-1_1_sol.csv",header = TRUE,stringsAsFactors = FALSE)
# df = read.csv("model2_g_1_0_-1_1_sol.csv",header = TRUE,stringsAsFactors = FALSE)
df = read.csv("model2_gtr_1_0_-1_1_sol.csv",header = TRUE,stringsAsFactors = FALSE)

df %>% head()
# dn();df %>% ggplot(aes(x=x1,y=x2))+geom_point()

t0 = df$t
dt = t0[2]-t0[1]

# x0 = df$x1
# y0 = df$x2
# x0 = df$x1/14
# y0 = df$x2/2.5
# x0 = 2*(df$x1-min(df$x1))/(max(df$x1)-min(df$x1)) - 1
# y0 = 2*(df$x2-min(df$x2))/(max(df$x2)-min(df$x2)) - 1
# x0 = 1.5*(df$x1-min(df$x1))/(max(df$x1)-min(df$x1)) -0.75
# y0 = 1.5*(df$x2-min(df$x2))/(max(df$x2)-min(df$x2)) -0.75
# x0 = (df$x1-min(df$x1))/(max(df$x1)-min(df$x1))
# y0 = (df$x2-min(df$x2))/(max(df$x2)-min(df$x2))
# x0 = 0.6*(df$x1-min(df$x1))/(max(df$x1)-min(df$x1)) + 0.2
# y0 = 0.6*(df$x2-min(df$x2))/(max(df$x2)-min(df$x2)) + 0.2

# const1 = 0.6*(0-min(df$x1))/(max(df$x1)-min(df$x1)) + 0.2
# const2 = 0.6*(0-min(df$x2))/(max(df$x2)-min(df$x2)) + 0.2
# gain1 = 0.6/(max(df$x1)-min(df$x1)) 
# gain2 = 0.6/(max(df$x2)-min(df$x2)) 
# const1 = 1.6*(0-min(df$x1))/(max(df$x1)-min(df$x1)) - 0.8
# const2 = 1.6*(0-min(df$x2))/(max(df$x2)-min(df$x2)) - 0.8
# gain1 = 1.6/(max(df$x1)-min(df$x1)) 
# gain2 = 1.6/(max(df$x1)-min(df$x1)) 
# const1 = 0
# const2 = 0
# gain1 = 0.8/max(abs(df$x1))
# gain2 = 0.8/max(abs(df$x2))
const1 = 0
const2 = 0
gain1 = 1/30
gain2 = 1/4

x0 = gain1 * df$x1 + const1
y0 = gain2 * df$x2 + const2

H1 = df$H[1]
tr1 = df$g11[1]/(gain1^2) + df$g22[1]/(gain2^2)


N = length(x0)

dn();matplot(t0,cbind(x0,y0))
dn();plot(x0,y0)

#時間微分の計算

xd = (x0[3:N]-x0[1:(N-2)])/(2*dt)
xdd = (x0[3:N]-2*x0[2:(N-1)]+x0[1:(N-2)])/(dt*dt)
yd = (y0[3:N]-y0[1:(N-2)])/(2*dt)
ydd = (y0[3:N]-2*y0[2:(N-1)]+y0[1:(N-2)])/(dt*dt)
x = x0[2:(N-1)]
y = y0[2:(N-1)]
t = t0[2:(N-1)]

X = matrix(c(x0[2:(N-1)],y0[2:(N-1)]),ncol=2)
Xd = matrix(c(xd,yd),ncol=2)
Xdd = matrix(c(xdd,ydd),ncol=2)

Tx = dim(X)[1]

# 基底関数の準備

phi = phi0
phid = phi0d
phi_name = phi0_name
# phi = phiA
# phid = phiAd
# phi_name = phiA_name
# phi = phiC
# phid = phiCd
# phi_name = phiC_name
# phi = bnst[[4]]
# phid = bnstd[[4]]
# phi_name = bnst_name[[4]]
# phi = nlgdr
# phid = nlgdrd
# phi_name = nlgdr_name

psi = phiB
psid = phiBd
psi_name = phiB_name
# psi = phiD
# psid = phiDd
# psi_name = phiD_name

Mm = 2 #6 10 15 21 28
Mp = 0 # 2 5 9 14
n = 2

na = n*(n+1)*Mm/2
nb = Mp


# index作成

idx <- array(0, dim = c(n, n, Mm))
jdx <- array(0, dim = Mp)
tdx <- array(0, dim = c(Tx,n))

kk=1
for(ii in 1:n){
  for(jj in ii:n){
    for(ll in 1:Mm){
      idx[ii,jj,ll] = kk
      idx[jj,ii,ll] = kk
      kk=kk+1
    }
  }
}
if(Mp>0){
  for(ll in 1:Mp){
    jdx[ll] = kk
    kk=kk+1
  }
}
vv=1
for(kk in 1:n){
  for(tt in 1:Tx){
    tdx[tt,kk] = vv
    vv=vv+1
  }
}

# 基底関数の名前
coef_name = rep("",na+nb)
for(ii in 1:n){
  for(jj in ii:n){
    for(ll in 1:Mm){
      coef_name[idx[ii,jj,ll]] = paste0("a",ii,jj,"_",phi_name[ll])
    }
  }
}
if(Mp>0){
  for(ll in 1:Mp){
    coef_name[ jdx[ll] ] = paste0("b","_",psi_name[ll])
  }
}


#二次元専用
B_list = as.list(1:Mm)
B2_list = as.list(1:Mm)
for(kk in 1:Mm){
  B_list[[kk]] = matrix(0,ncol=na+nb,nrow=3)
  B_list[[kk]][1,idx[1,1,kk]] = 1
  B_list[[kk]][2,idx[1,2,kk]] = 1
  B_list[[kk]][3,idx[2,2,kk]] = 1
  
  B2_list[[kk]] = array(0, dim = c(n, n, na+nb))
  B2_list[[kk]][1,1,idx[1,1,kk]] = 1
  B2_list[[kk]][1,2,idx[1,2,kk]] = 1
  B2_list[[kk]][2,1,idx[2,1,kk]] = 1
  B2_list[[kk]][2,2,idx[2,2,kk]] = 1
  
}


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
if(Mp>0){
  gpsi = array(0,dim=c(Mp,Tx))
  gpsid = array(0,dim=c(n,Mp,Tx))
  for(ll in 1:Mp){
    for(tt in 1:Tx){
      gpsi[ll,tt] = psi[[ll]](X[tt,1],X[tt,2])
      for(kk in 1:n){
        gpsid[kk,ll,tt] = psid[[kk]][[ll]](X[tt,1],X[tt,2])
      }
    }
  }
}


# 係数行列作成

# P = matrix(0,nrow=Tx*n, ncol=na+nb)
Pa = matrix(0,nrow=Tx*n, ncol=na+nb)
Pv = matrix(0,nrow=Tx*n, ncol=na+nb)
Pf = matrix(0,nrow=Tx*n, ncol=na+nb)
v = matrix(0,nrow=Tx*n, ncol=1)

# a
for(tt in 1:Tx){
  for(kk in 1:n){
    for(ii in 1:n){
      for(ll in 1:Mm){
        uu = idx[kk,ii,ll]
        vv = tdx[tt,kk]
        S = Xdd[tt,ii]*gphi[ll,tt]
        Pa[vv,uu] = Pa[vv,uu] + S
        
        for(jj in 1:n){
          uu = idx[ii,kk,ll]
          S = 0.5*Xd[tt,ii]*Xd[tt,jj]*gphid[jj,ll,tt]
          Pv[vv,uu] = Pv[vv,uu] + S
          
          uu = idx[kk,jj,ll]
          S = 0.5*Xd[tt,ii]*Xd[tt,jj]*gphid[ii,ll,tt]
          Pv[vv,uu] = Pv[vv,uu] + S
          
          uu = idx[ii,jj,ll]
          S = -0.5*Xd[tt,ii]*Xd[tt,jj]*gphid[kk,ll,tt]
          Pv[vv,uu] = Pv[vv,uu] + S
          
        }
      }
    }
  }
}

# b
if(Mp>0){
  for(tt in 1:Tx){
    for(kk in 1:n){
      for(ll in 1:Mp){
        vv = tdx[tt,kk]
        uu = jdx[ll]
        S = gpsid[kk,ll,tt]
        Pf[vv,uu] = Pf[vv,uu] + S
      }
    }
  }
}

P = Pa + Pv + Pf

# trace = fix
PT1 = matrix(0,nrow=1, ncol=na+nb)
KT = matrix(0,nrow=Tx, ncol=na+nb)
vT1 = matrix(tr1 ,nrow=1, ncol=1)
for(tt in 1:Tx){
  if(tt==1){
    for(ii in 1:n){
      for(ll in 1:Mm){
        uu = idx[ii,ii,ll]
        S = gphi[ll,tt]
        PT1[1,uu] = PT1[1,uu] + S
      }
    }
  }
  for(ii in 1:n){
    for(ll in 1:Mm){
      uu = idx[ii,ii,ll]
      S = gphi[ll,tt]
      KT[tt,uu] = KT[tt,uu] + S
    }
  }
}


# H = const
PHc = matrix(0,nrow=Tx-1, ncol=na+nb)
vHc = matrix(0,nrow=Tx-1, ncol=1)
for(tt in 1:(Tx-1)){
  for(ii in 1:n){
    for(jj in 1:n){
      for(ll in 1:Mm){
        uu = idx[ii,jj,ll]
        S = gphi[ll,tt]*Xd[tt,ii]*Xd[tt,jj]/2 - gphi[ll,tt+1]*Xd[tt+1,ii]*Xd[tt+1,jj]/2
        PHc[tt,uu] = PHc[tt,uu] + S
      }
    }
  }
  if(Mp>0){
    for(ll in 1:Mp){
      uu = jdx[ll]
      S = gpsi[ll,tt] - gpsi[ll,tt+1]
      PHc[tt,uu] = PHc[tt,uu] + S
    }
  }
}

# H = Fix
PH1 = matrix(0,nrow=1, ncol=na+nb)
KH = matrix(0,nrow=Tx, ncol=na+nb)
vH1 = matrix( H1 , nrow=1, ncol=1)
for(tt in 1:Tx){
  if(tt==1){
    for(ii in 1:n){
      for(jj in 1:n){
        for(ll in 1:Mm){
          uu = idx[ii,jj,ll]
          S = gphi[ll,tt]*Xd[tt,ii]*Xd[tt,jj]/2 
          PH1[tt,uu] = PH1[tt,uu] + S
        }
      }
    }
    if(Mp>0){
      for(ll in 1:Mp){
        uu = jdx[ll]
        S = gpsi[ll,tt] 
        PH1[tt,uu] = PH1[tt,uu] + S
      }
    }
  }
  for(ii in 1:n){
    for(jj in 1:n){
      for(ll in 1:Mm){
        uu = idx[ii,jj,ll]
        S = gphi[ll,tt]*Xd[tt,ii]*Xd[tt,jj]/2 
        KH[tt,uu] = KH[tt,uu] + S
      }
    }
  }
  if(Mp>0){
    for(ll in 1:Mp){
      uu = jdx[ll]
      S = gpsi[ll,tt] 
      KH[tt,uu] = KH[tt,uu] + S
    }
  }
}


P0 = rbind(PHc,PT1)
v0 = rbind(vHc,vT1)
nP0 = dim(P0)[1]

Q = rbind(PHc,PT1)
b = rbind(vHc,vT1)
M = rbind(P,Q)
u = rbind(v,b)

A1 = rbind(
  cbind( t(P)%*%P, t(PH1)),
  cbind( PH1,       0    ))
B1 = rbind(matrix(0,nrow=dim(P)[2],ncol=1),vH1)

A2 = rbind(
  cbind( t(P)%*%P, t(PT1)),
  cbind( PT1,       0    ))
B2 = rbind(matrix(0,nrow=dim(P)[2],ncol=1),vT1)

c_1 = solve(A1,B1,tol=0) #direct H
c_2 = solve(A2,B2,tol=0) #direct T

iPtP = solve(t(P)%*%P,tol=0)
Sc1 = PH1%*%iPtP%*%t(PH1)
Sc2 = PT1%*%iPtP%*%t(PT1)

c_3 = iPtP%*%t(PH1)%*%solve(Sc1,vH1,tol=0) # Schur complement H
c_4 = iPtP%*%t(PT1)%*%solve(Sc2,vT1,tol=0) # Schur complement T
c_5 = nullspace_sol(P,PH1,vH1 ) # null space H
c_6 = nullspace_sol(P,PT1,vT1 ) # null space T

c_7 = solve(t(rbind(P,PH1))%*%rbind(P,PH1),t(rbind(P,PH1))%*%rbind(v,vH1),tol = 0) # penalty H
c_8 = solve(t(rbind(P,PT1))%*%rbind(P,PT1),t(rbind(P,PT1))%*%rbind(v,vT1),tol = 0) # penalty T

c_9 = nullspace_sol(P,Q,b ) #null space ,T C

c_11 = solve(t(M)%*%M,t(M)%*%u,tol = 0) # penalty T C

sv <- svd(M) 
c_12 <- sv$v %*% diag(1/sv$d) %*% t(sv$u) %*% u # svd T C
sv <- svd(rbind(P,PH1)) 
c_12b <- sv$v %*% diag(1/sv$d) %*% t(sv$u) %*% rbind(v,vH1) # svd H
sv <- svd(rbind(P,PT1)) 
c_12c <- sv$v %*% diag(1/sv$d) %*% t(sv$u) %*% rbind(v,vT1) # svd T

res = eigen(t(P)%*%P)
scl = PT1%*%res$vectors[,length(res$values)]
c_13 = res$vectors[,length(res$values)]%*%(vT1/scl) #eigen T
scl = PH1%*%res$vectors[,length(res$values)]
c_13b = res$vectors[,length(res$values)]%*%(vH1/scl) #eigen H

eps13T = sqrt(sum((P %*% c_13)^2))
eps13H = sqrt(sum((P %*% c_13b)^2))
# adm = admm01_fast_m_large(A=P,Q=PT1,r=vT1,eps=sqrt(sum((P %*% c_13)^2)),max_iter = 1000,rho0 = 0.1)
adm = admm02(A=P,Q=PT1,r=vT1,eps=eps13T,max_iter = 1000,rho0 = 0.1)
c_14a = adm$x #admm T
adm = admm02(A=P,Q=PT1,r=vT1,eps=eps13T,max_iter = 1000,rho0 = 1.0)
c_14b = adm$x #admm T
adm = admm02(A=P,Q=PT1,r=vT1,eps=eps13T,max_iter = 1000,rho0 = 10.0)
c_14c = adm$x #admm T

adm2 = admm02(A=P,Q=PH1,r=vH1,eps=eps13H,max_iter = 1000,rho0 = 0.1)
c_14d = adm2$x #admm H
adm2 = admm02(A=P,Q=PH1,r=vH1,eps=eps13H,max_iter = 1000,rho0 = 1.0)
c_14e = adm2$x #admm H
adm2 = admm02(A=P,Q=PH1,r=vH1,eps=eps13H,max_iter = 1000,rho0 = 10.0)
c_14f = adm2$x #admm H

adm3 = admm_with_psd(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 0.1, rho_s = 0.1)
c_15a = adm3$x #admm PSD T
adm3 = admm_with_psd(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 1.0, rho_s = 1.0)
c_15b = adm3$x #admm PSD T
adm3 = admm_with_psd(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 10, rho_s = 10)
c_15c = adm3$x #admm PSD T

adm3 = admm_with_psd(A=P,Q=PH1,r=vH1,B_list=B_list,eps=eps13H, max_iter = 1000, rho = 0.1, rho_s = 0.1)
c_15d = adm3$x #admm PSD H
adm3 = admm_with_psd(A=P,Q=PH1,r=vH1,B_list=B_list,eps=eps13H, max_iter = 1000, rho = 1.0, rho_s = 1.0)
c_15e = adm3$x #admm PSD H
adm3 = admm_with_psd(A=P,Q=PH1,r=vH1,B_list=B_list,eps=eps13H, max_iter = 1000, rho = 10, rho_s = 10)
c_15f = adm3$x #admm PSD H

adm3 = admm_with_psd_mode(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 0.1, rho_s = 0.1, scale_mode = "transform")
c_16a = adm3$x #admm PSD T
adm3 = admm_with_psd_mode(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 1.0, rho_s = 1.0, scale_mode = "transform")
c_16b = adm3$x #admm PSD T
adm3 = admm_with_psd_mode(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 10, rho_s = 10, scale_mode = "transform")
c_16c = adm3$x #admm PSD T

adm3 = admm_with_psd_mode(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 0.1, rho_s = 0.1, scale_mode = "normalize")
c_16d = adm3$x #admm PSD T
adm3 = admm_with_psd_mode(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 1.0, rho_s = 1.0, scale_mode = "normalize")
c_16e = adm3$x #admm PSD T
adm3 = admm_with_psd_mode(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 10, rho_s = 10, scale_mode = "normalize")
c_16f = adm3$x #admm PSD T


adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 0.1, rho_s = 0.1, scale_mode = "transform")
c_17a = adm4$x #admm PSD T
adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 1.0, rho_s = 1.0, scale_mode = "transform")
c_17b = adm4$x #admm PSD T
adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 10, rho_s = 10, scale_mode = "transform")
c_17c = adm4$x #admm PSD T

adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 0.1, rho_s = 0.1, scale_mode = "normalize")
c_17d = adm4$x #admm PSD T
adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 1.0, rho_s = 1.0, scale_mode = "normalize")
c_17e = adm4$x #admm PSD T
adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 10, rho_s = 10, scale_mode = "normalize")
c_17f = adm4$x #admm PSD T

adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 0.1, rho_s = 0.1, scale_mode = "none")
c_17g = adm4$x #admm PSD T
adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 1.0, rho_s = 1.0, scale_mode = "none")
c_17h = adm4$x #admm PSD T
adm4 = admm_with_psd_mode_l2(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, max_iter = 1000, rho = 10, rho_s = 10, scale_mode = "none")
c_17i = adm4$x #admm PSD T

eta13T = sqrt(sum((PT1 %*% c_13 - vT1)^2))
eta13H = sqrt(sum((PH1 %*% c_13 - vH1)^2))

# c_18a = solve_admm_3dB(P=P, Q=PT1, r=vT1, B_list=B2_list, eps=eps13T, eta=eta13T, scale_mode = "none",max_iter = 2000)$x
# c_18b = solve_admm_3dB(P=P, Q=PT1, r=vT1, B_list=B2_list, eps=eps13T, eta=eta13T, scale_mode = "l2",max_iter = 2000)$x
# c_18c = solve_admm_3dB(P=P, Q=PT1, r=vT1, B_list=B2_list, eps=eps13T, eta=eta13T, scale_mode = "maxabs",max_iter = 2000)$x
# c_18d = solve_admm_3dB(P=P, Q=PT1, r=vT1, B_list=B2_list, eps=eps13T, eta=eta13T, scale_mode = "whiten",max_iter = 2000)$x
c_18a = admm_with_psd_mode_02(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, scale_mode = "normalize", adapt_method = "residual")$x
c_18b = admm_with_psd_mode_02(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, scale_mode = "transform", adapt_method = "residual")$x
c_18c = admm_with_psd_mode_02(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, scale_mode = "none", adapt_method = "residual")$x
c_18d = admm_with_psd_mode_02(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, scale_mode = "normalize", adapt_method = "adaptive")$x
c_18e = admm_with_psd_mode_02(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, scale_mode = "transform", adapt_method = "adaptive")$x
c_18f = admm_with_psd_mode_02(A=P,Q=PT1,r=vT1,B_list=B_list,eps=eps13T, scale_mode = "none", adapt_method = "adaptive")$x



dC = rbind(
  data.frame(solid=1,basis=coef_name,idx=1:(na+nb),coef=c_1[1:(na+nb)],method="direct",st="H"),
  data.frame(solid=2,basis=coef_name,idx=1:(na+nb),coef=c_2[1:(na+nb)],method="direct",st="T"),
  data.frame(solid=3,basis=coef_name,idx=1:(na+nb),coef=c_3[1:(na+nb)],method="Schur complement",st="H"),
  data.frame(solid=4,basis=coef_name,idx=1:(na+nb),coef=c_4[1:(na+nb)],method="Schur complement",st="T"),
  data.frame(solid=5,basis=coef_name,idx=1:(na+nb),coef=c_5[1:(na+nb)],method="Null space",st="H"),
  data.frame(solid=6,basis=coef_name,idx=1:(na+nb),coef=c_6[1:(na+nb)],method="Null space",st="T"),
  data.frame(solid=7,basis=coef_name,idx=1:(na+nb),coef=c_7[1:(na+nb)],method="penalty",st="H"),
  data.frame(solid=8,basis=coef_name,idx=1:(na+nb),coef=c_8[1:(na+nb)],method="penalty",st="T"),
  data.frame(solid=9,basis=coef_name,idx=1:(na+nb),coef=c_9[1:(na+nb)],method="Null space",st="T,C"),
  data.frame(solid=10,basis=coef_name,idx=1:(na+nb),coef=c_11[1:(na+nb)],method="penalty",st="T,C"),
  data.frame(solid=11,basis=coef_name,idx=1:(na+nb),coef=c_12[1:(na+nb)],method="svd",st="T,C"),
  data.frame(solid=12,basis=coef_name,idx=1:(na+nb),coef=c_12b[1:(na+nb)],method="svd",st="H"),
  data.frame(solid=13,basis=coef_name,idx=1:(na+nb),coef=c_12c[1:(na+nb)],method="svd",st="T"),
  data.frame(solid=14,basis=coef_name,idx=1:(na+nb),coef=c_13[1:(na+nb)],method="eigen",st="T"),
  data.frame(solid=15,basis=coef_name,idx=1:(na+nb),coef=c_13b[1:(na+nb)],method="eigen",st="H"),
  data.frame(solid=16,basis=coef_name,idx=1:(na+nb),coef=c_14a[1:(na+nb)],method="ADMM0.1",st="T"),
  data.frame(solid=17,basis=coef_name,idx=1:(na+nb),coef=c_14b[1:(na+nb)],method="ADMM1.0",st="T"),
  data.frame(solid=18,basis=coef_name,idx=1:(na+nb),coef=c_14c[1:(na+nb)],method="ADMM10",st="T"),
  data.frame(solid=19,basis=coef_name,idx=1:(na+nb),coef=c_14d[1:(na+nb)],method="ADMM0.1",st="H"),
  data.frame(solid=20,basis=coef_name,idx=1:(na+nb),coef=c_14e[1:(na+nb)],method="ADMM1.0",st="H"),
  data.frame(solid=21,basis=coef_name,idx=1:(na+nb),coef=c_14f[1:(na+nb)],method="ADMM10",st="H"),
  data.frame(solid=22,basis=coef_name,idx=1:(na+nb),coef=c_15a[1:(na+nb)],method="ADMMPSD-0.1",st="T"),
  data.frame(solid=23,basis=coef_name,idx=1:(na+nb),coef=c_15b[1:(na+nb)],method="ADMMPSD-1.0",st="T"),
  data.frame(solid=24,basis=coef_name,idx=1:(na+nb),coef=c_15c[1:(na+nb)],method="ADMMPSD-10",st="T"),
  data.frame(solid=25,basis=coef_name,idx=1:(na+nb),coef=c_15d[1:(na+nb)],method="ADMMPSD-0.1",st="H"),
  data.frame(solid=26,basis=coef_name,idx=1:(na+nb),coef=c_15e[1:(na+nb)],method="ADMMPSD-1.0",st="H"),
  data.frame(solid=27,basis=coef_name,idx=1:(na+nb),coef=c_15f[1:(na+nb)],method="ADMMPSD-10",st="H"),
  data.frame(solid=28,basis=coef_name,idx=1:(na+nb),coef=c_16a[1:(na+nb)],method="ADMMPSDt0.1",st="T"),
  data.frame(solid=29,basis=coef_name,idx=1:(na+nb),coef=c_16b[1:(na+nb)],method="ADMMPSDt1.0",st="T"),
  data.frame(solid=30,basis=coef_name,idx=1:(na+nb),coef=c_16c[1:(na+nb)],method="ADMMPSDt10",st="T"),
  data.frame(solid=31,basis=coef_name,idx=1:(na+nb),coef=c_16d[1:(na+nb)],method="ADMMPSDn0.1",st="T"),
  data.frame(solid=32,basis=coef_name,idx=1:(na+nb),coef=c_16e[1:(na+nb)],method="ADMMPSDn1.0",st="T"),
  data.frame(solid=33,basis=coef_name,idx=1:(na+nb),coef=c_16f[1:(na+nb)],method="ADMMPSDn10",st="T"),
  data.frame(solid=34,basis=coef_name,idx=1:(na+nb),coef=c_17a[1:(na+nb)],method="ADMMPSDL2t0.1",st="T"),
  data.frame(solid=35,basis=coef_name,idx=1:(na+nb),coef=c_17b[1:(na+nb)],method="ADMMPSDL2t1.0",st="T"),
  data.frame(solid=36,basis=coef_name,idx=1:(na+nb),coef=c_17c[1:(na+nb)],method="ADMMPSDL2t10",st="T"),
  data.frame(solid=37,basis=coef_name,idx=1:(na+nb),coef=c_17d[1:(na+nb)],method="ADMMPSDL2n0.1",st="T"),
  data.frame(solid=38,basis=coef_name,idx=1:(na+nb),coef=c_17e[1:(na+nb)],method="ADMMPSDL2n1.0",st="T"),
  data.frame(solid=39,basis=coef_name,idx=1:(na+nb),coef=c_17f[1:(na+nb)],method="ADMMPSDL2n10",st="T"),
  data.frame(solid=40,basis=coef_name,idx=1:(na+nb),coef=c_17g[1:(na+nb)],method="ADMMPSDL2-0.1",st="T"),
  data.frame(solid=41,basis=coef_name,idx=1:(na+nb),coef=c_17h[1:(na+nb)],method="ADMMPSDL2-1.0",st="T"),
  data.frame(solid=42,basis=coef_name,idx=1:(na+nb),coef=c_17i[1:(na+nb)],method="ADMMPSDL2-10",st="T"),
  data.frame(solid=43,basis=coef_name,idx=1:(na+nb),coef=c_18a[1:(na+nb)],method="ADMMPSDnR",st="T"),
  data.frame(solid=44,basis=coef_name,idx=1:(na+nb),coef=c_18b[1:(na+nb)],method="ADMMPSDtR",st="T"),
  data.frame(solid=45,basis=coef_name,idx=1:(na+nb),coef=c_18c[1:(na+nb)],method="ADMMPSD-R",st="T"),
  data.frame(solid=46,basis=coef_name,idx=1:(na+nb),coef=c_18d[1:(na+nb)],method="ADMMPSDnA",st="T"),
  data.frame(solid=47,basis=coef_name,idx=1:(na+nb),coef=c_18e[1:(na+nb)],method="ADMMPSDtA",st="T"),
  data.frame(solid=48,basis=coef_name,idx=1:(na+nb),coef=c_18f[1:(na+nb)],method="ADMMPSD-A",st="T")
) %>% group_by(method,st) %>% mutate(scale=max(abs(coef)),pol=sign(mean(max(coef)+min(coef)))) %>% ungroup()

dn();dC %>% ggplot_bw(aes(x=idx,y=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point()+geom_line(aes(linetype=st))+facet_wrap(~method)

dn();dC %>% ggplot_bw(aes(y=basis,x=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point(alpha=0.5)+geom_path()+facet_wrap(~method,nrow=1)+scale_y_discrete(limits = coef_name)


dn();dC %>% separate(basis,into=c("aij","phi"),sep="_",remove = FALSE) %>% filter(aij=="a11") %>% 
  ggplot_bw(aes(x=basis,y=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point(alpha=0.5)+geom_line(aes(linetype=st))+geom_hline(yintercept = 0,linetype = "dotted")+
  facet_wrap(~method,scales="free_y",ncol=4)+theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0))

dn();dC %>% separate(basis,into=c("aij","phi"),sep="_",remove = FALSE) %>% filter(aij=="a12") %>% 
  ggplot_bw(aes(x=basis,y=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point(alpha=0.5)+geom_line(aes(linetype=st))+geom_hline(yintercept = 0,linetype = "dotted")+
  facet_wrap(~method,scales="free_y",ncol=4)+theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0))

dn();dC %>% separate(basis,into=c("aij","phi"),sep="_",remove = FALSE) %>% filter(aij=="a22") %>% 
  ggplot_bw(aes(x=basis,y=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point(alpha=0.5)+geom_line(aes(linetype=st))+geom_hline(yintercept = 0,linetype = "dotted")+
  facet_wrap(~method,scales="free_y",ncol=4)+theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0))

if(Mp>0){
  dn();dC %>% separate(basis,into=c("aij","phi"),sep="_",remove = FALSE) %>% filter(aij=="b") %>% 
    ggplot_bw(aes(x=basis,y=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
    geom_point(alpha=0.5)+geom_hline(yintercept = 0,linetype = "dotted")+
    facet_wrap(~method,scales="free_y")+theme(axis.text.x = element_text(angle = -45, hjust = 0))
}

dG = NULL

sollist = dC$solid %>% unique()
for(cc in sollist){
  
  c_n = dC %>% filter(solid==cc) %>% .$coef
  method_n = dC %>% filter(solid==cc) %>% .$method %>% unique()
  st_n = dC %>% filter(solid==cc) %>% .$st %>% unique()
  print(c(method_n,st_n))
  
  lam1 = (1:Tx)*0
  lam2 = (1:Tx)*0
  g11 = (1:Tx)*0
  g12 = (1:Tx)*0
  g22 = (1:Tx)*0
  for(tt in 1:Tx){
    for(ll in 1:Mm){
      philt = gphi[ll,tt]
      g11[tt] = g11[tt] + c_n[ idx[1,1,ll] ] * philt
      g12[tt] = g12[tt] + c_n[ idx[1,2,ll] ] * philt
      g22[tt] = g22[tt] + c_n[ idx[2,2,ll] ] * philt
    }
    g = matrix(c(g11[tt],g12[tt],g12[tt],g22[tt]),nrow=2)
    eg = eigen(g,symmetric = TRUE, only.values = TRUE)
    lam1[tt] = max(eg$values)
    lam2[tt] = min(eg$values)
  }
  
  dG0 = data.frame(t = t,trG = KT%*%c_n[1:(na+nb)] ,H = KH%*%c_n[1:(na+nb)] ,
                   e1 = P[1:Tx,]%*%c_n[1:(na+nb)] , e2 = P[Tx+1:Tx,]%*%c_n[1:(na+nb)] ,
                   lam1 = lam1, lam2=lam2,g11=g11,g12=g12,g22=g22,
                   method=method_n,st=st_n)
  
  dG = rbind(dG,dG0)
}


dG2 = dG %>% group_by(method,st) %>% 
  summarise(e_r=mean(sqrt(e1*e1+e2*e2))/trG[1], Hsd_r=sd(H)/trG[1], slam1=min(sign(lam1)), slam2=min(sign(lam2)) ) %>% 
  ungroup() %>% arrange(-slam1-slam2,e_r) 

print(dG2)
dG2 %>% view()
dn();dC %>% inner_join(y=dG2 %>% filter(slam1+slam2 ==2),by=c("method","st")) %>% 
  ggplot_bw(aes(y=basis,x=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point(alpha=0.5)+geom_path()+facet_wrap(~method,nrow=1)+scale_y_discrete(limits = coef_name)+
  geom_vline(xintercept = 0,linetype = "dotted")

# dn();colSums(P^2) %>% plot()


dn();  ggplot_bw()+geom_point(
  data= dG %>% inner_join(y=dG2 %>% filter(slam1+slam2 ==2),by=c("method","st")) %>% 
    filter((1000*t) %in% seq(from=1,to=20000,by=100)),
  aes(x=t,y=g11,shape=st,colour = st,group = paste(method,st)))+geom_path()+facet_wrap(~method,scale="free_y")+
  geom_line(mapping = aes(x=t,y=g11/(gain1*gain1)),data = df %>% filter((1000*t) %in% seq(from=1,to=20000,by=100)),alpha=0.5)

dn();  ggplot_bw()+geom_point(
  data= dG %>% inner_join(y=dG2 %>% filter(slam1+slam2 ==2),by=c("method","st")) %>% 
    filter((1000*t) %in% seq(from=1,to=20000,by=100)),
  aes(x=t,y=g12,shape=st,colour = st,group = paste(method,st)))+geom_path()+facet_wrap(~method,scale="free_y")+
  geom_line(mapping = aes(x=t,y=g12/(gain1*gain2)),data = df %>% filter((1000*t) %in% seq(from=1,to=20000,by=100)),alpha=0.5)

dn();  ggplot_bw()+geom_point(
  data= dG %>% inner_join(y=dG2 %>% filter(slam1+slam2 ==2),by=c("method","st")) %>% 
    filter((1000*t) %in% seq(from=1,to=20000,by=100)),
  aes(x=t,y=g22,shape=st,colour = st,group = paste(method,st)))+geom_path()+facet_wrap(~method,scale="free_y")+
  geom_line(mapping = aes(x=t,y=g22/(gain2*gain2)),data = df %>% filter((1000*t) %in% seq(from=1,to=20000,by=100)),alpha=0.5)

dn();  ggplot_bw()+geom_point(
  data= dG %>% inner_join(y=dG2 %>% filter(slam1+slam2 ==2),by=c("method","st")) %>% 
    filter((1000*t) %in% seq(from=1,to=20000,by=100)),
  aes(x=t,y=H,shape=st,colour = st,group = paste(method,st)))+geom_path()+facet_wrap(~method,scale="free_y")+
  geom_line(mapping = aes(x=t,y=H1),data = df %>% filter((1000*t) %in% seq(from=1,to=20000,by=100)),alpha=0.5)

