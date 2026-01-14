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
df = read.csv("model2_g_1_0_-1_1_sol.csv",header = TRUE,stringsAsFactors = FALSE)

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
x0 = (df$x1-min(df$x1))/(max(df$x1)-min(df$x1)) 
y0 = (df$x2-min(df$x2))/(max(df$x2)-min(df$x2)) 

const = 2*(0-min(df$x1))/(max(df$x1)-min(df$x1)) - 1
gain1 = 2/(max(df$x1)-min(df$x1)) 
gain2 = 2/(max(df$x2)-min(df$x2)) 



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

# phi = phiC
# phid = phiCd
# phi_name = phiC_name
phi = bnst[[2]]
phid = bnstd[[2]]
phi_name = bnst_name[[2]]

psi = phiB
psid = phiBd
psi_name = phiB_name

Mm = 9 #6 10 15 21 28
Mp = 2 # 2 5 9 14
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
for(kk in 1:Mm){
  B_list[[kk]] = matrix(0,ncol=na+nb,nrow=3)
  B_list[[kk]][1,idx[1,1,kk]] = 1
  B_list[[kk]][2,idx[1,2,kk]] = 1
  B_list[[kk]][3,idx[2,2,kk]] = 1
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

# trace = 1, trace
PT1 = matrix(0,nrow=1, ncol=na+nb)
KT = matrix(0,nrow=Tx, ncol=na+nb)
vT1 = matrix(1,nrow=1, ncol=1)
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

# H = 1
PH1 = matrix(0,nrow=1, ncol=na+nb)
KH = matrix(0,nrow=Tx, ncol=na+nb)
vH1 = matrix(1,nrow=1, ncol=1)
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
  data.frame(solid=22,basis=coef_name,idx=1:(na+nb),coef=c_15a[1:(na+nb)],method="ADMMPSD0.1",st="T"),
  data.frame(solid=23,basis=coef_name,idx=1:(na+nb),coef=c_15b[1:(na+nb)],method="ADMMPSD1.0",st="T"),
  data.frame(solid=24,basis=coef_name,idx=1:(na+nb),coef=c_15c[1:(na+nb)],method="ADMMPSD10",st="T"),
  data.frame(solid=25,basis=coef_name,idx=1:(na+nb),coef=c_15d[1:(na+nb)],method="ADMMPSD0.1",st="H"),
  data.frame(solid=26,basis=coef_name,idx=1:(na+nb),coef=c_15e[1:(na+nb)],method="ADMMPSD1.0",st="H"),
  data.frame(solid=27,basis=coef_name,idx=1:(na+nb),coef=c_15f[1:(na+nb)],method="ADMMPSD10",st="H"),
  data.frame(solid=28,basis=coef_name,idx=1:(na+nb),coef=c_16a[1:(na+nb)],method="ADMMPSDt0.1",st="T"),
  data.frame(solid=29,basis=coef_name,idx=1:(na+nb),coef=c_16b[1:(na+nb)],method="ADMMPSDt1.0",st="T"),
  data.frame(solid=30,basis=coef_name,idx=1:(na+nb),coef=c_16c[1:(na+nb)],method="ADMMPSDt10",st="T"),
  data.frame(solid=31,basis=coef_name,idx=1:(na+nb),coef=c_16d[1:(na+nb)],method="ADMMPSDn0.1",st="T"),
  data.frame(solid=32,basis=coef_name,idx=1:(na+nb),coef=c_16e[1:(na+nb)],method="ADMMPSDn1.0",st="T"),
  data.frame(solid=33,basis=coef_name,idx=1:(na+nb),coef=c_16f[1:(na+nb)],method="ADMMPSDn10",st="T")
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
  for(tt in 1:Tx){
    g = matrix(c(0,0,0,0),nrow=2)
    for(ll in 1:Mm){
      philt = gphi[ll,tt]
      g11 = c_n[ idx[1,1,ll] ] * philt
      g21 = c_n[ idx[2,1,ll] ] * philt
      g12 = g21
      g22 = c_n[ idx[2,2,ll] ] * philt
      
      g = g + matrix(c(g11,g21,g12,g22),nrow=2)
    }
    eg = eigen(g,symmetric = TRUE, only.values = TRUE)
    lam1[tt] = max(eg$values)
    lam2[tt] = min(eg$values)
  }
  
  dG0 = data.frame(t = t,trG = KT%*%c_n[1:(na+nb)] ,H = KH%*%c_n[1:(na+nb)] ,
                   e1 = P[1:Tx,]%*%c_n[1:(na+nb)] , e2 = P[Tx+1:Tx,]%*%c_n[1:(na+nb)] ,
                   lam1 = lam1, lam2=lam2,
                   method=method_n,st=st_n)
  
  dG = rbind(dG,dG0)
}


dG2 = dG %>% group_by(method,st) %>% 
  summarise(e_r=mean(sqrt(e1*e1+e2*e2))/trG[1], Hsd_r=sd(H)/trG[1], slam1=min(sign(lam1)), slam2=min(sign(lam2)) ) %>% 
  ungroup() %>% arrange(-slam1-slam2,e_r) 

print(dG2)

dn();dC %>% inner_join(y=dG2 %>% filter(slam1+slam2 ==2),by=c("method","st")) %>% 
  ggplot_bw(aes(y=basis,x=coef/scale/pol,shape=st,colour = st,group = paste(method,st)))+
  geom_point(alpha=0.5)+geom_path()+facet_wrap(~method,nrow=1)+scale_y_discrete(limits = coef_name)

# dn();colSums(P^2) %>% plot()
