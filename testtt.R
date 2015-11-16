whos <- function(){
  l <- ls(env=parent.env(environment()))
  str0 <- 1:length(l)
  eval(parse(text=paste("w__d <- data.frame(matrix(0,ncol=4,nrow=",
                        length(l),"))",sep="")),env=parent.env(environment()))
  for (ii in 1:length(l)){
    #str0[ii] <- paste("w__d<-data.frame(name=",l[1],",class=class(",l[1],"), mode=mode(",l[1],"), length=length(",l[1],"))", sep="")
    str0[ii] <- paste("w__d[",ii,", ]<-c(ls()[",ii,"],class(",l[ii],"), mode(",l[ii],"), length(",l[ii],"))", sep="")
  }
  eval(parse(text=str0),env=parent.env(environment()))
  
  eval(parse(text="names(w__d) <- c(\"name\", \"class\", \"mode\", \"length\")"),env=parent.env(environment()))
  #eval(parse(text=str1))
  #eval(parse(text=str2))
  #eval(parse(text=str3))
  #print(parent.env(environment()))
  
  eval(parse(text=paste("head(w__d, n=",length(l),")",sep="")),env=parent.env(environment()))
  #eval(parse(text="rm(w__d)"),env=parent.env(environment()))
}


