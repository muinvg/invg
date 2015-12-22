install.packages("sankey")
library("sankey")

n<-15

V1 <- sample(letters[1:24],n,replace=T)
V2 <- sample(letters[1:24],n,replace=T)

G=NULL;for (i in seq(V1)){ G = paste(G,V1[i],V2[i],"\n",sep=" ")}

e <- read.table(stringsAsFactors=F,textConnection(G))

PED <- make_sankey(edges=u)

sankey(PED)

cat(G)



S2e <- function(S){
	eS <- NULL
	for( jj in seq(S)){
		aS <- rev(strsplit(gsub(")","",S[jj]),sp="\\(")[[1]])
		
		ran <- 1:(length(aS)-1)
		if(is.null(eS)){
			eS <- data.frame(V1=I(aS[ran]), V2=I(aS[ran+1]) )
		}else{
			eS <- rbind(eS,data.frame(V1=aS[ran], V2=aS[ran+1] ) )
		}
	
		
	}
	return(eS)
}

e1=S2e(c("max(abs(HI))","mean(abs(HI))"))
