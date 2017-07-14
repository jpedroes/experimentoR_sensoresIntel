dados = converte(read.table("subsfin.txt"))	
locs = converte(read.table("mote_locs.txt")[,2:3])
locs=locs[c(-5,-15),]


#q1
erromedio = function(result, dados){
	rmse_epoca = c()
	re2_min_ep = c()
	re2_max_ep = c()
	for(i in 1:nrow(result)){
		a = (dados[i,] - result[i,])^2
		rmse_epoca[i] = sqrt(sum(a))/ncol(result)
	}
	return (rmse_epoca)
}

erromin = function(result, dados){
	rmse_epoca = c()
	re2_min_ep = c()
	re2_max_ep = c()
	for(i in 1:nrow(result)){
		a = (dados[i,] - result[i,])^2
		re2_min_ep[i] = min(sqrt(a))
	}
	return (re2_min_ep)
}

erromax = function(result, dados){
	re2_max_ep = c()
	for(i in 1:nrow(result)){
		a = (dados[i,] - result[i,])^2
		re2_max_ep[i] = max(sqrt(a))
	}
	return (re2_max_ep)
}

converte = function(dados) {
 	c<-NULL
 	for (i in 1:ncol(dados))
		c<-cbind(c,as.numeric(dados[,i]))

	return (c)
}


result_regLin = resultado_regLin(locs, dados)
result_regRBF = resultado_rbf(locs, dados)
result_regSmooth = resultado_smooth(locs, dados)
result_krig = resultado_krigen(locs, dados)

#erros lin
plot(erromedio(result_regLin, dados), type="l", col="blue", ylim=c(0,25))

lines(erromin(result_regLin, dados), type="l", pch=22, lty=2, col="red")

lines(erromax(result_regLin, dados), type="l", pch=22, lty=2, col="green")

#erros rbf

plot(erromedio(result_regRBF, dados), type="l", col="blue", ylim=c(0,30))

lines(erromin(result_regRBF, dados), type="l", pch=22, lty=2, col="red")

lines(erromax(result_regRBF, dados), type="l", pch=22, lty=2, col="green")

#erros smooth

plot(erromedio(result_regSmooth, dados), type="l", col="blue", ylim=c(0,25))

lines(erromin(result_regSmooth, dados), type="l", pch=22, lty=2, col="red")

lines(erromax(result_regSmooth, dados), type="l", pch=22, lty=2, col="green")

#erros krig
plot(erromedio(result_krig, dados), type="l", col="blue", ylim=c(0,25))

lines(erromin(result_krig, dados), type="l", pch=22, lty=2, col="red")

lines(erromax(result_krig, dados), type="l", pch=22, lty=2, col="green")


resultado_regLin = function(sensores ,dados){
	result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	X = matrix(1, ncol = 3, nrow = nrow(locs))
	X[,2] = locs[,1]
	X[,3] = locs[,2]
	for(i in 1:nrow(dados)){
		for(j in 1:nrow(sensores)){
			x = X[-j,]		
			v = solve(t(x)%*%x)%*%t(x)%*%dados[i,-j]
			result[i,j] = (v[1,1] + v[2,1]*sensores[j,1] + v[3,1]*sensores[j,2])
		}
	}
	return(result)
}

resultado_rbf = function(sensores ,dados){
	result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	H = matrix(1, ncol = nrow(sensores)+1, nrow = nrow(sensores))
	#H = array(1, dim= nrow(locs)-1)
	s = 0.01

	for(i in 1:nrow(H)){
		for (j in 1:nrow(H)) {
			H[i,j+1] = exp(-s*((locs[j,1] - locs[i,1])^2 + (locs[j,2] - locs[i,2])^2)) 
		}
	}

	for(i in 1:nrow(dados)){
		for(j in 1:nrow(sensores)){			
			Ha= H[-(j),-(j+1)]
			alphaj= MASS::ginv(t(Ha)%*%Ha)%*%t(Ha)%*%dados[i,-j]
			Hj = H[j,-(j+1)]
			result[i,j] = Hj %*% alphaj
		}	
		print(i*j)
	}	
	return(result)
}


resultado_krigen = function(sensores, dados){
	library("SpatialExtremes")

	result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))

	for(i in 1:nrow(dados)){
		krig = kriging(dados[i,], sensores, t(sensores[j,]), cov.mod = "powexp", sill = 1, range = 10, smooth = 0.75)
	}

	for(i in 1:nrow(dados)){
		for(j in 1:nrow(sensores)){
			krig = kriging(dados[i,-j], sensores[-j,], t(sensores[j,]), cov.mod = "powexp", sill = 1, range = 10, smooth = 0.75)
			result[i,j] = krig$krig.est
			
		}
		print(i*j)
	}

	return(result)
}

