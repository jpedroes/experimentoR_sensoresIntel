dados = converte(read.table("subsfin.txt"))	
locs = converte(read.table("mote_locs.txt")[,2:3])
locs=locs[c(-5,-15),]

#result_regLin1 = eco_p(0.1, 1)
result_regLin2 = eco_p(0.2, 1)
result_regLin3 = eco_p(0.3, 1)
result_regLin4 = eco_p(0.4, 1)
result_regLin5 = eco_p(0.5, 1)
result_regLin6 = eco_p(0.6, 1)
result_regLin7 = eco_p(0.7, 1)
result_regLin8 = eco_p(0.8, 1)
#result_regLin9 = eco_p(0.9, 1)
#result_regLin = eco_p(1, 1)

plot(erromedio(result_regLin2, dados), type="l", col="blue", ylim=c(0,01))
lines(erromedio(result_regLin3, dados), type="l", pch=22, lty=2, col="yellow")
lines(erromedio(result_regLin4, dados), type="l", pch=22, lty=2, col="black")
lines(erromedio(result_regLin5, dados), type="l", pch=22, lty=2, col="red")
lines(erromedio(result_regLin6, dados), type="l", pch=22, lty=2, col="green")
lines(erromedio(result_regLin7, dados), type="l", pch=22, lty=2, col="pink")
lines(erromedio(result_regLin8, dados), type="l", pch=22, lty=2, col="orange")

result_kri2 = eco_p(0.2, 4)
result_kir3 = eco_p(0.3, 4)
result_kri4 = eco_p(0.4, 4)
result_kri5 = eco_p(0.5, 4)
result_kri6 = eco_p(0.6, 4)
result_kri7 = eco_p(0.7, 4)
result_kri8 = eco_p(0.8, 4)
print(1)


result_rbf2 = eco_p(0.2, 2)
result_rbf3 = eco_p(0.3, 2)
result_rbf4 = eco_p(0.4, 2)
result_rbf5 = eco_p(0.5, 2)
result_rbf6 = eco_p(0.6, 2)
result_rbf7 = eco_p(0.7, 2)
result_rbf8 = eco_p(0.8, 2)
print(1)

#==============================================================

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

eco_p = function(p, type){
	dados = converte(read.table("subsfin.txt"))	
	sensores = converte(read.table("mote_locs.txt")[,2:3])
	sensores=sensores[c(-5,-15),]
	result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	for(i in 1:nrow(dados)){
		dt_epoca = c()
		count = 0
		for(j in 1:ncol(dados)){
			if(length(dt_epoca) > 48){
				j = ncol(dados)+1
			}
			if(runif(1) < p){
				count = count + 1
				dt_epoca[count] = j
			}
			if(length(dt_epoca)<1){
				dt_epoca[1] = sample(1:52, 1)
			}

		}
		if(type == 1){
			result[i,] = regLin_p(dt_epoca, sensores, dados[i,])
			
		}
		if(type == 2){
			result[i,] = regRbf_p(dt_epoca, sensores, dados[i,])
		}
		if(type == 3){
			result[i,] = regSmo_p(dt_epoca, sensores, dados[i,])	
		}
		if(type == 4){
			result[i,] = regKri_p(dt_epoca, sensores, dados[i,])
		}
	}
	return(result)
}

regLin_p = function(dt_epoca, sensores, dados){

	result = array(0, dim = length(dados))
	X = matrix(1, ncol = 3, nrow = nrow(sensores))
	X[,2] = sensores[,1]
	X[,3] = sensores[,2]

	if(length(dt_epoca) == length(dados)){
		return (dados)
	}
	else{
		for(i in 1:length(result)){
			if(i %in% dt_epoca){
				result[i] = dados[i]
			}
			else{
				x = X[-dt_epoca,]		
				v = solve(t(x)%*%x)%*%t(x)%*%dados[-dt_epoca]
				result[i] = (v[1,1] + v[2,1]*sensores[i,1] + v[3,1]*sensores[i,2])
			}
		}
	}	
	return(result)
}

regRbf_p = function(dt_epoca, sensores, dados){
	#result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	result = array(0, dim = length(dados))
	H = matrix(1, ncol = nrow(sensores)+1, nrow = nrow(sensores))
	s = 0.01

	for(i in 1:nrow(H)){
		for (j in 1:nrow(H)) {
			H[i,j+1] = exp(-s*((sensores[j,1] - sensores[i,1])^2 + (sensores[j,2] - sensores[i,2])^2)) 
		}
	}

	for(i in 1:length(result)){
		if(i %in% dt_epoca){
			result[i] = dados[i]
		}
		else{
			Ha= H[-dt_epoca,-dt_epoca]
			alphaj= MASS::ginv(t(Ha)%*%Ha)%*%t(Ha)%*%dados[-dt_epoca]
			nt = dt_epoca + 1
			Hj = H[i, -(nt)]
			result[i]= Hj %*% alphaj
		}
	}

	return(result)
}

regSmo_p = function(dt_epoca, sensores, dados){
	dst = matrix(0, nrow = nrow(sensores), ncol = nrow(sensores))
	dst = dist(sensores, diag=TRUE, upper=TRUE)
	dst = as.matrix(dst)
	#result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	result = array(0, dim = length(dados))

	for(i in 1:length(result)){
		if(i %in% dt_epoca){
			result[i] = dados[i]
		}
		else{
			result[i] = sum(dst[i,-dt_epoca]*dados[-dt_epoca]/sum(dst[i,-dt_epoca]))
		}
	}

	return (result)
}

regKri_p = function(dt_epoca, sensores, dados){
	#library("SpatialExtremes")
	#result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	result = array(0, dim = length(dados))

	

	for(i in 1:length(result)){
		if(i %in% dt_epoca){
			result[i] = dados[i]
		}
		else{
			krig = kriging(dados[-dt_epoca], sensores[-dt_epoca,], t(sensores[i,]),
			 cov.mod = "powexp", sill = 1, range = 10, smooth = 0.75)
			result[i] = krig$krig.est
		}
	}

	return(result)
}

