dados = converte(read.table("subsfin.txt"))	
locs = converte(read.table("mote_locs.txt")[,2:3])
locs=locs[c(-5,-15),]

#++++++++++++++++++++++++++++++++++++++++

sensores = matrix(1, ncol = ncol(locs), nrow = (nrow(locs)+100))

sensores[1:nrow(locs), 1:ncol(locs)] = locs

for(i in (nrow(locs)+1):nrow(sensores)){
	sensores[i,1] = sample(min(locs[,1]):max(locs[,1]), 1) 
	sensores[i,2] = sample(min(locs[,2]):max(locs[,2]), 1) 
}

dados_rgl = complete_rgl(sensores, dados, locs)
dados_rbf = complete_rbf(sensores, dados, locs)
dados_smo = complete_smo(sensores, dados, locs)
dados_krig = complete_krig(sensores, dados, locs)


	
result_regLin2 = eco_p(0.2, 1, dados_rgl)
result_regLin3 = eco_p(0.3, 1, dados_rgl)
result_regLin4 = eco_p(0.4, 1, dados_rgl)
result_regLin5 = eco_p(0.5, 1, dados_rgl)
result_regLin6 = eco_p(0.6, 1, dados_rgl)
result_regLin7 = eco_p(0.7, 1, dados_rgl)
result_regLin8 = eco_p(0.8, 1, dados_rgl)

result_rbf2 = eco_p(0.2, 2,sensores, dados_rbf)
result_rbf3 = eco_p(0.3, 2,sensores, dados_rbf)
result_rbf4 = eco_p(0.4, 2,sensores, dados_rbf)
result_rbf5 = eco_p(0.5, 2,sensores, dados_rbf)
result_rbf6 = eco_p(0.6, 2,sensores, dados_rbf)
result_rbf7 = eco_p(0.7, 2, sensores,dados_rbf)
result_rbf8 = eco_p(0.8, 2, sensores,dados_rbf)
print(1)

#++++++++++++++++++PREENCHENDO OS 100 PONTOS ALEATORIOS++++++++

complete_rgl = function(sensores, dados, locs){

	result = matrix(0, nrow = nrow(dados), ncol = nrow(sensores))

	result[1:nrow(dados), 1:ncol(dados)] = dados

	X = matrix(1, ncol = 3, nrow = nrow(locs))
	X[,2] = locs[,1]
	X[,3] = locs[,2]

	for(i in 1:nrow(result)){
		for(j in (nrow(locs)+1):nrow(sensores)){
					
			v = solve(t(X)%*%X)%*%t(X)%*%dados[i,]
			result[i,j] = (v[1,1] + v[2,1]*sensores[j,1] + v[3,1]*sensores[j,2])
		}
	}
	
	return(result)
}


complete_rbf = function(sensores, dados, locs){
	#result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	result = matrix(0, nrow = nrow(dados), ncol = nrow(sensores))
	result[1:nrow(dados), 1:ncol(dados)] = dados

	H = matrix(1, ncol = nrow(sensores)+1, nrow = nrow(sensores))
	s = 0.01

	for(i in 1:nrow(H)){
		for (j in 1:nrow(H)) {
			H[i,j+1] = exp(-s*((sensores[j,1] - sensores[i,1])^2 + (sensores[j,2] - sensores[i,2])^2)) 
		}
	}

	for(i in 1:nrow(result)){
		for(j in (nrow(locs)+1):nrow(sensores)){
			Ha = H[1:nrow(locs), 1:(nrow(locs)+1)]
			alpha= MASS::ginv(t(Ha)%*%Ha)%*%t(Ha)%*%dados[i,]
			Hj = H[j, 1:53]
			result[i,j]= Hj %*% alpha
		}
		print(i*j)
	}

	return(result)
}

complete_kri = function(sensores, dados, locs){
	#library("SpatialExtremes")
	#result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	result = matrix(0, nrow = nrow(dados), ncol = nrow(sensores))
	result[1:nrow(dados), 1:ncol(dados)] = dados

	for(i in 1:nrow(result)){
		for(j in (nrow(locs)+1):nrow(sensores)){	
			krig = kriging(dados[i,], locs, t(sensores[j,]),
			 cov.mod = "powexp", sill = 1, range = 10, smooth = 0.75)
			result[i,j] = krig$krig.est
		}
		print(i*j)
	}

	return(result)
}

resultado_smooth = function(sensores, dados, locs){
	dst = matrix(0, nrow = nrow(sensores), ncol = nrow(sensores))
	dst = dist(sensores, diag=TRUE, upper=TRUE)
	dst = as.matrix(dst)
	#result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	result = matrix(0, nrow = nrow(dados), ncol = nrow(sensores))
	result[1:nrow(dados), 1:ncol(dados)] = dados

	for(i in 1:nrow(result)){
		for(j in (nrow(locs)+1):nrow(sensores)){	

			result[i,j] = sum(dst[j,1:52]*dados[i,]/sum(dst[j,1:52]))
		}
	}

	return (result)
}


#++++++++++++++FUNÇÕES+++++++++++++++++++++++
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

eco_p = function(p, type, sensores, dados){
	result = matrix(0, nrow = nrow(dados), ncol = ncol(dados))
	for(i in 1:nrow(dados)){
		dt_epoca = c()
		count = 0
		for(j in 1:ncol(dados)){
			if(length(dt_epoca) > 148){
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
	X = matrix(1, ncol = 3, nrow = nrow(locs))
	X[,2] = locs[,1]
	X[,3] = locs[,2]

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
		print(i*j)
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