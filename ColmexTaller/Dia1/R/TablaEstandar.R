# Tim Riffe 14-06-2015
# Dia 1, tabla de vida, extensiones, interpretaciones.

####################################################################
# instalamos unos packages 
####################################################################

###############################################################
# todo lo que escribes ewn este lado es invisible al R
###############################################################

# regla 1) escribir apuntes directamente en el codigo

# R es una calculadora:"
a <- 1+1
a+ 1

# mas basica = vector
a <- 1:10
a <- 0:110

# objecto muy comun es data.frame:
a <- data.frame(col1 = 1:3, col2 = rnorm(3))
a

# cuales dimensiones tiene:
dim(a)
# buscar columnas por nombre:
a$col1
#a[,1]
#a[,"col1"]
#a[["col1"]]

# data.frame != matriz
b <- as.matrix(a)

t(b)
b * 3:1

# bucle basica-----
z <- NULL
for (i in 1:3){
	z[i] <- i ^ 2
}
#############################
# alguna funciones basicas de R:
sum(c(2,4,6))
cumsum(c(2,4,6))
prod(c(2,4,6))
cumprod(c(2,4,6))
##############################
# definimos una funcion:

MiFun <- function(x,y){
	
	# escribimos todo legible:
	(x ^ y) / exp(y)
	
	# la funcion es compartible
	
	# la funcion es permanente
	# trasferible dentro de codigo
}
MiFun(x,y)
##################

# devtools nos ayuda bajar packages que esten en github.com (lo que uso yo)
install.packages("devtools")
library(devtools)

# ahora instalamos un package que nos hace tablas de vida:
install_github("timriffe/LifeTable/LifeTable")
# otro package que nos vincula con los bases de datos HMD, HFD, y otros:
install.packages("XML")
install_github("timriffe/TR1/TR1/HMDHFDplus")
# otro package que nos hace superficies Lexis baratas
install_github("timriffe/LexisUtils/LexisUtils")
# buscamos el package MortalitySmooth de Giancarlo Camarda:
install.packages("MortalitySmooth")

# tal vez se puede usar este package:
#devtools::load_all("/home/tim/git/DistributionTTD/DistributionTTD/R/DistributionTTD")
#install_github("DistributionTTD/DistributionTTD/R/DistributionTTD")
# cargar los packages
library(HMDHFDplus)
library(LifeTable) # LT()
library(LexisUtils)      # por si acaso
library(MortalitySmooth)

############################################################
# cargar datos de defunciones y exposiciones de algun pais #
#----------------------------------------------------------#
# !! Registrar como usuario del HMD: www.mortality.org !!  #
############################################################

# ahora, estando registrado, puedes usar tu contraseña y nombre de usuario
# para bajar datos directo de la pagina web del HMD

# haz esto primero
us <- userInput() # esto es interactivo
pw <- userInput()

Pais <- "SWE"
Def <- readHMDweb(Pais,"Deaths_1x1",username = us, password = pw)
Exp <- readHMDweb(Pais,"Exposures_1x1",username = us, password = pw)

# habra que cambiar esto:
#setwd("/home/tim/git/ColmexTaller/ColmexTaller/Dia1")

# guardar para poder trabajar sin internet:
#save(Def, file = "Datos/Def.Rdata")
#save(Exp, file = "Datos/Exp.Rdata")
# Def <- local(get(load("Datos/Def.Rdata")))
# Exp <- local(get(load("Datos/Exp.Rdata")))
# hecharle un vistazo
head(Def)
head(Exp)

# elige hombres o mujers y algun año
sexo <- "Male"
año  <- 2000

# sabemos que va de la edad 0 a 110+
Dx <- Def[Def$Year == año, sexo]
Ex <- Exp[Exp$Year == año, sexo]

# nuestro Mx 
Mx <- Dx / Ex

# define algunos edades:
edades <- 0:110

# miramos que hay: bastante ruido, pero parece normalito
plot(edades, Mx, log = 'y', type = 'l', xlab = "edad", ylab = "tasa mortalidad",
		main = paste(Pais, sexo, año, sep = ", " ))

Ex

Ex[Ex == 0] <- .1
Dx[Dx == 0] <- .1
Tabla <- LT(Nx = Ex, Dx = Dx, sex = tolower(sexo), mxsmooth = FALSE)
names(Tabla)
head(Tabla$LT)


# miramos las columnas distintas:
plot(edades, Tabla$dx, type = 'l', main = paste("d(x)", Pais, sexo, año, sep = ", " ))
plot(edades, Tabla$lx, type = 'l', main = paste("l(x)", Pais, sexo, año, sep = ", " ))
plot(edades, Tabla$qx, type = 'l', main = paste("q(x)", Pais, sexo, año, sep = ", " ),log='y')

# e(x) es bastante lineal al principio, pero luego no
plot(edades, Tabla$ex, type = 'l', main = paste("e(x)", Pais, sexo, año, sep = ", " ), 
		asp = 1,xlim=c(0,110),ylim=c(0,110),xaxs="i",yaxs='i',lwd=2)
abline(a=Tabla$ex[1],b=-1,col="red")


# vamos a usar dx para calcular todas las formulas en la lectura, y para eso
# hace falta definir un monton de funciones aqui:

# un par de utilidades:
Mna0 <- function(M){
	M[is.na(M)]  <- 0
	M[is.nan(M)] <- 0
	M
}
Minf0 <- function(M){
	M[is.infinite(M)]  <- 0
	M
}

# esta funcion hace lo de la equacion (2) de la lectura. 
# se puede replicar de manera mas intuitivo si se quiere
da2fya <- function(da, stagger = FALSE){
	N       <- length(da)
	ay      <- 1:N - 1
	
	da      <- Mna0(da)   # remove NAs if any       
	da      <- c(da, da * 0) / sum(da) # pad out with 0s
	fya     <- matrix(da[col(matrix(nrow = N, 
									ncol = N)) + ay], 
			nrow = N, 
			ncol = N, 
			dimnames = list(Ex = ay, 
					Age = ay)
	)
	if (stagger){
		fya <- (fya + cbind(fya[, 2:ncol(fya)], 0)) / 2
	}
	fya <- Minf0(Mna0(fya / rowSums(fya)))
	fya
}

# otra manera de calcular e(x), lo de eqn (3)
getex <- function(dx, ax = rep(.5, length(dx))){
	fya <- da2fya(dx)
	rowSums((col(fya) - (1 - ax)) * fya)
}

# funcion de momentos, equacion (4) de la lectura
momentN <- function(dx, n, ax = rep(.5, length(dx))){
	fya <- da2fya(dx)
	ex  <- rowSums((col(fya) - (1 - ax)) * fya)
	rowSums(((col(fya) - .5) - ex)^n * fya)
}

# funcion de asimetria, eqn (5)
getSkewst <- function(dx, ax){
	# formula from http://en.wikipedia.org/wiki/Skewness
	momentN(dx, n = 3, ax = ax) / (momentN(dx, n = 2, ax = ax) ^ (3/2)) 
}

# funcion de curtosis, eqn (6)
getKurtst <- function(dx, ax){
	# formula from http://en.wikipedia.org/wiki/Kurtosis
	momentN(dx, n = 4, ax = ax) / (momentN(dx, n = 2, ax = ax) ^ 2) - 3
}

# tambien se puede calcular quantiles sencillos sobre la distribucion
# condicionada de f(y|x) ... no lo he escrito mucho en la lectura, 
# pero esto es igual de util o mas
getQuantile <- function(dx, probs = c(.1,.25,.5,.75,.9)){
	fya     <- da2fya(dx)
	# needs age at clean breaks, like lx
	CDF     <- cbind(0, t(apply(fya, 1, cumsum)))
	Age     <- 1:ncol(CDF) - 1
	# monotonic avoids negatives. Last value zero for ease, but not comparable with e_\omega
	Quantiles <- sapply(probs, function(p, CDF, Age){
				medAges <- apply(CDF, 1, function(cdf, Age, p){
							if (length(unique(cdf)) > 2){
								return(splinefun(Age~cdf, method = "monoH.FC")(p))
							} else {
								# default value, didn't think much on this
								return(0.5)
							}
						}, Age = Age, p = p)     
			}, CDF = CDF, Age = Age)
	Quantiles[Quantiles < 0] <- 0 
	if (length(probs) == 1){
		return(Quantiles)
	}
	colnames(Quantiles) <- probs
	Quantiles 
}

##############################
logit <- function(x){
	log(x/(1+x))
}

# una idea general de que forma tiene f(y|x):
image(0:110,0:110,logit(da2fya(Tabla$dx)),xlab="chrono",ylab = "thano")
contour(0:110,0:110,logit(da2fya(Tabla$dx)),add=TRUE,n=10)

# pero tal vez mas intuitivo seria:
matplot(0:110,t(da2fya(Tabla$dx)[c(seq(1,101,by=10)),]),type='l',
		lty=1,col=gray(seq(0,.5,length=11)),lwd=seq(1,2,length=11),
		main = "f(y|x), edades seleccionadas", xlab = "años restantes",ylab="f(y|x)")
# f(y|x) es la distribucion que queremos resumir de varias formas por edad.
# como la distribucion de y varia POR edad, resulatara de qualquier indice que
# usamos para resumirlo tendra su propia pauta por edad. La question es si estas
# pautas por edad nos resultan utiles o no. Si son de utilidad, probablemente sera
# desde el perspectivo individual mirando hacia el futuro.

###############################################################
# Varianza y desviacion estandar [sigma^2(y|x) y sigma(y|x)]  #
###############################################################
plot(edades, momentN(Tabla$dx,2,Tabla$ax), type = 'l', xlab = "edad", ylab = "Varianza",
		main = "Varianza de duracion restante de vida")
plot(edades, sqrt(momentN(Tabla$dx,2,Tabla$ax)), type = 'l', xlab = "edad", ylab = "DS",
		main = "Desviacion estandar de duracion restante de vida")

###############################################################
# Skewness = asimetria de f(y|x)                              #
###############################################################
plot(edades, getSkewst(Tabla$dx,Tabla$ax), type = 'l', xlab = "edad", ylab = "Asimetria",
		main = "Asimetria de duracion restante de vida", ylim = c(-2,5),xlim=c(0,100))
abline(h=0,col = "red")
text(20,-.5,"f(y|x) tira a edades jovenes")
text(80,2,"f(y|x) tira a edades mayores")

###############################################################
# Kurtosis = curtosis ~= apuntamiento de f(y|x)               #
###############################################################
plot(edades, getKurtst(Tabla$dx,Tabla$ax), type = 'l', xlab = "edad", ylab = "apuntamiento",
		main = "apuntamiento de duracion restante de vida", ylim=c(-2,7),xlim=c(0,100))
abline(h=0,col = "red")

###############################################################
# Kurtosis = curtosis ~= apuntamiento de f(y|x)               #
###############################################################
CV <-  sqrt(momentN(Tabla$dx,2,Tabla$ax)) / getex(Tabla$dx,Tabla$ax)
plot(edades, CV, type = 'l', xlab = "edad", ylab = "CV",
		main = "coeficiente de variacion de duracion restante de vida", xlim=c(0,100))
text(25,.9,"menos informacion = e(x) menos util...")
text(25,.1,"mas informacion = e(x) util")
# este es mi favorito, porque es totalmente comparable sobre la edad:
# piensa: si tienes 10 años y te dicen que e(10) = 60 +/- 5 es mucho mas informativo
# que si tienes 80 años y te dicen que e(80) = 5 +/- 5 (valores inventados...)
