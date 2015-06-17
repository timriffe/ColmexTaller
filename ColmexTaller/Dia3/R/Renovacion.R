
Nac    <- dget("HFDSWE.txt")

##############################
# vamos a hacer Lotka en R!
##############################

# 1) necesitamos datos mortalidad y fecundidad
library(HMDHFDplus)
# 2) para buenas tablas de vida
library(LifeTable) # si no lo tienes, mira en el codigo de lunes
##############################
# si hace falta, define username & password
us <- userInput()
pw <- userInput()
Pais <- "SWE"
Def <- readHMDweb(Pais,"Deaths_1x1",us,pw)
Exp <- readHMDweb(Pais,"Exposures_1x1",us,pw)

# Tendria que funcionar, pero no me carga la pagina web hoy, no se porque...
#Nac <- readHFDweb(Pais,"birthsRR",us,pw)
# plan B: cargar datos localmente (mira en carpeta /Dia3/datos/


# los datos de HMD y HFD son de dimensiones distintas

# tienen años distintos
range(Nac$Year)
range(Def$Year)

# y los nacimientos de HFD no estan clasificados por sexo:
# y recuerde que queremos hijas nacidas a madres.
# asi que tendremos que dividir los nacimientos segun algun suposito
# lo normal es 1.05 ... 
# relacion de sexo al nacer: (razon de masculinidad)
RSN <- 1.05
# proporcion feminina:
(PF <- 1 / (1 + RSN)) # si lo pones en parentesis se ve el resultado tb
# multiplicamos los nacimientos para tener hijas solo.
# NB: se podria tener el RSN agregado desde el HMD de hecho. Verias
# que se cambia un poco cada año y que a veces parece haber tendencias.
# Seguimos con el 1.05 ...
##########################################################################
Nac$Births <- Nac$Births * PF
#Nac$Total <- Nac$Total * PF # claro, solo haz esto una vez
##########################################################################

# herramientas para seleccionar datos:
sexo <- "Female"
año  <- 1900

# buscamos datos de 1900 porque si:
Dx <- Def[Def$Year == año, sexo]
Ex <- Exp[Exp$Year == año, sexo]
# N son nacimientos para nosotros, pero nota que muchas veces en
# la literatura se usa para estocs de poblacion (Px)

Nx <- Nac$Births[Nac$Year == año][1:111]
# Nx <- Nac$Total[Nac$Year == año] # para datos bajados del HFD directamente
# los Nx van de la edad 12 al 55 solo.

#### necesario si hemos bsucado datos del HFD, pero los hemos copiado localmente
# añadimos 0s antes y despues para alinear edades
#Nx <- c(rep(0,12),Nx,rep(0,55))
# calculamos tasas de fecundidad:
fx <- Nx / Ex

# limpieza, sino tenemos casos de NaN:
fx[Ex == 0] <- 0

#plot(0:110, fx, type = 'l')

# vector de edades:
edades <- 0:110

# calculamos Lx, etc
Dx[Dx == 0] <- .1
Ex[Ex == 0] <- .1
Tabla <- LT(Ex,Dx, mxsmooth=FALSE)

Lx <- Tabla$Lx
dx <- Tabla$dx

# miramos lx por ver que tenemos
plot(edades, Tabla$lx, type = 'l', main = paste("supervivencia",Pais,sexo,año))

######################################################
# R0 de Lotka:
(R0 <- sum(Lx*fx))

# interpretacion: es el ISF feminino con un descuento 
# para la mortalidad de las madres. No garantiza que 
# las hijas llegan a edades reproductivas, solo que nacen.
# en este caso es > 1, y tenemos una poblacion que crece
# de forma natural.

######################################################
# la equacion Lotka parece a esto:
# 1 == sum(exp(-x*r)*Lx*fx)
# pero tenemos que buscar el 'r', la tasa intrinsica de crecimiento
#########
# una funcion que podemos optimizar:

LotkaMin <- function(r, Lx, fx, x = .5:110.5){
	abs(1 - sum(exp(-x*r)*Lx*fx))
}
# positivo, y grande en mi opnion:
(r <- optimize(LotkaMin, c(-.05,.05), Lx = Lx, fx = fx, tol = 1e-15)$minimum)
# comprobamos:
1 - sum(exp(-(edades+.5)*r)*Lx*fx) # un numero como 2e-9 es R intentando decir 0
                                   # quiere decir que 'r' esta bien.

# ya con este r sabemos la estructura estable de la poblacion:

plot(edades, exp(-(edades+.5)*r)*Lx , type = 'l')
cx <- exp(-(edades+.5)*r)*Lx / sum(exp(-(edades+.5)*r)*Lx)
plot(edades, cx, type = 'l')
# compare con el Lx tal cual:
lines(edades, Lx/sum(Lx), col = "blue")
legend("topright", lty=1, col=c("black","blue"), legend = c("c(x) estable","l(x) escalada"))

cx2 <- exp(-(edades+.5)*-r)*Lx / sum(exp(-(edades+.5)*-r)*Lx)
lines(edades, cx2, col = "red")

# duracion media de una generacion
(TG <- log(R0) / r) # TG = tiempo de generacion

R0

# es un poco diferente que la edad promedia de fx:
(Mbar <- sum((edades + .5) * fx) / sum(fx))
(Abar <- sum((edades + .5) * fx * Lx) / sum(fx * Lx))
# pero muy parecido!

# esto es conveniente porque tal vez no hace falta optimizar r tan exactamente,
# si tenemos una relacion entre r, R0, y TG
log(R0) / TG
log(R0) / Mbar # diferente en cuarto digito, no mal!
log(R0) / Abar
#####################################################################
# ahora a repitir para poblaciones tanatologicas

# hace falta una matriz de dx parecida a la de f(y|x),
# pero no escalada:

# Todo depende del dx!
dx <- Tabla$dx
N  <- length(dx)
# importamos la funcion de fya del lunes
da2fya <- function(da){
	
	N <- length(da)
	# definimos una matriz vacia:
	fya <- matrix(0,N,N)
	# ahora iteramos para rellenarlo con da:
	for (i in 1:N){
		# busca la fila i
		# despues ponemos el vector desde la columna 1 hasta su final
		fya[i,1:(N-i+1)] <- da[i:N]
	}
	# lo unico que hace falta es que las filas suman a 1
	# nota que rowSums(fya) = l(a) = colSums(fya) en este momento
	fya <- fya / rowSums(fya)
	# el ultimo objecto esta devuelto por defecto
	fya
}

# el fya que usamos
fya <- da2fya(dx)

# se usa para reproportionar los eventos y los estocs clasificadas
# por ead cronologica en edad tanatologica
Ny  <- colSums(Nx * fya) # filas crono | columna tano
Ey  <- colSums(Ex * fya)

# calculamos las tasa de fecundidad tanatologicas
fy  <- Ny / Ey

# el ISF tanatatologica ''siempre'' es mas alta que el ISF normal
sum(fy)
# miramos
plot(edades, fy, type = 'l')

# ahora ponemos dx en una matriz no escalada
dxy <- matrix(0,N,N)
for (i in 1:N){
	# filas crono | columnas tano
	dxy[i, 1:(N-i+1)] <- dx[i:N]
}

# decimos que las filas son cronologicas y las columnas son tanatologicas
sum(colSums(dxy * exp(-x*r)) * fy)
exp(-x*r)*Lx
x <- .5:110.5
plot(edades, colSums(dxy * exp(-x*-r)) / sum(dxy * exp(-x*-r)), 
		type = 'l', ylim = c(0,.025))
lines(edades, exp(-x*-r)*Lx / sum(exp(-x*-r)*Lx), 
		col = "blue")

# ahora queremos optimizar el 'r' de la equacion de renovacion
# tanatologica:

ThanoMin <- function(rt,dxy,fy,x=.5:110.5){
	abs(1 - sum(colSums(dxy * exp(-x*rt)) * fy))
}

rt <- optimize(ThanoMin, c(-.05,.05), dxy= dxy, fy=fy, tol = 1e-15)$minimum

r ; rt
sum(colSums(dxy * exp(-x*.011)) * fy)

library(devtools)
install_github("timriffe/Pyramid/Pyramid")

library(Pyramid)

args(Pyramid)

Exh <- Exp[Exp$Year == año, "Male"]

Pyramid(Exh, Ex, grid = FALSE,
		fill.males="royalblue",fill.females="royalblue",
		coh.axis = TRUE,year = 1900, prop = FALSE)
segments(0,0,0,110,col = 'white')
args(Pyramid)

Dxh <- Def[Def$Year == año,"Male"]

Dxh[Dxh == 0] <- .1
Exh[Exh == 0] <- .1

Tablah <- LT(Exh, Dxh, mxsmooth = FALSE)

fyxh <- da2fya(Tablah$dx)
fyxf <- da2fya(Tabla$dx)

Pxyf <- Ex * fyxf
Pxyh <- Exh * fyxh

Pyramid(Pxyh, Pxyf, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = FALSE, 
		border.males = NA, border.females = NA
		)
segments(0,0,0,110,col = 'white')

library(RColorBrewer)
display.brewer.all()
cols <- brewer.pal(9, "YlGn")
colRampa <- colorRampPalette(cols, space = "Lab")
colRampa(111)

Pyramid(Pxyh, Pxyf, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = FALSE, 
		border.males = NA, border.females = NA, 
		fill.males = colRampa(111), fill.females = colRampa(111)
)

aggVec <- function(x,Int){
	ages <- 1:length(x) - 1
	groups <- ages - ages %% Int
	tapply(x,groups,sum, na.rm = TRUE)
}
aggVec(Ex,10)

PxyhAgg <- t(apply(Pxyh,1,aggVec,Int = 10))
PxyfAgg <- t(apply(Pxyf,1,aggVec,Int = 10))

Pyramid(PxyhAgg, PxyfAgg, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = TRUE, 
		border.males = NA, border.females = NA, 
		fill.males = colRampa(11), fill.females = colRampa(11)
)

PyramidOutline(Exh,Ex,scale = 100)




PxyhAggT <- t(apply(Pxyh,2,aggVec,Int = 10))
PxyfAggT <- t(apply(Pxyf,2,aggVec,Int = 10))

Pyramid(PxyhAggT, PxyfAggT, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = TRUE, 
		border.males = NA, border.females = NA, 
		fill.males = rev(colRampa(11)), fill.females = rev(colRampa(11))
)

PyramidOutline(rowSums(PxyhAggT),rowSums(PxyfAggT),scale = 100)

head(Nac)
Nac <- Nac[!is.na(Nac$Births), ]

Nac$Cohort <- Nac$Year - Nac$Age

library(reshape2)

NacAP <- acast(Nac[Nac$Age >= 13 & Nac$Age <= 55, ], 
		       Age~Year, value.var = "Births")
NacAC <- acast(Nac[Nac$Age >= 13 & Nac$Age <= 55, ], 
		Age~Cohort, value.var = "Births")
nrow(NacAP)
barplot(NacAP, border = NA, col = colRampa(43), space = 0)

NacAPAgg <- apply(NacAP,2,aggVec,Int = 5)
NacACAgg <- apply(NacAC,2,aggVec,Int = 5)
nrow(NacAPAgg)
dev.new()
barplot(NacAPAgg, border = NA, col = colRampa(9), space = 0)
barplot(NacACAgg, border = NA, col = colRampa(9), space = 0)




