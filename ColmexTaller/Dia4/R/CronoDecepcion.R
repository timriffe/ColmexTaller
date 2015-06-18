
# Author: tim
###############################################################################
# Quiero demostrar como una curva de edad de alguna caracteristica puede ser  #
# deceptiva por fines de proyeccion.                                          #
###############################################################################

###############################################################################
# PD: este codigo viene du un proyecto mio en curso, pero lo pongo aqui       #
#     solo para hacer un punto en la clase de jueves.                         #
###############################################################################

###############################################################################
# 1) definimos algunas funciones:

# el modelo Gompertz de mortalidad: nos devuelve m(x)
mxgomp <- function(a,b,x){
	a*exp(b*x)
}
# metodo barato de covertir m(x) en l(x)
mx2lx <- function(mx){
	N <- length(mx)
	mx <- c(mx[-c(N-1,N)],sum(mx[(N-1):N]))
	c(1, exp(-cumsum(mx)))
}
# convertir l(x) en d(x)
lx2dx <- function(lx){
	-diff(c(lx,0))
}
# convertir l(x) en e(0) [lineal]
lx2e0 <- function(lx){
	sum((lx + c(lx[-1],0)) / 2)
}

# el antilogit = expit
expit <- function(g){
	exp(g) / (1+exp(g))
}

######################################################
# ahora, definimos dos regimenes de mortalidad. Uno alto, el otro mas bajo.

# x = la edad
x   <- 0:100
# hacemos m(x)
mx1 <- mxgomp(a=.0005,b=.075,x)
mx2 <- mxgomp(a=.0005,b=.065,x) # el azul es mas bajo

# mirar las dos m(x)
#plot(mx1,log='y', type = 'l')
#lines(mx2,col = "blue")

# hacemos l(x) desde m(x)
lx1 <- mx2lx(mx1)
lx2 <- mx2lx(mx2)
# cuales son los e(0) implicados?
lx2e0(lx1);lx2e0(lx2) # diferencia de unos 7 años

# miramos l(x) nada raro aqui
par(mai=c(.9,.9,.1,.1))
plot(x,lx1,type='l', xlab = "chrono age", ylab = "l(x)")
lines(x,lx2,col = "blue")
legend("topright",col=c("black","blue"),lty=1,legend=c("l(x) A","l(x) B"),bty="n")

# hacemos d(x) para los dos casos
dx1 <- lx2dx(lx1)
dx2 <- lx2dx(lx2)

############################################################
# gran paso: imaginamos una pauta de alguna caracteristica
# que varie solo por edad cronologica y otra funcion de edad
# que es tanatologica.
# chrono and thano patterns:

# g(x) es alguna caracteristica cronologica (dolor de espalda?)
gx <- expit(seq(-6,4,length.out=101))
# t(y) es alguna caracteristica tanatologica (algun ADL?)
ty <- expit(seq(1,-80,length.out = 101))

# miramos g(x)
par(mai=c(.9,.9,.1,.1))
plot(x,gx,type='l', xlab = "edad cronologica", ylab = "proporcion con condicion g")

# como cambia la prevalencia en la poblacion estacionaria
# si solo cambiamos el regimen de mortalidad. especificamente,
# si bajamos la mortalidad?
par(mai=c(.9,.9,.1,.1))
plot(x,lx1,type='l', xlab = "chrono age", ylab = "l(x) and l(x)C(x)")
lines(x,gx,col="red")
polygon(c(x,rev(x)),c(rep(0,101),rev(gx*lx1)),col="#FF000050")
lines(x,lx2,col = "blue")
polygon(c(x,rev(x)),c(rep(0,101),rev(gx*lx2)),col="#0000FF50")
legend(73,.8,col=c("black","blue","red"),
		lty=c(1,1,1),
		legend=c("l(x) A","l(x) B","g(x)"),bty="n")
legend(73,.65,col=c("#FF000050","#0000FF50"),
		fill=c("#FF000050","#0000FF50"),
		legend=c("G(x) con mort 1","G(x) con mort 2"),bty="n")
# la respuesta: sube la prevalencia de condicion 'g' si es que la condicion es
# cronologica.

##########################################
# unos indices para el caso cronologico: #
##########################################
# Crono: aumenta de esperanza de vida aumenta la prevalencia de la condicion g
sum(gx * lx1) # esperanza de morbididad 1
sum(gx * lx2) # esperanza de morbididad 2 (mas grande en caso de variacion cronologica)
lx2e0(lx1) - sum(gx*lx1) # e0 saludabre 1
lx2e0(lx2) - sum(gx*lx2) # e0 saludabre 2 (ha subido, pero la morbididad ha subido mas!)

############################################################################
# ahora miramos que pasa si es que la caracteristica es funcion de la edad
# tanatologica SOLO.
############################################################################

# miramos t(y): solo sube al final de la vida: recuerda que el tiempo mueve a la izquierda!
par(mai=c(.9,.9,.1,.1))
plot(x, ty, type = 'l', xlab = "edad tanatologica",ylab =  "proporcion con condicion T")

##################################3
# para calcular, mejor tener t(y) en una matriz:
Thx <- matrix(0,101,101)
# columnas = crono ; filas = duracion de vida. despues ponderamos con d(x)
for (i in 1:101){
	Thx[i,1:(102-i)] <- ty[(102-i):1]
}
# dibuja como es la pauta en la matriz para que se vea

######################################
# como cambia la prevalencia de T bajo los dos regimenes de mortalidad?
par(mai=c(.9,.9,.1,.1))
plot(x,lx1,type='l', xlab = "edad cronologica", ylab = "l(x), T(x), T(x)/l(x)")
polygon(c(x,rev(x)),c(rep(0,101),rowSums(Thx*rev(dx1))),col="#FF000050")
lines(x,lx2,col="blue")
polygon(c(x,rev(x)),c(rep(0,101),rowSums(Thx*rev(dx2))),col="#0000FF50")
lines(x,colSums(Thx*rev(dx1))/lx1, col = "black",lty=2)
lines(x, colSums(Thx * rev(dx2))/lx2, col = "blue",lty=2)

legend(75,1,col=c("black","blue","black","blue"),
		lty=c(1,1,2,2),
		legend=c("l(x) mort 1","l(x) mort 2","t'(x) mort 1","t'(x) mort 2"),bty="n")
legend(75,.83,col=c("#FF000050","#0000FF50"),
		fill=c("#FF000050","#0000FF50"),
		legend=c("T(x) mort 1","T(x) mort 2"),bty="n")
# la prevalencia baja cuando la mortalidad baja!
####################################################
par(mai=c(.9,.9,.1,.1))
plot(x,colSums(Thx*rev(dx1))/lx1, type = 'l',xlab = "edad cronologica", ylab = "proporcion con condicion T")
lines(x, colSums(Thx * rev(dx2))/lx2, col = "blue")
legend("topleft",1,col=c("black","blue",NA,"#FF000050","#0000FF50"),
		lty=c(1,1),
		legend=c("t'(x) mort 1","t'(x) mort 2"),bty="n")

####################################################
# indices de resumen:
# cambio trivial en la esperanza de morbididad
sum(Thx*rev(dx1)) # esperanza de morbididad 1
sum(Thx*rev(dx2)) # esperanza de morbididad 2
lx2e0(lx1) - sum(Thx*rev(dx1)) # esperanza de salud 1
lx2e0(lx2) - sum(Thx*rev(dx2)) # esperanza de salud 2 # caso 2 es una mejora muy fuerte!

100 * (1 - sum(Thx*rev(dx1)) / lx2e0(lx1)) # esperanza % saludabre
100 * (1 - sum(Thx*rev(dx2)) / lx2e0(lx2)) # algo mejor ---
                                           # tampoco puedes esperar mucho mas aqui, estamos acercando a 100...

