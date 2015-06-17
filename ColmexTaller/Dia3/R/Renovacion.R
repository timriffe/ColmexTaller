
##############################
# vamos a hacer Lotka en R!
##############################

# 1) necesitamos datos mortalidad y fecundidad
library(HMDHFDplus)
# 2) para buenas tablas de vida
library(LifeTable) # si no lo tienes, mira en el codigo de lunes
##############################
# si hace falta, define username & password
#us <- userInput()
#pw <- userInput()
Pais <- "SWE"
Def <- readHMDweb(Pais,"Deaths_1x1",us,pw)
Exp <- readHMDweb(Pais,"Exposures_1x1",us,pw)

Nac <- readHFDweb(Pais,"birthsRR",us,pw)

# los datos de HMD y HFD son de dimensiones distintas

# tienen años distintos
range(Nac$Year)
range(Def$Year)

# y los nacimientos de HFD no estan clasificados por sexo:
head(Nac)
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
Nac$Total <- Nac$Total * PF # claro, solo haz esto una vez
##########################################################################

# herramientas para seleccionar datos:
sexo <- "Female"
año  <- 1900

# buscamos datos de 1900 porque si:
Dx <- Def[Def$Year == año, sexo]
Ex <- Exp[Exp$Year == año, sexo]
# N son nacimientos para nosotros, pero nota que muchas veces en
# la literatura se usa para estocs de poblacion (Px)
Nx <- Nac$Total[Nac$Year == año]
# los Nx van de la edad 12 al 55 solo. 
# añadimos 0s antes y despues para alinear edades
Nx <- c(rep(0,12),Nx,rep(0,55))
# calculamos tasas de fecundidad:
fx <- Nx / Ex
# limpieza, sino tenemos casos de NaN:
fx[Ex == 0] <- 0

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
cx <- exp(-(edades+.5)*r)*Lx / sum(exp(-(edades+.5)*r)*Lx)
plot(edades, cx, type = 'l')
# compare con el Lx tal cual:
lines(edades, Lx/sum(Lx), col = "blue")
legend("topright", lty=1, col=c("black","blue"), legend = c("c(x) estable","l(x) escalada"))

# duracion media de una generacion
(TG <- log(R0) / r) # TG = tiempo de generacion

# es un poco diferente que la edad promedia de fx:
(Mbar <- sum((edades + .5) * fx) / sum(fx))
# pero muy parecido!

# esto es conveniente porque tal vez no hace falta optimizar r tan exactamente,
# si tenemos una relacion entre r, R0, y TG
log(R0) / TG
log(R0) / Mbar # diferente en cuarto digito, no mal!

#####################################################################
# ahora a repitir para poblaciones tanatologicas

