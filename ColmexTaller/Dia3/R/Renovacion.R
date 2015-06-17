
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

################################################################
# aqui hemos hecho una pausa, y yo he empezado tal cual de programar
# sin explicarme bien. Notas añadida después
################################################################
# ahora queremos optimizar el 'r' de la equacion de renovacion
# tanatologica:
# esta funcion sigue la pauta de la de antes para Lotka,
# adentro se ve una version discreta 
ThanoMin <- function(rt,dxy,fy,x=.5:110.5){
	abs(1 - sum(colSums(dxy * exp(-x*rt)) * fy)) # valor absoluto pk estamos minimizando
}
# esto funciona como lo de antes, que explique en la pizarra
rt <- optimize(ThanoMin, c(-.05,.05), dxy= dxy, fy=fy, tol = 1e-15)$minimum

# comparamos el 'r' crono con el 'r' tano
r ; rt
# debe de ser igual a 1
sum(colSums(dxy * exp(-.5:110.5*rt)) * fy)


# aqui he preguntado por alguno deseo, y se ha pedido piramides
# instalamos otro package desde github entonces para eso:
library(devtools)
install_github("timriffe/Pyramid/Pyramid")
# cargar package
library(Pyramid)
# tiene demasiadas argumentos
args(Pyramid)

# claro, hasta ahora solo mujers, pero aqui necesitamos los dos sexos
Exh <- Exp[Exp$Year == año, "Male"]

# esto hemos ido modificando iterativamente para ver como cambiar cosas
Pyramid(Exh, Ex, grid = FALSE,
		fill.males="royalblue",fill.females="royalblue",
		coh.axis = TRUE,year = 1900, prop = FALSE)
# aun asi queda un poco feita
segments(0,0,0,110,col = 'white') 
args(Pyramid)

# mus guai es descomponer el piramide por años restantes,
# luego dibujarlo de forma mas impactante.
# hace falta mortalidad de hombres 'h'

# defunciones
Dxh <- Def[Def$Year == año,"Male"]
# limpieza
Dxh[Dxh == 0] <- .1
Exh[Exh == 0] <- .1
# calculamos la tabla de vida para los hombres
Tablah <- LT(Exh, Dxh, mxsmooth = FALSE)

# calculamos el dx reescalada para poder descomponer = f(y|x)
fyxh <- da2fya(Tablah$dx) # hombres
fyxf <- da2fya(Tabla$dx)  # mujeres

# Ex = poblaciones mujers | Exh = pob hombres
# Pxy son matrices de poblaciones por años vividos y restantes
# sea conciente de que esto es muy aproximado y tiene unos 
# supositos fuertes, pero lo hacemos porque aun asi da informacion
Pxyf <- Ex * fyxf
Pxyh <- Exh * fyxh

# ahora: tiene demasiados colores, feito
Pyramid(Pxyh, Pxyf, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = FALSE, 
		border.males = NA, border.females = NA
		)
segments(0,0,0,110,col = 'white')

# aqui los colores buenos (sabios) de Cynthia Brewer
library(RColorBrewer)
display.brewer.all() # elijimos paleta
cols <- brewer.pal(9, "YlGn") # desde amarillos hasta verdes
colRampa <- colorRampPalette(cols, space = "Lab") # construimos rampa de interpolacion entre estos colors
colRampa(111) # tenemos 111 categorias de años restantes, nos hace 111 colores...

# mejor que antes, pero aun no legible, hay que agrupar edades
Pyramid(Pxyh, Pxyf, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = FALSE, 
		border.males = NA, border.females = NA, 
		fill.males = colRampa(111), fill.females = colRampa(111)
)

# una funcion que agrupa las edades:
# x: el variable que agrupamos
# Int = el intervalos de edad de salida
aggVec <- function(x,Int){
	# suponemos que las edades son sencillas y empiezan en 0
	ages <- 1:length(x) - 1
	# determinar en cual grupo va cada edad
	groups <- ages - ages %% Int
	# una manera tipica de agrupar (aggregar) cosas en R
	tapply(x,groups,sum, na.rm = TRUE)
}
aggVec(Ex,10) # por ejemplo

# un bucle interno sobre las filas.
# esto quiere decir que operamos sobre cada edad cronologico por separado.
# cada edad cronologica tiene un vector de edades tanatologicas que queremos
# agrupar. 
# Pxy = la matriz ; 1 = filas ; aggVec = la funcion que le aplicamos ; 
# Int = 10 = un argumento que pasamos a la funcion
PxyhAgg <- t(apply(Pxyh,1,aggVec,Int = 10)) # luego el transpuesto por conveniencia
PxyfAgg <- t(apply(Pxyf,1,aggVec,Int = 10))

# mucho mejor!
Pyramid(PxyhAgg, PxyfAgg, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = TRUE, 
		border.males = NA, border.females = NA, 
		fill.males = colRampa(11), fill.females = colRampa(11)
)
# una linea de perfil!
PyramidOutline(Exh,Ex,scale = 100)
# le falta leyenda y quidar los detalles. Otra vez: un grafico artisanal.
# tampoco ha sido mucho trabajo. Pero si quieres hacerlo mas a menudo, mas
# vale simplificar el trabajo tras definir una funcion que los hace todo, o
# preferiblemente en pasos modulares. Es decir, 1) una funcion come dx, y Px
# y te devuelve la matriz descompuesta. de heco, ya tengo uno aqui:
# https://gist.github.com/timriffe/4948828 usalo si quieres!

# 2) una funcion de agrupacion que tambien te hace la parte de apply() seria util
# 3) finalmente un modulo que te dibuja la piramida con leyenda y todo para no tener
# que reinventar la rueda cada vez


###############################################
# ahora invertimos lo de arriba: la misma matriz
# tambien nos permite hacer la piramide por tiempo
# restante descompuesto por la edad.

# solo hace falta agrupar, esta vez sobre columnas en vez de filas (el '2')
PxyhAggT <- t(apply(Pxyh,2,aggVec,Int = 10))
PxyfAggT <- t(apply(Pxyf,2,aggVec,Int = 10))

# y ya nos sale asi de facil:
Pyramid(PxyhAggT, PxyfAggT, grid = FALSE,
		coh.axis = TRUE,year = 1900, prop = TRUE, 
		border.males = NA, border.females = NA, 
		fill.males = rev(colRampa(11)), fill.females = rev(colRampa(11))
)
PyramidOutline(rowSums(PxyhAggT),rowSums(PxyfAggT),scale = 100)

###########################################
# Ojo: estos dos ultimas piramides son 'multiestados' y 
# y funciona programaticamente igual representar datos de
# poblaciones en categorias distintas: estado civil, nivel educativo,
# etc etc. filas = edades, columnas = estados.
###########################################

#############################################
# aqui he preguntado por algun otro deseo, y nadie dijo nada
# asi que he decidio jugar con los nacimientos
#############################################
Nac        <- Nac[!is.na(Nac$Births), ]

# definimos cohorte de nacimiento, pero toma en cuenta
# que es super aproximado, porque los datos son AP, y cada
# AP sencillo tiene 2 cohortes adentros (los 2 triangulos Lexis).
Nac$Cohort <- Nac$Year - Nac$Age

# quiero transformarlo en matrices AP y AC
library(reshape2)

# recomiendo aprender esta funcion acast(). Lo uso casi cada vez que uso R!
NacAP <- acast(Nac[Nac$Age >= 13 & Nac$Age <= 55, ], 
		       Age~Year, value.var = "Births") # edad por periodo
NacAC <- acast(Nac[Nac$Age >= 13 & Nac$Age <= 55, ], 
		Age~Cohort, value.var = "Births") # edad por cohorte
nrow(NacAP)
# otra vez demasiados colores
barplot(NacAP, border = NA, col = colRampa(43), space = 0)
# agrupamos como antes (no hace falta transpuesto aqui)
NacAPAgg <- apply(NacAP,2,aggVec,Int = 5)
NacACAgg <- apply(NacAC,2,aggVec,Int = 5)
nrow(NacAPAgg)

# miramos los dos graficos uno encima del otro. Intenta reescalar la
# ventana grafica. hmmm. En RStudio por defecto te mete los graficos
# en una ventanita, pero busca como sacar el grafico en una ventana libre
# para poder poner los dos al lado. Busca los booms y echos, y intenta
# relacionar la experiencia de los cohortes con lo de los periodos
barplot(NacAPAgg, border = NA, col = colRampa(9), space = 0)
dev.new()
barplot(NacACAgg, border = NA, col = colRampa(9), space = 0)


# fin!

