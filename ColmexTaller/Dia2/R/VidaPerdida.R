
############################################################
# suponemos que las instalaciones de ayer ahora funcionan!
# si no es asi, busca Tim!
############################################################

# 1) carga algun package:

library(HMDHFDplus)

# otra vez a meter contraseña, etc. Definelo como variable
# para no dejar tu info guardado en el fichero de R...
us <- userInput() # esto es interactivo
pw <- userInput()

Pais <- "SWE"
Def <- readHMDweb(Pais,"Deaths_1x1",username = us, password = pw)
Exp <- readHMDweb(Pais,"Exposures_1x1",username = us, password = pw)
# Def <- local(get(load("Dia1/Datos/Def.Rdata")))
# Exp <- local(get(load("Dia1/Datos/Exp.Rdata")))

#############################################################################
# hoy vamos a definir unas funciones sencillas. demasiada sencillas de hecho
#############################################################################

# 1) convertir Mx a lx rapidamente
# fingiendo que estamos en tiempo continuo
mx2lx <- function(mx){
	c(1,exp(-cumsum(mx)))[1:length(mx)]
}
# se puede comparar el error implicado en otro momento
# 2) covertir mx en Lx rapidamente:
mx2Lx <- function(mx){
	lx1 <- mx2lx(mx)
	# convertir en Lx, suponiendo linear entre edades (mal para infantes)
	(lx1 + c(lx1[-1],0)) / 2 # Lx
}
# 3) convertir mx a ex rapidamente:
mx2ex <- function(mx){
	Lx  <- mx2Lx(mx)
	# la parte de la izquierda es Tx, luego condicionar con lx
	(rev(cumsum(rev(Lx))) / lx1)[1:length(mx)]
}

##############################################
# primero PYLL:

sexo <- "Male"
año  <- 2000

# sabemos que va de la edad 0 a 110+
Dx <- Def[Def$Year == año, sexo]
Ex <- Exp[Exp$Year == año, sexo]

# un poco limpieza...
Mx <- Dx / Ex
Mx[is.nan(Mx) | is.infinite(Mx)] <- 0

# el PYLL estandar se define como Dx * ex
PYLLx <- Dx * mx2ex(Mx)
barplot(PYLLx, main = "mucha vida perdido tras la mortalidad infantil!")
# 14392 PYLL0 tras 182 defuniones de infantes!

# ahora, las edades perdidos pero estos infantes se extienden sobre todas las edades superiores..
# por lo tanto, la barra grande de la edad 0 es facil de malinterpretar

# queremos lx condicionada por si mismo, asi:
# 1) definimos una matriz vacia con las dimensiones corectas
N  <- length(Dx)
LX <- matrix(0, nrow = N, ncol = N)
lx <- mx2lx(Mx)
Lx <- mx2Lx(Mx)
# ahora lo rellenamos iterativamente:
for (i in 1:N){
	LX[i,i:N] <- Lx[i:N]
}
# ahora condicionado por supervivencia, lx
LX <- LX / lx

# compare sumas: uau son iguales!
sum(Dx * LX) ; sum(PYLLx)

barplot(colSums(Dx * LX), main = "misma cantidad, pero redistribuida", border = NA, space = 0)

# lo mismo pero desagregado por x de Dx ( variantes de gris)
barplot(Dx * LX, main = "misma cantidad, pero redistribuida", border = NA, space = 0)

##########################################################
# paramos aqui: sugiero intentarlo con causes el viernes
# y entonces tomar mas tiempo para visualizarlo bien
##########################################################



