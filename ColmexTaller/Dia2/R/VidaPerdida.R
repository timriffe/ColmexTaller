
############################################################

# Datos = 2010 Hombres, Mexico, gracias Victor!
Dx <- c(17219, 1178, 653, 482, 407, 371, 353, 347, 350, 364, 394, 444, 
		521, 630, 773, 946, 1140, 1344, 1544, 1732, 1899, 2043, 2163, 
		2259, 2334, 2391, 2434, 2465, 2486, 2501, 2515, 2531, 2553, 2582, 
		2617, 2653, 2690, 2724, 2756, 2786, 2815, 2846, 2882, 2925, 2978, 
		3042, 3117, 3203, 3297, 3400, 3509, 3622, 3735, 3849, 3961, 4074, 
		4187, 4301, 4418, 4537, 4658, 4780, 4902, 5027, 5154, 5287, 5424, 
		5567, 5711, 5853, 5991, 6120, 6236, 6326, 6373, 6372, 6353, 6346, 
		6311, 6214, 6083, 5943, 5785, 5603, 5405, 5188, 4942, 4667, 4367, 
		4046, 3712, 3375, 3039, 2710, 2390, 2080, 1780, 1493, 1222, 972, 
		750, 558, 399, 274, 179, 111, 65, 35, 18, 8)
Ex <- c(1155746, 1146185, 1140057, 1141653, 1148602, 1154777, 1157785, 
		1156630, 1150848, 1142343, 1133299, 1124918, 1119408, 1118061, 
		1119274, 1119380, 1114658, 1103012, 1084576, 1061014, 1034319, 
		1006327, 978592, 952250, 928124, 906693, 887792, 870999, 856095, 
		843161, 832438, 824167, 818345, 814362, 810808, 806085, 798981, 
		788633, 774578, 757033, 736719, 714435, 691045, 667449, 644309, 
		621974, 600474, 579646, 559346, 539424, 519655, 499809, 479753, 
		459473, 439131, 418969, 399173, 379924, 361398, 343600, 326452, 
		309880, 293951, 278736, 264286, 250618, 237744, 225570, 213933, 
		202702, 191807, 181163, 170669, 160090, 149156, 137940, 127213, 
		117544, 108162, 98556, 89305, 80771, 72794, 65300, 58349, 51832, 
		45638, 39798, 34340, 29313, 24750, 20681, 17100, 13987, 11306, 
		9006, 7049, 5395, 4027, 2920, 2051, 1388, 903, 560, 331, 187, 
		99, 49, 22, 9)
Mx <- Dx / Ex
#############################################################################
# hoy vamos a definir unas funciones sencillas. demasiada sencillas de hecho
#############################################################################

#plot(0:109, Dx/Ex, type = 'l', log = 'y')
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
	# lx es el divisor final para calcular ex
	lx1 <- mx2lx(mx)
	# usamos Lx para calcular Tx sobre la marcha
	Lx  <- mx2Lx(mx)
	# la parte de la izquierda es Tx, luego condicionar con lx
	(rev(cumsum(rev(Lx))) / lx1)
}

##############################################
# primero PYLL:
# un poco limpieza...
Mx <- Dx / Ex

edades <- 0:109
plot(edades, Mx, type = 'l', log = 'y', 
		main = "suavizado fuertamente?\n *no importa, pero mas vale saber los supositos\n Parece parametrico casi")
########################################################
#
# el APP estandar se define como Dx * ex
APPx <- Dx * mx2ex(Mx)
barplot(APPx, main = "mucha vida perdida tras la mortalidad infantil!")
# 14392 PYLL0 tras 182 defuniones de infantes!

# ahora, las edades perdidos pero estos infantes se extienden sobre todas las edades superiores..
# por lo tanto, la barra grande de la edad 0 es facil de malinterpretar

# queremos lx condicionada por si mismo, asi:
# 1) definimos una matriz vacia con las dimensiones corectas

hazmeLX <- function(mx){
	N  <- length(mx)
	# caja vacia para asignar Lx en bucle abajo
	LX <- matrix(0, nrow = N, ncol = N)
	# vector de lx, para reescalar Tx al final
	lx <- mx2lx(mx)
	# vector de Lx suponiendo lineal entre edades
	Lx <- mx2Lx(mx)
    # ahora lo rellenamos iterativamente:
	for (i in 1:N){
		LX[i,i:N] <- Lx[i:N]
	}
	#######
	 ######
	  #####
	   ####
	    ###
		 ##
		  #
	# los valors en cada fila suman a Tx
    # ahora condicionado por supervivencia, lx
	LX <- LX / lx
	# ahora los valors de cada fila suman a ex
	LX
}
# compare sumas: uau son iguales!
LX <- hazmeLX(Mx)
sum(Dx * LX) ; sum(APPx)

barplot(colSums(Dx * LX), main = "misma cantidad, pero redistribuida", border = NA, space = 0)

# lo mismo pero desagregado por x de Dx ( variantes de gris)
barplot(Dx * LX, main = "misma cantidad, pero redistribuida", border = NA, space = 0)

PPX <- Dx * LX

###############################################3
# escrito en clase: queremos un grafico con mejor definicion:
# agrupar edades. Hay muchas maneras, y esta es la que se me
# ha ocurrido en clase, asi que anoto aqui los pasos
# nuestras edades sencillas que queremos agrupar: estan en filas
Sencillas <- row(PPX) - 1
# el factor de agrupacion: es decir el intervalo de edad que queremos
Fac <- 10
# %% es el modulo, 
# como cuando tenias 7 años y querias dividir 6/10 = 1 con 4 de resto
Grup <- Sencillas - Sencillas %% Fac
# un vector que usaremos en el bucle para seleccionar...
unicas <- unique(Grup[,1])
# iniciamos ppxout como objeto vacio
 ppxout <- NULL
 # una iteracion por cada grupo de edad
for (i in 1:11){
	# una matriz hecha solo de los del grupo de edad
	PPXi <- matrix(PPX[Grup == unicas[i]],ncol = ncol(PPX))
	# aggregar cada columna: acabamos un un vector
	ppxi <- colSums(PPXi)
	# pegar fila nueva
	ppxout <- rbind(ppxout, ppxi)
}

barplot(ppxout, main = "misma cantidad, pero redistribuida", border = NA, space = 0)

## Ahora, si queremos unos colores mas impactantes, pero tambien
## sin violar mejores practicas con el uso de color en areas, lo
## mas facil siempre es referir al package RColorBrewer. Tambien
## vale para mapas! mira http://colorbrewer2.org/ y otros sitios asi

## 1) installa si hace falta
##install.packages("RColorBrewer")
library(RColorBrewer)

# mira las paletas que tiene:
# display.brewer.all()
# queremos uno del primer grupo en este caso!
BrewCol <- brewer.pal(9,"GnBu") # primer argumento = cuantos colores, segundo cual paleta
# Problema: tenemos 11 clases de edad y solo 9 colores
# solucion: interpolamos colores.
# colorRampPalette() devuelve una funcion para tu paleta. Si le das mas colores la
# interpolacion puede salir mejor. Luego elige uno del los 'espacios' de color.
# RGB es el defecto... y es defectuoso como espacio de color....
# Yo suelo usar Lab, pero existen otros.
MiRampa <- colorRampPalette(BrewCol, space = "Lab")
miscols <- MiRampa(11) # se usa de esta forma
barplot(ppxout, main = "misma cantidad, pero redistribuida", 
		border = NA, space = 0, col = miscols)
# o con los colores al reves:
barplot(ppxout, main = "misma cantidad, pero redistribuida", 
		border = NA, space = 0, col = rev(miscols))
##################################################################
# ahora los colores ligeros se pueden perder contra el blanco.
# quiero ponerle una linea delgada de perfil para darle forma.
# cambio de estrategia: barplot() nos ha hecho lo que puede,
# pero cuando se quiere controlar bien los detalles hay que 
# usar funciones mas primitivos a veces....
##################################################################
# Todo lo de abajo es lo que llamo artisanal. Implica trabajo!
# Solo vale la pena hacerlo para presentaciones, entradas de blog
# articulos, y tal vez tesis. Este tipo de trabajo primitivo NO
# es recomendado si estas explorando dattos... para explorar, usa
# base, lattice o ggplot2. Construimos de esta forma para que se ve
# lo que hace falta para tener control total sobre una figura.
##################################################################
# 1) abrir plot() vacio pero con las dimensiones correctas:
plot(NULL, type = "n", xlim = c(0,110), ylim = c(0, 170000), 
		axes = FALSE, xlab = "", ylab = "")
# 2) usare poligonos para las barras de color. No un poligon
#    por barra vertical, como en barplot(), sino un poligono unico
#    para cada region de color. Es mas complicado, pero sale mejor
#    calidad siempre si lo guardas... 

# estrategia: hacemos la suma cumulativa de cada columna de ppxout:
# apply es un bucle interno | 
# '2' dice que operamos sobre columnas | 
# cumsum es la funcion que applicamos
ppxcum <- apply(ppxout,2,cumsum)
# 3) ahora, cada valor y en esa matriz tiene 2 valores x...
# (porque estamos trabajando con valores discretos, y se usa
# mejor barras en este caso...)
x      <- c(0,rep(1:109,each=2),110)
# 4) añadir fila de ceros en pppxcum, porque el fondo es plano:
ppxcum <- rbind(0,ppxcum)
# 5) listos para un bucle para dibjuar los poligonos!
for ( i in 1:11){
	# el poligono tiene un fondo y un techo! dibuja siempre
	# en la misma direccion, por eso repitimos los x asi.
	polygon(x = c(x,rev(x)), y = rep(c(ppxcum[i,],rev(ppxcum[i+1,])),each=2),
			border = NA, col = rev(miscols)[i])
}
# 6) ya esta pintando mejor!
# ahora ponemos lineas de referencias que no nos molesten!
abline(v = seq(10,100,by=10), col = "#FFFFFF50", lwd = .5) # blaco, semitransparente :-)
abline(h = seq(25000,150000,by = 25000),, col = "#FFFFFF50", lwd = .5)
# ojo: siempre hay que andar con quidado con las lineas de referencia.
# si los hace y piensas que distraen, quitelos! siempre!

# 7) la linea de perfil
lines(x,rep(ppxcum[nrow(ppxcum),],each=2))
# 8) ahora ejes
# segments es super basico.
# 8.1) eje y:
segments(0,0,0,150000)            
segments(0,seq(0,150000,by = 25000), # xpd nos da permiso de dibujar fuera de la area
		-1,seq(0,150000,by = 25000), xpd = TRUE)
options(scipen=6) # para que no nos escribe 100000 en vez de 1e5 ...
text(0,seq(0,150000,by = 25000),        # pos = 2 dice que alinea a la izquierda.
		labels=seq(0,150000,by = 25000),pos=2,xpd=TRUE)
# 8.2) eje x:
segments(seq(0,110,by=10),0,seq(0,110,by=10),-2000,xpd=TRUE)
text(seq(0,110,by=10),-1000,seq(0,110,by=10),pos=1,xpd=TRUE) # pos = 1 = abajo

# Todavia hace falta etiquetas de eje, y una leyenda explicativa con color,
# tendrias que luego quitar el espacio blanco de los margenes. El titulo
# lo pones siempre en el paper y mejor no en el mismo grafico.

# finalmente: gaurdalo como pdf, asi se queda en formato vector. Para Word or PPT,
# se puede investigar como guardar en formato EMF (pero lastima, graficos en formato
# vector de Windows no apoyan la transperencia, y lo hemos usado). NUNCA usa jpg
# para este tipo de grafico. Lo dire mañana en clase.

##########################################################
# paramos aqui: sugiero intentarlo con causes el viernes
# y entonces tomar mas tiempo para visualizarlo bien
##########################################################



