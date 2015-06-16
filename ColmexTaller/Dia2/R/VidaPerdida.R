
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
# un poco limpieza...
Mx <- Dx / Ex
Mx[is.nan(Mx) | is.infinite(Mx)] <- 0

edades <- 0:109
plot(edades, Mx, type = 'l', log = 'y', 
		main = "suavizado fuertamente?\n *no importa, pero mas vale saber los supositos\n Parece parametrico casi")

# el PYLL estandar se define como Dx * ex
PYLLx <- Dx * mx2ex(Mx)
barplot(PYLLx, main = "mucha vida perdida tras la mortalidad infantil!")
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



