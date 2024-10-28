####################################

#           Data Quesos

#####################################
#
#   ANALISIS DE MEDIDAS REPETIDAS
library(MASS)
library(stats)
library(rstatix)
library(tidyverse)
library(ggplot2)
library(datarium)
library(purrr)
library(plyr)
library(dplyr)
library(ggpubr)
base_quesos2 #Base de datos

#Para peso
quesos2_medR <- base_quesos2[, -(4:10)] # Conservar solo variables necesarias
quesos2_medR

quesos2_long <- quesos2_medR %>% # Formato largo
  gather(key = "time", value = "score", p0,p1,p2) %>%
  convert_as_factor(id,sustrato)
print(quesos2_long)

#        SUPUESTOS
# 1ro. OUTLIERS
quesos2_long %>%
  group_by (sustrato, time) %>%
  identify_outliers (score)
# No hay outliers extremos

# 2do. NORMALIDAD
quesos2_long %>%
  group_by (sustrato, time) %>%
  shapiro_test (score)
ggqqplot (quesos2_long, "score", facet.by = "time")
# P>0.05 Por lo tanto se cumple.

# ANÁLISIS
rep_quesos2v <- anova_test (quesos2_long, dv = score, wid = id,
                          within = c(sustrato, time)) #Analisis de interac/2 vías
get_anova_table (rep_quesos2v) 
# Todas las P<0.05
# Tanto la interaccion como los factores de manera independiente cambian con el t. 

#### POST-HOC
prueba_post2 <- quesos2_long %>%
  pairwise_t_test (score ~ time, paired = T, p.adjust.method = "bonferroni")
prueba_post2
# Todas las P<=.05 por lo tanto, el peso cambia signif. a lo largo del t.

#### INTERACCIONES
#Agrupación por tratamiento para ver su efecto en el t
inter_2 <- quesos2_long %>%
  group_by (sustrato) %>%
  anova_test (dv = score, wid = id, within = time) %>%
  get_anova_table ( ) %>%
  adjust_pvalue (method = "bonferroni")
inter_2
# El peso cambia a lo largo del t con eltratamiento

####    BOXPLOT
library(ggplot2)
quesos2_long #Base

box_quesos <- ggplot(quesos2_long, aes(x=time, y=score, fill = sustrato))+
  geom_boxplot()+
  ggtitle("Analisis de medidas repetidas para Peso (p)")
box_quesos


#------------------------------------------------------------------------------

##############################

#     PCA
library(permute)
library(lattice)
library(vegan)

quesos_Vc <- base_quesos[,-c(1,2,3,4,5,7,8,10,11)] # Vcont
quesos_Vc
str(quesos_Vc)
matx_qVc <- as.matrix(quesos_Vc) # Matriz

anal_qVc <- prcomp(matx_qVc, scale. = TRUE) 
anal_qVc
summary(anal_qVc) # % de varc de mis CP's

anal_qVc$x #Componenetes del analisis

pc1_x <- anal_qVc$x [ , 1]
pc2_y <- anal_qVc$x [ , 2]

##  ANOVA
analisis_pc1 <- aov(pc1_x ~ sustrato, base_quesos) #ANOVA del PC1

summary.aov(analisis_pc1) # Rechazo la Ho

TukeyHSD(analisis_pc1)# Que Ha acepto

library(ggplot2)
boxp_PC1 <- ggplot(analisis_pc1, aes(x = sustrato, y = pc1_x))+ 
  geom_boxplot()+
  ggtitle(" PC1 ~ sustrato (P≤0.001)")
boxp_PC1

#-------------------------------------------------------------------------------

####################################
#
#      Análisis de SIMILITUD

library(permute)
library(lattice)
library(vegan)
base_qSim <- base_quesos[, c(6,9,10,12)]

analisis_simil <- anosim(base_qSim, grouping = base_quesos$sustrato, 
                         distance = "manhattan")
summary(analisis_simil)
# Significance<0.05
# La variacion dentro del grupo es menor a la que existe entre grupos.

## Hacer ELIPSE
library(ggplot2)

plot(pc1_x, pc2_y, xlab = "PC1 (46%)", ylab = "PC2 (34%)",
     xlim = c(-4, 4), ylim = c(-4, 4), 
     abline(h = 0, v = 0, lty = "dotted"))
ordiellipse( anal_qVc , base_quesos$substrate , conf = 0.95, col = 1:4)
text(x=-3.5, y=3.5, label = "Anosim R=0.58; P<0.05", cex = 0.8)

#-------------------------------------------------------------------------------

#################################

#        DISCRIMINANTES
library(MASS)
base_quesos2

base_quesos2.1 <- base_quesos2[,-c(1,3,4,5,7,8,11)]

quesos_discrimts <- lda( sustrato ~ . , data = base_quesos2.1)
quesos_discrimts #Analisis

val_pred_coord <- predict(quesos_discrimts)
val_pred_coord # Extrac de LD1 y LD2

# GRAFICAR $x
crd_x <- val_pred_coord$x [ , 1]   
crd_y <- val_pred_coord$x [ , 2] # Coordenadas 

tipo_q <- base_quesos2.1$sustrato # Obj. con el sustrato
etiq_color <- function(variable)
{
  if( length( grep("oaxaca" , variable)))
    ("#9ACD32")
  else if( length( grep("panela" , variable)))
    ("#8B008B")
  else
    ("#00CED1")
}
etiq_colq <- unlist( lapply( tipo_q , etiq_color) )

#Gráfico
graf_disc <- plot(crd_x , crd_y , xlab = "1er Discriminante (96%)",
                  ylab = "2do Discriminante (3%)", type = "n",
                  main="Análisis de discriminantes")
text(crd_x, crd_y, base_quesos2.1$sustrato, cex = 0.7, col = etiq_colq)

#-------------------------------------------------------------------------------

################################
#     HOTELLING para t_MO
library(corpcor)
library(Hotelling)
library(ggplot2)

base_quesos1.2 <- base_quesos[-c(1:26),] # Ver el efecto del tp de MO
base_quesos1.2

Vq_resp <- cbind(base_quesos1.2$area_o,base_quesos1.2$pH_2,base_quesos1.2$p2) #Comb. de Vr 

hot_met2 <- hotelling.test(Vq_resp~t_MO, base_quesos1.2) # Análisis
hot_met2

#      PRUEBAS DE T
# p2 - sustrato
test1 <- t.test(p2~t_MO, base_quesos1.2)
test1 #No singinf

# pH_2 - Sustrato
test2 <- t.test(pH_2~t_MO, base_quesos1.2) 
test2 # P<0.05
# Unico signif. 
boxp_t2 <- boxplot(pH_2~t_MO, base_quesos1.2,
                   main = "T-test P≤0.001",
                   xlab = "Tipo de microorganismo",
                   ylab = "pH final")

# area_o - Sustrato
test3 <- t.test(area_o~t_MO, base_quesos1.2)
test3 # P=0.05 

#-------------------------------------------------------------------------------

################
#   MANOVA
base_quesos

# No se cumplen los supuestos por lo que no lo considere...

base_quesos$sustrato <- as.factor(base_quesos$sustrato)
base_quesos$p2 <- as.numeric(base_quesos$p2)
str(base_quesos)

quesos_bact <- base_quesos[-c(27,29,30,31,39) , ] #Solo bact
quesos_bact

quesos_mano <- quesos_bact[ ,-c(1,3,4,5,7,8,10,11)] #Modif. de la base  
quesos_mano

Vq_resp <- cbind(quesos_mano$area_o,quesos_mano$pH_2,quesos_mano$p2) #Comb. de Vr 
Vq_resp

matriz_Vq <- as.matrix(Vq_resp) #Matris de Vr

mano_queso <- manova(matriz_Vq~sustrato, quesos_mano) #Análisis

summary.manova(mano_queso)
summary.aov(mano_queso)

###################
# SUPUESTOS
library(rstatix)

# 1er MULTI-NORMALIDAD
mshapiro_test(Vq_resp) #No se cumple.
#NO se puede realizar el manova o mancov

#-------------------------------------------------------------------------------

