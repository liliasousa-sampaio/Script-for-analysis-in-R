#Scripts to reproduce the multi-model, PCA, and NMDS analyses performed in the article
“Human Disturbances And Woody Flora Diversity Along A Caatinga–Brejo De Altitude Gradient In Northeastern Brazil”
Sousa et.,2026


library(MASS)       
library(rJava)      
library(nnet)       
library(car)        
library(qcc)        
library(glmulti)    
library(AICcmodavg) 
library(aod)        
library(pscl)       
library(leaps)      
library(splines)   
library(survival)   
library(lattice)    
library(Matrix)     
library(lme4)       

setwd(choose.dir())
dir()


##############MULTI-MODELS ###########################
#import the table containing the response and predictor variables "Alpha diversity, Taxonomic beta diversity, or Phylogenetic beta diversity"

Diversidade_alfa.csv<-read.table(file.choose(),row.names=1,header=T,sep=";") 
View(Diversidade_alfa.csv) #ver a tabela e se está OK
str(Diversidade_alfa.csv)


#MUTMODEL##
glmulti(Q_0~elevacao+diversidade+agropecuaria+cad, data=Diversidade_alfa.csv, intercept = TRUE, marginality = FALSE, method = "h", level = 1, crit = "aicc", confsetsize = 1024, family=gaussian)->result1
result1 
summary(result1) 
weightable(result1) 
resultadonti=weightable(result1)
write.csv(resultadonti, "Q_0.csv")


# Here, we obtain the set of possible models, their respective AICc values ​​(Corrected Akaike Information Criterion),
# and the Akaike weights (model weights – wi). These values ​​are useful for identifying:
# - The set of plausible models (models with ΔAICc < 2).
# - The set of models that best approximate the true model (sum of weights wi ≥ 0.95).


weightable(result1) 

# We extracted the coefficient (slope) values ​​to evaluate the effect and direction
# of the explanatory variables in each model.

coef(result1, select="all", newdata=NA, se.fit=FALSE, varweighting="Johnson", icmethod="Burnham", alphaIC=0.05)
coef.glmulti(result1)
write.csv(coef.glmulti(result1), "Estimat_Q_0.csv")



# This step provides the estimators (coefficients) and the unconditional variance.
# The sign of the estimator indicates the effect of the independent variable on the response variable.
# If the unconditional variance is greater than the estimator value, caution is required in interpretation,
# as there may be significant variation—both upward and downward—relative to this value.



mod1 <- glm(Q_0 ~ 1, data=Diversidade_alfa.csv, family = gaussian)
mod2 <- glm(Q_0~elevacao+diversidade+agropecuaria+cad, data=Diversidade_alfa.csv, family = gaussian)
anova(mod1, mod2, test = "Chi")  
write.csv(anova(mod1, mod2, test = "Chi"), "Anova_Q_0.csv")  




#PCA Analysis

library(stats)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggfortify)
library(ade4)
library(vegan)
library(usdm)


setwd(choose.dir())
dir()

Pertubacao.csv<-read.table(file.choose(),row.names=1,header=T,sep=";") #Spreadsheet with the environmental and spatial filter variables 
View(Pertubacao.csv) 
str(Pertubacao.csv)
vif(Pertubacao.csv)


#ANALISE 1 PCA
PCAa <- prcomp(Pertubacao.csv, scale=TRUE) 
summary(PCAa)


# Contributions of variables to eixo 1
p <- fviz_contrib(PCAa, choice = "var", axes = 1, fill = "darkgray",color = "black", top = 17)
p + theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 14, color = "black"))

# Contributions of variables to PC2
p2<- fviz_contrib(PCAa, choice = "var", axes = 2, fill = "darkgray",color = "black", top = 17)
p2 + theme(axis.title = element_text(size = 18),
           axis.text = element_text(size = 14, color = "black"))


# Results variaveis
res.var <- get_pca_var(PCAa)
res.var$coord          

res.var$contrib        

# Results for individuals - 
res.ind <- get_pca_ind(PCAa)

res.ind$coord         
res.ind$contrib        
res.ind$cos2 


autoplot(PCAa)


#edit the chart to highlight the desired elements
#re-import a table containing the elements to be highlighted in the chart

parcelas<- read.table(file.choose(),row.names=1,header=T,sep=",") #Group spreadsheet
attach(parcelas)
View(parcelas)

autoplot(PCAa)
autoplot(PCAa, data = parcelas, colour = 'Group', size = 5) #Group é o nome da coluna da planilha fator 
autoplot(PCAa, data = parcelas, colour = 'Group', label = TRUE, label.size = 3)
autoplot(PCAa, data = parcelas, colour = 'Group', label = TRUE, shape = FALSE, label.size = 5)
autoplot(PCAa, data = parcelas, colour = 'Group', label = TRUE, loadings = TRUE)
autoplot(PCAa, data = parcelas, colour = 'Group', 
         loadings = TRUE, loadings.colour = 'black', size = 5,
         loadings.label = TRUE, loadings.label.colour = "black", loadings.label.size = 4) 
#theme_classic() 



autoplot(PCAa, data = parcelas, colour = 'Group', 
         loadings = TRUE, loadings.colour = 'black', size = 5,
         loadings.label = TRUE, loadings.label.colour = "black", loadings.label.size = 4) +
  scale_color_manual(values = c(
    "CAA"   = "#FF6F61",  # Vermelho suave
    "BREJO" = "#00B050",  # Verde forte
    "TRANS" = "#5DA5E7"   # Azul claro
  )) +
  labs(color = "Ambiente")





#NMDS Analysis: Three tables are imported for this analysis.
Abundancia.csv<-read.table(file.choose(),row.names=1,header=T,sep=";") #Spreadsheet 1 Abundance
View(Abundancia.csv) 
attach(Abundancia.csv)

factor.csv<- read.table(file.choose(),row.names=1,header=T,sep=",") #Spreadsheet 2 Group
View(factor.csv)
attach(factor.csv) 

Planilha_Elev.csv<- read.table(file.choose(),row.names=1,header=T,sep=";") #spreadsheet 3 Environmental and spatial predictor variables
View(Planilha_Elev.csv)
attach(Planilha_Elev.csv)


####metodo BRAY-CURTIS
NMDSB <- metaMDS(Abundancia.csv, distance = "bray", k = 2, trymax=30) 
ordiplot(NMDSB,type="n")
orditorp(NMDSB,display="sites",col="red",air=0.01)
orditorp(NMDSB,display="species",cex=0.5,air=0.01)



MDS1 = NMDSB$points[,1] 
MDS2 = NMDSB$points[,2] 
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Group = Group) #Group is the name of the coloring worksheet.
NMDS

#NMDS 

ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(aes(color = Group, size=4))

#NMDS Grafico 2
ggplot(data = NMDS, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Group), size = 4) +
  scale_color_manual(values = c(
    "CAA"   = "#FF6F61",
    "BREJO" = "#00B050",
    "TRANS" = "#5DA5E7"
  )) +
  labs(color = "Ambiente")  

#Chart templates for Envifit
#Common NMDS plot in R

plot(NMDSB, display = c("sites", "species"), choices = c(1, 2))

plot(NMDSB, display = c("sites"), type="t", choices = c(1, 2))

#analyzing the significance of the NMDS using ANOSIM

floradist<- vegdist(Abundancia.csv, method = "bray") 
simple.results.anosim <- anosim(x = floradist, grouping = Group, permutations = 999)
simple.results.anosim #see R2 and P results

#analyzing the effect of predictor variables on species composition #ENVIFT

fit3 <- envfit(NMDSB,Planilha_Elev.csv, perm=1000)
fit3

plot(fit3, factorp.max = 0.05, col = "blue")





