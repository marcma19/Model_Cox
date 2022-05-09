# | --------------------------------------------------------------
# | Auteurs : Gary Bloch, Ornella Fettaya, Mathias Marciano
# | Description : Code utilisé pour mémoire sur le Modèle de Cox
# | Commentaire : nous conseillons l'utilisation de Rstudio cloud
# | pour compiler le code sinon il sera difficile d'importer la
# | librairie timereg
# | --------------------------------------------------------------

#---------------------Librairies utilisées-------------------
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
library(corrplot)
library(timereg)
set.seed(80)

#---------------------Codes pour generer les figures-------------------

#Courbe de Kaplan-Meier (page 11)
km_estimate = survfit(Surv(time, status) ~ 1, data = lung, type="kaplan-meier")
ggsurvplot(
  fit = km_estimate,
  xlab = "Jours",
  ylab = "Probabilité de survie")

#Courbe de Kaplan-Meier en fonction du sexe (page 11)
ggsurvplot(
  fit = survfit(Surv(time, status) ~ sex, data = lung, type="kaplan-meier"),
  xlab = "Jours",
  ylab = "Probabilité de survie",
  legend.title="Sexe",
  font.main = c(12, "bold", "darkred"),
  legend.labs=c("Homme", "Femme") )

#Courbe de Nelson-Aalen du risque cumule (page 13)
courbe<-summary(survfit(Surv(time,status)~1, data=lung, type="fh"))
risque=cumsum(courbe$n.event/courbe$n.risk) #Calcul du risque d'apres la formule de N.A
plot(courbe$time, risque, type="s", xlab="Jours", ylab = "Risque cumulé",xlim=c(0,800))

#Courbe de Breslow du risque cumule (page 14)
plot(km_estimate$time, -log(km_estimate$surv), xlab="Jours", ylab =
       "Risque cumulé", type='l', xlim=c(0,800))

#Tableau des correlations de lung (page 24)
corrplot(cor(lung,use = "complete.obs"), type="upper", tl.col="black",
         tl.srt=45)

#Variation du coefficient β de wt.loss au cours du temps (page 26)
model<-coxph(Surv(time,status) ~ age+factor(sex)+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, data = lung)
test_proportional<-cox.zph(model)
plot(test_proportional[7],lwd=2,ylim=c(-0.1,0.1))
abline(h= model$coef[7], col=3, lwd=2, lty=2) #Rajout de la droite valeur moyenne de Beta
legend("topleft", legend=c('Variation de Beta au cours
du temps', 'Valeur moyenne de Beta'), lty=c(3,2,1), col=c(1,3,1), lwd=2)

#Variation du coefficient β de ph.karno au cours du temps (page 27)
model<-coxph(Surv(time,status) ~ age+factor(sex)+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, data = lung)
test_proportional<-cox.zph(model)
plot(test_proportional[4],lwd=2,ylim=c(-0.1,0.1))
abline(h= model$coef[4], col=3, lwd=2, lty=2) #Rajout de la droite valeur moyenne de Beta
legend("bottomright", legend=c('Variation de Beta au cours
du temps', 'Valeur moyenne de Beta'), lty=c(3,2,1), col=c(1,3,1), lwd=2)


#---------------------Codes pour generer les modeles-------------------

#Une premiere selection grossiere (page 21)
model_1<-coxph(Surv(time,status) ~ age+factor(sex)+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, data = lung)
cox.zph(model_1)

#Une selection plus poussee (page 22)
model_1<-coxph(Surv(time,status) ~ age+factor(sex)+ph.ecog+pat.karno+wt.loss, data = lung)
cox.zph(model_1)

#Un premier modele de Cox (page 22)
model_1<-coxph(Surv(time,status) ~ age+factor(sex)+ph.ecog+wt.loss, data = lung)
cox.zph(model_1)

#Estimation des coefficients Beta du premier modele (page 23)
summary(model_1)$coef

#AIC du premier modele (page 23)
step(model_1, data = lung, direction=c("forward"))$Start

#Selection en fonction des correlations (page 24)
model_2<-coxph(Surv(time,status) ~ age+factor(sex)+pat.karno+meal.cal+wt.loss, data = lung)
cox.zph(model_2)

#Deuxieme modele de Cox (page 25)
model_2<-coxph(Surv(time,status) ~ age+factor(sex)+pat.karno+wt.loss, data = lung)
cox.zph(model_2)

#Estimation des coefficients Beta du second modele (page 25)
summary(model_2)$coef

#AIC du second modele (page 25)
step(model_2, data = lung, direction=c("forward"))$Start

#Modele avec meal.cal qui depend du temps (page 28)
#La fonction const considere la variable en argument comme constante
set.seed(80)
model_time_meal<-timecox(Surv(time, status) ~ const(age) +const(ph.karno) +meal.cal+ const(pat.karno)+ const(ph.ecog) + const(factor(sex)) + const(wt.loss),data=lung,n.sim=500, max.time=700)
summary(model_time_meal)

#Modele avec ph.karno qui depend du temps (page 29)
#Bien executer le code avec la graine aleatoire
set.seed(80)
model_time_karno<-timecox(Surv(time, status) ~ const(age) + ph.karno +const(meal.cal) +const(pat.karno)+ const(ph.ecog) + const(factor(sex)) + const(wt.loss),data=lung,n.sim=500, max.time=700)
summary(model_time_karno)

#Troisieme modele de cox avec une dependance au temps (page 30)
model_finale<-coxph(formula = Surv(time, status) ~ age + tt(ph.karno) + factor(sex) + wt.loss, data = lung, tt = function(x, t, ...) x * log(t + 20))
summary(model_finale)$coef

#AIC du troisieme modele (page 30)
step(model_finale, direction=c("forward"))$Start

#Modele avec ajout d'une colonne aleatoire (page 30)
random<-runif(228, 0, 1000) #228 realisations entre 0 et 1000 de facon uniforme
data_new<-cbind(lung,random) #ajout de la colonne aleatoire a lung
model<-coxph(Surv(time,status) ~ age+factor(sex)+pat.karno+wt.loss+random, data = data_new)
summary(model)

#AIC du modele avec colonne aleatoire (page 31)
step(model, direction=c("forward"))$Start

#Interpretation modele selectionne (page 32)
model_select<-coxph(Surv(time,status) ~ age+factor(sex)+pat.karno+wt.loss, data = lung)
summary(model_select)
