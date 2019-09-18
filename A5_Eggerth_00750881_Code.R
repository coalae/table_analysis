# AUFGABENBLATT 5
# Cordula Eggerth (00750881)

# Verwendete Literaturquellen:
# .	Folien und R-Codes zu den bisher vorgetragenen Kapiteln aus UK Erweiterungen des linearen Modells (Prof. Wilfried Grossmann, 2019).
# .	Chi-Square-Test Assumptions (2019): http://www.simafore.com/blog/bid/56480/2-key-assumptions-to-be-aware-of-before-applying-the-chi-square-test.
# .	Chi-Quadrat-Test (2019): https://us.sagepub.com/sites/default/files/upm-binaries/820 20_Cha pter_11.pdf. 
# .	John McDonald (2014): Small Numbers in Chi-Square and G-Tests, Handbook of Biological Statistics. http://www.biostathandbook.com/small.html.  


rm(list=ls())

#install.packages("lme4")
library(lme4)
#install.packages("faraway")
library(faraway)


#***********************************************************************************************
# AUFGABE 1 (Zweidimensionale Tabellen)
#***********************************************************************************************
# Ein Hersteller von Büroartikel fertigt Aktenordner in den Farben gelb, rot und 
# blau an. Zur Analyse des Kaufverhaltens werden für die drei wichtigsten Absatzmärkte 
# A, B und C die Anzahlen der in einem bestimmten Zeitintervall georderten 
# Bestelleinheiten ermittelt. 
# Die Werte sind folgender Tabelle zu entnehmen: siehe Angabe.

# a) Beantworte mit einem loglinearen Modell die Frage ob in den drei Märkten die 
#    Ordnerfarben unterschiedlich beliebt sind.
# b) Stelle die Daten mit einem Mosaicplot dar.
# c) Stelle den Zusammenhang mittels Korrespondenzanalyse dar und interpretiere 
#    den Plot.

# tabelle erstellen
y_ordner <- c(564,672,611, 309,198,307, 448,299,425)
markt <- gl(3,3,labels=c("A","B","C"))
farbe <- gl(3,1,9,labels=c("gelb","rot","blau"))
df_ordner <- data.frame(y_ordner, markt, farbe)
tabelle_ordner <- xtabs(y_ordner ~ markt + farbe)

# DESKRIPTIVE STATISTIKEN
summary(df_ordner$y_ordner)

# randhaeufigkeiten (absolut, relativ)
# absolut
addmargins(tabelle_ordner)
# relativ
addmargins(round(prop.table(tabelle_ordner)*100,
                 digits=1))

# LOG-LINEARES MODELL
# (annahme: basierend auf poissonverteilungsmodell)
# (haupteffekte)
mod.ordner <- glm(y_ordner ~ markt + farbe, 
                  family=poisson, df_ordner)
summary(mod.ordner)

drop1(mod.ordner, test="Chi")

# p-value & chi-quadrat krit. wert
1 - pchisq(mod.ordner$deviance,2)
qchisq(0.95,2)

# saturiertes modell
mod.ordner.sat <- glm(y_ordner ~ markt * farbe, 
                      family=poisson, df_ordner)
summary(mod.ordner.sat)

# chi-quadrat-test auf die ordneranzahlen
tab.mod.ordner <- matrix(y_ordner,nrow=3,byrow=TRUE)
chisq.test(tab.mod.ordner) # ergebnis: chi-qu.-test signif.

# BARCHART
colors <- c("darkseagreen1","darkseagreen3","darkseagreen4")
percentage <- prop.table(tabelle_ordner)*100

barplot(percentage, main="Betrachtung in %", ylab="Ordnerfarbe", 
        col=colors, horiz=TRUE)
legend("topright", legend=rownames(percentage), 
       title="Markt", col=colors, lwd=3, lty=3)

# MOSAICPLOT
mosaicplot(tabelle_ordner, color=colors, 
           main="Mosaicplot zu Ordnerbestelleinheiten",
           xlab="Markt", ylab="Ordnerfarbe")

# KORRESPONDENZANALYSE
# residuenmatr. 
residual_matr <- xtabs(residuals(mod.ordner, type="pearson") ~
                         farbe + markt, df_ordner)
residual_matr

# SVD (singular value decomposition)
# fuer erste 2 komponenten
svd_res <- svd(residual_matr,2,2)
svd_res

# zerlege in u- und v-komponenten
sv1 <- svd_res$u %*% diag(sqrt(svd_res$d[1:2]))
sv2 <- svd_res$v %*% diag(sqrt(svd_res$d[1:2]))

# inertia 
inertia <- svd_res$d[1]^2+svd_res$d[2]^2
inertia

# vgl. der inertia mit chi-quadrat
# (vorauss.: erwartetete haeuf. groesser 5) >> hier erfuellt
res <- matrix(numeric(ncol(tabelle_ordner)*nrow(tabelle_ordner)),
              ncol=ncol(tabelle_ordner))
for(i in 1:nrow(tabelle_ordner)){
  for(j in 1:ncol(tabelle_ordner)){
    res[i,j] <- rowSums(tabelle_ordner)[i] * 
      colSums(tabelle_ordner)[j]/sum(tabelle_ordner)
  }
}
all(res>5)

chisq.test(tabelle_ordner)
# ergebnis: chi-quadr.-test signifikant 
# >> zshg. zwischen markt und ordnerfarbe

gesamt <- t(svd_res$d) %*% svd_res$d
gesamt # ergibt selbiges wie chisq teststat. bzw. inertia
       # anzeichen fuer korrekte skalierung 

# korrespondenzanalyse plot
aa <-1.1 * max(abs(sv1),abs(sv2)) 
plot(rbind(sv1,sv2), asp=1, 
     xlim=c(-aa,aa),ylim=c(-aa,aa), 
     xlab="SV1",ylab="SV2",type="n") 
abline(h=0,v=0) 
text(sv2,c("A","B","C"), cex=0.7, col="mediumslateblue") 
text(sv1,c("gelb","rot","blau"), cex=0.7, col="lightsalmon3")
legend("topright", cex=0.8, legend=c("sv2","sv1"), 
       col=c("mediumslateblue","lightsalmon3"), 
       lwd=2, lty=1)

tt <- rbind(sv2,sv1)
tt


#*********************************************************************************************
# AUFGABE 2 (Zweidimensionale Tabellen)
#***********************************************************************************************
# Die folgende Tabelle gibt die Beziehungen zwischen den Leistungen von Studenten in
# Mathematik und Statistik an.
# Tabelle: siehe Angabe.

# a) Beantworte mit einem loglinearen Modell die Frage ob es einen Zusammenhang 
#    zwischen den Noten in Mathematik und Statistik gibt.
# b) Stelle die Daten mit einem Mosaicplot dar und mit geeigneten Prozentwerte in
#    einem Barplot dar.
# c) Stelle den Zusammenhang mittels Korrespondenzanalyse dar und interpretiere 
#    den Plot.

# tabelle erstellen
y_leistungen <- c(56,71,12, 37,163,38, 24,42,85)
labels_leistungen <- c("Sehr gut", "Durchschnitt", "Schlecht")
statistik <- gl(3,3,labels=labels_leistungen)
mathematik <- gl(3,1,9,labels=labels_leistungen)
df_leistungen <- data.frame(y_leistungen, mathematik, statistik)
tabelle_leistungen <- xtabs(y_leistungen ~ statistik + mathematik)
 
# DESKRIPTIVE STATISTIKEN 
summary(df_leistungen$y_leistungen)

# randhaeufigkeiten (absolut, relativ)
  # absolut
addmargins(tabelle_leistungen)
  # relativ
addmargins(round(prop.table(tabelle_leistungen)*100,
                   digits=1))

# LOG-LINEARES MODELL
# (annahme: basierend auf poissonverteilungsmodell)
# (haupteffekte)
mod.leistungen <- glm(y_leistungen ~ statistik + mathematik, 
                      family=poisson, df_leistungen)
summary(mod.leistungen)

drop1(mod.leistungen, test="Chi")

  # p-value
1 - pchisq(mod.leistungen$deviance,2)

  # krit. wert
qchisq(0.95,2)

# ergebnis: zshg. zwischen leistungen in math. und stat., 
#           weil p value signifikant 

# saturiertes modell
mod.leistungen.sat <- glm(y_leistungen ~ statistik * mathematik, 
                      family=poisson, df_leistungen)
summary(mod.leistungen.sat)

# analyse mit hilfe von chi-quadrat-test
tab.mod <- matrix(y_leistungen,nrow=3,byrow=TRUE)
chisq.test(tab.mod)

# ergebnis: chi-qu.-test ist signifikant


# BARCHART
colors <- c("darkseagreen1","darkseagreen3","darkseagreen4")
percentage <- prop.table(tabelle_leistungen)*100

barplot(percentage, main="Betrachtung in %", ylab="Mathematik", 
        col=colors, horiz=TRUE)
legend("topright", legend=rownames(percentage), 
       title="Statistik", col=colors, lwd=3, lty=3)

# MOSAICPLOT
mosaicplot(tabelle_leistungen, color=colors, 
           main="Mosaicplot zu Leistungen der Studenten",
           xlab="Statistik", ylab="Mathematik")

# ergebnis: wenn laenge und breite des plots in vierecken
#           gleich lang, dann sind die anteile ca. gleich


# KORRESPONDENZANALYSE
# residuenmatr. 
residual_matr <- xtabs(residuals(mod.leistungen, type="pearson") ~
                 mathematik + statistik, df_leistungen)
residual_matr

# SVD (singular value decomposition)
# fuer erste 2 komponenten
svd_res <- svd(residual_matr,2,2)
svd_res

# zerlege in u- und v-komponenten
sv1 <- svd_res$u %*% diag(sqrt(svd_res$d[1:2]))
sv2 <- svd_res$v %*% diag(sqrt(svd_res$d[1:2]))

# inertia 
inertia <- svd_res$d[1]^2+svd_res$d[2]^2
inertia

# vgl. der inertia mit chi-quadrat
# (vorauss.: erwartetete haeuf. groesser 5) >> hier erfuellt
res <- matrix(numeric(ncol(tabelle_leistungen)*nrow(tabelle_leistungen)),
                      ncol=ncol(tabelle_leistungen))
for(i in 1:nrow(tabelle_leistungen)){
  for(j in 1:ncol(tabelle_leistungen)){
    res[i,j] <- rowSums(tabelle_leistungen)[i] * 
                colSums(tabelle_leistungen)[j]/sum(tabelle_leistungen)
  }
}
all(res>5)

chisq.test(tabelle_leistungen)
# ergebnis: chi-quadr.-test signifikant 
# >> zshg. zwischen math.- und stat.-leistungen

gesamt <- t(svd_res$d) %*% svd_res$d
gesamt # ergibt selbiges wie chisq teststat. bzw. inertia
       # d.h. korrekte skalierung 

# korrespondenzanalyse plot
aa <-1.1 * max(abs(sv1),abs(sv2)) 
plot(rbind(sv1,sv2), asp=1, 
     xlim=c(-aa,aa),ylim=c(-aa,aa), 
     xlab="SV1",ylab="SV2",type="n") 
abline(h=0,v=0) 
text(sv2,labels_leistungen, cex=0.7, col="mediumslateblue") 
text(sv1,labels_leistungen, cex=0.7, col="lightsalmon3")
legend("topright", cex=0.8, legend=c("sv2","sv1"), 
       col=c("mediumslateblue","lightsalmon3"), 
       lwd=2, lty=1)

tt <- rbind(sv2,sv1)
tt


#***********************************************************************************************
# AUFGABE 3 (Zweidimensionale Tabellen)
#***********************************************************************************************
# In einem Experiment sollen zwei Testmethoden verglichen werden. Jeder Test soll dabei
# von je 10 Personen durchgeführt werden. Die Ergebnisse des Tests sind dabei nur bestanden 
# und nicht bestanden. Man untersuche die folgenden vier Testergebnisse:
# Tabellen: siehe Angabe.

# Bei welchen Ergebnissen liefert der Chi-Quadrat-Test signifikante Ergebnisse?
# Wie ändern sich die Ergebnisse, wenn man anstelle von 10 Personen je Test jeweils 20 
# oder 40 Personen die Tests durchführen lässt, also die Werte in den Tabellen mit dem 
# Faktor 2 bzw. 4 multipliziert? Begründe die Änderungen.

 
# daten >> tabellen erstellen
  # factor levels
tests <- gl(2,2,labels=c("Test1","Test2"))
bestehen <- gl(2,1,4,labels=c("bestanden","nicht bestanden"))

  # ERGEBNIS 1
y_ergebnisse1 <- c(2,8, 8,2)
df_ergebnisse1 <- data.frame(y_ergebnisse1, tests, bestehen)
tabelle_ergebnisse1 <- xtabs(y_ergebnisse1 ~ tests + bestehen)

 # ERGEBNIS 2
y_ergebnisse2 <- c(2,8, 3,7)
df_ergebnisse2 <- data.frame(y_ergebnisse2, tests, bestehen)
tabelle_ergebnisse2 <- xtabs(y_ergebnisse2 ~ tests + bestehen)

 # ERGEBNIS 3
y_ergebnisse3 <- c(2,8, 4,6)
df_ergebnisse3 <- data.frame(y_ergebnisse3, tests, bestehen)
tabelle_ergebnisse3 <- xtabs(y_ergebnisse3 ~ tests + bestehen)

 # ERGEBNIS 4
y_ergebnisse4 <- c(3,7, 7,3)
df_ergebnisse4 <- data.frame(y_ergebnisse4, tests, bestehen)
tabelle_ergebnisse4 <- xtabs(y_ergebnisse4 ~ tests + bestehen)

# chi-quadrat-tests fuer 10-personen-pro-test-setting
tab.ergebnisse1 <- matrix(y_ergebnisse1,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse1)

tab.ergebnisse2 <- matrix(y_ergebnisse2,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse2)

tab.ergebnisse3 <- matrix(y_ergebnisse3,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse3)

tab.ergebnisse4 <- matrix(y_ergebnisse4,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse4)
  # ergebnisinterpretation:
  # - keiner der tests ist signifikant auf 0.05 level
  # - ausserdem ist die voraussetzung, dass mind. 5 
  #   beobachtungen pro zelle in jeder tabelle sind,  
  #   nicht erfuellt

# chi-quadrat-tests fuer 20-personen-pro-test-setting
# annahme: gleichbleibende relative aufteilung der gruppen

 # ERGEBNIS 1
y_ergebnisse1.pers20 <- c(2,8, 8,2)*2
df_ergebnisse1.pers20 <- data.frame(y_ergebnisse1.pers20, tests, bestehen)
tabelle_ergebnisse1.pers20 <- xtabs(y_ergebnisse1.pers20 ~ tests + bestehen)

tab.ergebnisse1.pers20 <- matrix(y_ergebnisse1.pers20,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse1.pers20) # signifikant
                                   # vorauss. nicht erfuellt
 
 # ERGEBNIS 2
y_ergebnisse2.pers20 <- c(2,8, 3,7)*2
df_ergebnisse2.pers20 <- data.frame(y_ergebnisse2.pers20, tests, bestehen)
tabelle_ergebnisse2.pers20 <- xtabs(y_ergebnisse2.pers20 ~ tests + bestehen)

tab.ergebnisse2.pers20 <- matrix(y_ergebnisse2.pers20,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse2.pers20) # nicht signifikant
                                   # vorauss. nicht erfuellt

 # ERGEBNIS 3
y_ergebnisse3.pers20 <- c(2,8, 4,6)*2
df_ergebnisse3.pers20 <- data.frame(y_ergebnisse3.pers20, tests, bestehen)
tabelle_ergebnisse3.pers20 <- xtabs(y_ergebnisse3.pers20 ~ tests + bestehen)

tab.ergebnisse3.pers20 <- matrix(y_ergebnisse3.pers20,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse3.pers20) # nicht signifikant
                                   # vorauss. nicht erfuellt

 # ERGEBNIS 4
y_ergebnisse4.pers20 <- c(3,7, 7,3)*2
df_ergebnisse4.pers20 <- data.frame(y_ergebnisse4.pers20, tests, bestehen)
tabelle_ergebnisse4.pers20 <- xtabs(y_ergebnisse4.pers20 ~ tests + bestehen)

tab.ergebnisse4.pers20 <- matrix(y_ergebnisse4.pers20,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse4.pers20) # signifikant
                                   # vorauss. erfuellt


# chi-quadrat-tests fuer 40-personen-pro-test-setting
# annahme: gleichbleibende relative aufteilung der gruppen

# ERGEBNIS 1
y_ergebnisse1.pers40 <- c(2,8, 8,2)*4
df_ergebnisse1.pers40 <- data.frame(y_ergebnisse1.pers40, tests, bestehen)
tabelle_ergebnisse1.pers40 <- xtabs(y_ergebnisse1.pers40 ~ tests + bestehen)

tab.ergebnisse1.pers40 <- matrix(y_ergebnisse1.pers40,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse1.pers40) # signifikant
                                   # vorauss. erfuellt

# ERGEBNIS 2
y_ergebnisse2.pers40 <- c(2,8, 3,7)*4
df_ergebnisse2.pers40 <- data.frame(y_ergebnisse2.pers40, tests, bestehen)
tabelle_ergebnisse2.pers40 <- xtabs(y_ergebnisse2.pers40 ~ tests + bestehen)

tab.ergebnisse2.pers40 <- matrix(y_ergebnisse2.pers40,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse2.pers40) # nicht signifikant
                                   # vorauss. erfuellt

# ERGEBNIS 3
y_ergebnisse3.pers40 <- c(2,8, 4,6)*4
df_ergebnisse3.pers40 <- data.frame(y_ergebnisse3.pers40, tests, bestehen)
tabelle_ergebnisse3.pers40 <- xtabs(y_ergebnisse3.pers40 ~ tests + bestehen)

tab.ergebnisse3.pers40 <- matrix(y_ergebnisse3.pers40,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse3.pers40) # nicht signifikant
                                   # vorauss. erfuellt

# ERGEBNIS 4
y_ergebnisse4.pers40 <- c(3,7, 7,3)*4
df_ergebnisse4.pers40 <- data.frame(y_ergebnisse4.pers40, tests, bestehen)
tabelle_ergebnisse4.pers40 <- xtabs(y_ergebnisse4.pers40 ~ tests + bestehen)

tab.ergebnisse4.pers40 <- matrix(y_ergebnisse4.pers40,nrow=2,byrow=TRUE)
chisq.test(tab.ergebnisse4.pers40) # signifikant
                                   # vorauss. erfuellt


#***********************************************************************************************
# AUFGABE 4 (Dreidimensionale Tabellen)
#***********************************************************************************************
# Die Daten femsmoke in der Library faraway zeigen die Ergebnisse einer Studie 
# über das Rauchen bei Frauen in den Jahren 1972 - 1974. Die Variable y gibt die 
# Anzahl der Fälle in den Gruppen an, die durch Raucher, Tot und Altersgruppe 
# gebildet werden. 
# Die kleinen Fallzahlen in manchen Altersgruppen ergeben sich dadurch, dass Personen 
# im Laufe der Untersuchung ausgeschieden werden mussten. Beachte bei der Modellierung, 
# dass dieser Datensatz ein Beispiel für Simpsons Paradoxon ist. In den einzelnen 
# Altersgruppen sind die Ergebnisse anders als das Gesamtergebnis über alle Altersgruppen.

# DATEN UND DESKRIPTIVES: 
# daten laden und tabelle erstellen
data(femsmoke)
head(femsmoke, n=5)
summary(femsmoke)

is.data.frame(femsmoke)    # ist schon data.frame
is.factor(femsmoke$smoker) # ist schon factor
is.factor(femsmoke$dead)   # ist schon factor
is.factor(femsmoke$age)    # ist schon factor

tab.fem <- xtabs(femsmoke$y ~ femsmoke$smoker + femsmoke$dead + femsmoke$age)
tab.fem
summary(tab.fem)

# MOSAICPLOT:
par(mfrow=c(1,1))
mosaicplot(tab.fem, color=c("dodgerblue","dodgerblue3",
                            "dodgerblue4","darkslateblue",
                            "lightblue","lightcyan4","navajowhite3"),
           cex=0.5)

# LOG-LINEARES MODELL / MODELLE PRUEFEN:
# zshg. dead - smoker 
mod1.fem <- xtabs(y ~ smoker + dead, femsmoke)
summary(mod1.fem) # signifikanter zshg.
round(prop.table(mod1.fem,1), digits=3)

# zshg. age - smoker 
mod2.fem <- xtabs(y ~ age + smoker, femsmoke)
summary(mod2.fem) # signifikanter zshg.
round(prop.table(mod2.fem,1), digits=3)

# zshg. age - dead 
mod3.fem <- xtabs(y ~ age + dead, femsmoke)
summary(mod3.fem) # signifikanter zshg.
round(prop.table(mod3.fem,1), digits=3)

# zshg. age - dead - smoker
mod4.fem <- xtabs(y ~ age + dead + smoker, femsmoke)
summary(mod4.fem) # signifikanter zshg.
round(prop.table(mod4.fem,1),digits=3)

# MODELL "TOTALE UNABHAENGIGKEIT":
m1.fem <- glm(y ~ age + dead + smoker, data = femsmoke, family = poisson) 
summary(m1.fem)
c(deviance(m1.fem), df.residual(m1.fem))
qchisq(0.95,df.residual(m1.fem))

# MODELL DER 2-FACH-INTERAKTIONEN:
m2.fem <- glm(y ~ (age + dead + smoker)^2, data=femsmoke, family=poisson) 
summary(m2.fem)
c(deviance(m2.fem), df.residual(m2.fem))
qchisq(0.95,df.residual(m2.fem))
drop1(m2.fem, test="Chisq")

# MODELL D. BEDINGTEN UNABH. (von smoker und age gegeben dead)  
m3.fem <- glm(y ~ dead*smoker + dead*age, data=femsmoke, family=poisson) 
summary(m3.fem)
c(deviance(m3.fem), df.residual(m3.fem))
qchisq(0.95,df.residual(m3.fem))
drop1(m3.fem, test="Chisq")
par(mfrow=c(2,2))
plot(m3.fem)

# SUCHE NACH MOEGLICHKEITEN EINFACHERER MODELLE
m4.fem <- glm(y ~ age + dead*smoker, data=femsmoke, family=poisson) 
summary(m4.fem)
c(deviance(m4.fem), df.residual(m4.fem))
qchisq(0.95,df.residual(m4.fem))
drop1(m4.fem, test="Chisq")

# MODELLVERGLEICH
anova(m3.fem, m4.fem, test="Chisq")

# ZSFG. UEBER VARIABLE smoker 
  # smoker variable gesamt
tab.smoke <- xtabs(y ~ dead + age, data=femsmoke) 
tab.smoke
chisq.test(tab.smoke)
round(prop.table(tab.smoke), digits=2)

  # smoker=="yes"
tab.smoke.yes <- xtabs(y ~ dead + age, data=femsmoke,
                       subset=(smoker=="yes")) 
tab.smoke.yes
fisher.test(tab.smoke.yes) # voraussetzung nicht erfuellt
chisq.test(tab.smoke.yes) # voraussetzung nicht erfuellt  
summary(tab.smoke.yes)
addmargins(tab.smoke.yes)
round(prop.table(tab.smoke.yes), digits=2)

  # smoker=="no"
tab.smoke.no <- xtabs(y ~ dead + age, data=femsmoke,
                       subset=(smoker=="no")) 
tab.smoke.no
fisher.test(tab.smoke.no) # voraussetzung nicht erfuellt
chisq.test(tab.smoke.no)  
summary(tab.smoke.no)
addmargins(tab.smoke.no)
round(prop.table(tab.smoke.no), digits=2)


# ZSFG. UEBER VARIABLE dead 
  # dead variable gesamt
tab.dead <- xtabs(y ~ smoker + age, data=femsmoke) 
tab.dead
chisq.test(tab.dead)
round(prop.table(tab.dead,1), digits=2)

  # dead=="yes"
tab.dead.yes <- xtabs(y ~ smoker + age, data=femsmoke,
                       subset=(dead=="yes")) 
tab.dead.yes
fisher.test(tab.dead.yes) # voraussetzung nicht erfuellt
chisq.test(tab.dead.yes) # voraussetzung nicht erfuellt  
summary(tab.dead.yes)
addmargins(tab.dead.yes)
round(prop.table(tab.dead.yes,1), digits=2)

  # dead=="no"
tab.dead.no <- xtabs(y ~ smoker + age, data=femsmoke,
                      subset=(dead=="no")) 
tab.dead.no
fisher.test(tab.dead.no) # voraussetzung nicht erfuellt
chisq.test(tab.dead.no)  
summary(tab.dead.no)
addmargins(tab.dead.no)
round(prop.table(tab.dead.no,1), digits=2)


# ZSFG. UEBER VARIABLE age 
# age variable gesamt
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke) 
tab.age
chisq.test(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="18-24"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                     subset=(age=="18-24")) 
tab.age
fisher.test(tab.age) # voraussetzung nicht erfuellt
chisq.test(tab.age)  
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="25-34"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                      subset=(age=="25-34")) 
tab.age
fisher.test(tab.age)
chisq.test(tab.age)   # voraussetzung nicht erfuellt
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="35-44"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                 subset=(age=="35-44")) 
tab.age
fisher.test(tab.age) # voraussetzung nicht erfuellt
chisq.test(tab.age)  
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="45-55"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                 subset=(age=="45-55")) 
tab.age
fisher.test(tab.age) # voraussetzung nicht erfuellt
chisq.test(tab.age)  
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="55-64"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                 subset=(age=="55-64")) 
tab.age
fisher.test(tab.age) # voraussetzung nicht erfuellt
chisq.test(tab.age)  
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="65-74"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                 subset=(age=="65-74")) 
tab.age
fisher.test(tab.age) # voraussetzung nicht erfuellt
chisq.test(tab.age)  
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)

# age=="75+"
tab.age <- xtabs(y ~ smoker + dead, data=femsmoke,
                 subset=(age=="75+")) 
tab.age
fisher.test(tab.age) # voraussetzung nicht erfuellt
chisq.test(tab.age)  
summary(tab.age)
addmargins(tab.age)
round(prop.table(tab.age,1), digits=2)


#***********************************************************************************************
# AUFGABE 5 (Dreidimensionale Tabellen)
#***********************************************************************************************
# Im Datensatz suicide in der Library faraway findet man die Ergebnisse von Selbstmorden 
# in Großbritannien. Die Kreuzklassifizierungsvariablen sind:
#   . Cause = Methode des Selbstmordes
#   . Age = Altersgruppe (y = young, m = middle, o = old)
#   . Sex = Geschlecht (m = male, f = female)
# Lassen sich bestimmte "Vorlieben" für Methoden in bestimmten Kombinationen von Alter 
# und Geschlecht erkennen?
# In diesem Beispiel könnte man auch eine Korrespondenzanalyse für die zweidimensionale 
# Tabelle durchführen, die sich aus Methode und der Kombination von Altersgruppe und 
# Geschlecht ergibt.


# DATEN LADEN UND DESKRIPTIVES:
data(suicide)
head(suicide, n=5)

summary(suicide) 

is.data.frame(suicide)   # ist schon data.frame
is.factor(suicide$cause) # ist schon factor
is.factor(suicide$age)   # ist schon factor
is.factor(suicide$sex)   # ist schon factor

tab.sui <- xtabs(suicide$y ~ suicide$cause + suicide$age + suicide$sex)
tab.sui
summary(tab.sui)

# MOSAICPLOT:
par(mfrow=c(1,1))
mosaicplot(tab.sui, color=c("dodgerblue","dodgerblue3",
                            "dodgerblue4","darkslateblue",
                            "lightblue","lightcyan4","navajowhite3"),
           cex=0.7)

# LOG-LINEARES MODELL / MODELLE PRUEFEN:
# zshg. cause - age
mod1.sui <- xtabs(y ~ cause + age, suicide)
summary(mod1.sui) # signifikanter zshg.
round(prop.table(mod1.sui,1), digits=3)

# zshg. cause - sex 
mod2.sui <- xtabs(y ~ cause + sex, suicide)
summary(mod2.sui) # signifikanter zshg.
round(prop.table(mod2.sui,1), digits=3)

# zshg. age - sex 
mod3.sui <- xtabs(y ~ age + sex, suicide)
summary(mod3.sui) # signifikanter zshg.
round(prop.table(mod3.sui,1), digits=3)

# zshg. age - sex - cause
mod4.sui <- xtabs(y ~ age + sex + cause, suicide)
summary(mod4.sui) # signifikanter zshg.
round(prop.table(mod4.sui,1),digits=3)

# MODELL DER TOTALEN UNABHAENGIGKEIT:
m1.sui <- glm(y ~ age + sex + cause, data = suicide, family = poisson) 
summary(m1.sui)
c(deviance(m1.sui), df.residual(m1.sui))
qchisq(0.95,df.residual(m1.sui))

# MODELL DER 2-FACH-INTERAKTIONEN:
m2.sui <- glm(y ~ (age + sex + cause)^2, data=suicide, family=poisson) 
summary(m2.sui)
c(deviance(m2.sui), df.residual(m2.sui))
qchisq(0.95,df.residual(m2.sui))
drop1(m2.sui, test="Chisq")

# MODELL D. BEDINGTEN UNABH. (von cause und age gegeben sex):  
m3.sui <- glm(y ~ sex*age + sex*cause, data=suicide, family=poisson) 
summary(m3.sui)
c(deviance(m3.sui), df.residual(m3.sui))
qchisq(0.95,df.residual(m3.sui))
drop1(m3.sui, test="Chisq")
par(mfrow=c(2,2))
plot(m3.sui)

# SUCHE NACH EINFACHEREN MODELLEN:
m4.sui <- glm(y ~ age + sex*cause, data=suicide, family=poisson) 
summary(m4.sui)
c(deviance(m4.sui), df.residual(m4.sui))
qchisq(0.95,df.residual(m4.sui))
drop1(m4.sui, test="Chisq")

# MODELLVERGLEICH
anova(m3.sui, m4.sui, test="Chisq")


# ZSFG. UEBER VARIABLE sex 
# sex variable gesamt
tab.sex <- xtabs(y ~ age + cause, data=suicide) 
tab.sex 
chisq.test(tab.sex) # signifikant
round(prop.table(tab.sex), digits=2)

# sex=="m"
tab.sex.m <- xtabs(y ~ age + cause, data=suicide,
                       subset=(sex=="m")) 
tab.sex.m
fisher.test(tab.sex.m) # voraussetzung nicht erfuellt
chisq.test(tab.sex.m)  # signifikant
summary(tab.sex.m)
addmargins(tab.sex.m)
round(prop.table(tab.sex.m), digits=2)

# sex=="f"
tab.sex.f <- xtabs(y ~ age + cause, data=suicide,
                   subset=(sex=="f")) 
tab.sex.f
fisher.test(tab.sex.f) # voraussetzung nicht erfuellt
chisq.test(tab.sex.f)  # signifikant
summary(tab.sex.f)
addmargins(tab.sex.f)
round(prop.table(tab.sex.f), digits=2)


# ZSFG. UEBER VARIABLE age 
# age variable gesamt
tab.age <- xtabs(y ~ sex + cause, data=suicide) 
tab.age 
chisq.test(tab.age) # signifikant
round(prop.table(tab.age), digits=2)

# age=="y"
tab.age.y <- xtabs(y ~ age + cause, data=suicide,
                   subset=(age=="y")) 
tab.age.y
fisher.test(tab.age.y) # voraussetzung nicht erfuellt
chisq.test(tab.age.y)  # signifikant
summary(tab.age.y)
addmargins(tab.age.y)
round(prop.table(tab.age.y), digits=2)
 
# age=="o"
tab.age.o <- xtabs(y ~ age + cause, data=suicide,
                   subset=(age=="o")) 
tab.age.o
fisher.test(tab.age.o) # voraussetzung nicht erfuellt
chisq.test(tab.age.o)  # signifikant
summary(tab.age.o)
addmargins(tab.age.o)
round(prop.table(tab.age.o), digits=2)

# age=="m"
tab.age.m <- xtabs(y ~ age + cause, data=suicide,
                   subset=(age=="m")) 
tab.age.m
fisher.test(tab.age.m) # voraussetzung nicht erfuellt
chisq.test(tab.age.m)  # signifikant
summary(tab.age.m)
addmargins(tab.age.m)
round(prop.table(tab.age.m), digits=2)


# ZSFG. UEBER VARIABLE cause 
# cause variable gesamt
tab.cause <- xtabs(y ~ sex + age, data=suicide) 
tab.cause 
chisq.test(tab.cause) # signifikant
round(prop.table(tab.cause), digits=2)

# cause=="drug"
tab.cause.drug <- xtabs(y ~ sex + age, data=suicide,
                   subset=(cause=="drug")) 
tab.cause.drug
fisher.test(tab.cause.drug)
chisq.test(tab.cause.drug)  # signifikant
summary(tab.cause.drug)
addmargins(tab.cause.drug)
round(prop.table(tab.cause.drug), digits=2)

# cause=="gas"
tab.cause.gas <- xtabs(y ~ sex + age, data=suicide,
                           subset=(cause=="gas")) 
tab.cause.gas
fisher.test(tab.cause.gas)
chisq.test(tab.cause.gas)  # signifikant
summary(tab.cause.gas)
addmargins(tab.cause.gas)
round(prop.table(tab.cause.gas), digits=2)

# cause=="gun"
tab.cause.gun <- xtabs(y ~ sex + age, data=suicide,
                       subset=(cause=="gun")) 
tab.cause.gun
fisher.test(tab.cause.gun)
chisq.test(tab.cause.gun)  # signifikant
summary(tab.cause.gun)
addmargins(tab.cause.gun)
round(prop.table(tab.cause.gun), digits=2)

# cause=="hang"
tab.cause.hang <- xtabs(y ~ sex + age, data=suicide,
                        subset=(cause=="hang")) 
tab.cause.hang
fisher.test(tab.cause.hang)
chisq.test(tab.cause.hang)  # signifikant
summary(tab.cause.hang)
addmargins(tab.cause.hang)
round(prop.table(tab.cause.hang), digits=2)

# cause=="jump"
tab.cause.jump <- xtabs(y ~ sex + age, data=suicide,
                        subset=(cause=="jump")) 
tab.cause.jump
fisher.test(tab.cause.jump)
chisq.test(tab.cause.jump)  # signifikant
summary(tab.cause.jump)
addmargins(tab.cause.jump)
round(prop.table(tab.cause.jump), digits=2)

# cause=="other"
tab.cause.other <- xtabs(y ~ sex + age, data=suicide,
                         subset=(cause=="other")) 
tab.cause.other
fisher.test(tab.cause.other)
chisq.test(tab.cause.other)  # signifikant
summary(tab.cause.other)
addmargins(tab.cause.other)
round(prop.table(tab.cause.other), digits=2)


# ALTERSGRUPPE UND GESCHLECHT ZUSAMMENLEGEN IN 1 FACTOR: 
  # kombinationen geschlecht-altersgruppe:
  # m-m, m-o, m-y, f-m, f-o, f-y
m_m <- suicide[suicide$sex=="m" & suicide$age=="m", ]
m_o <- suicide[suicide$sex=="m" & suicide$age=="o", ]
m_y <- suicide[suicide$sex=="m" & suicide$age=="y", ]
f_m <- suicide[suicide$sex=="f" & suicide$age=="m", ]
f_o <- suicide[suicide$sex=="f" & suicide$age=="o", ]
f_y <- suicide[suicide$sex=="f" & suicide$age=="y", ]

combined.df <- data.frame(y=m_m$y, cause=m_m$cause, sex_age=rep("mm",6))
df_mo <- data.frame(y=m_o$y, cause=m_o$cause, sex_age=rep("mo",6))
df_my <- data.frame(y=m_y$y, cause=m_y$cause, sex_age=rep("my",6))
df_fm <- data.frame(y=f_m$y, cause=f_m$cause, sex_age=rep("fm",6))
df_fo <- data.frame(y=f_o$y, cause=f_o$cause, sex_age=rep("fo",6))
df_fy <- data.frame(y=f_y$y, cause=f_y$cause, sex_age=rep("fy",6))
combined.df <- rbind(combined.df, df_mo, df_my, df_fm, df_fo, df_fy)

head(combined.df, n=10)


# KORRESPONDENZANALYSE
# basis
mod.sui <- glm(y ~ cause + sex_age,
               family=poisson, combined.df)
summary(mod.sui)

tabelle_sui <- xtabs(combined.df$y ~ combined.df$cause + combined.df$sex_age)
tabelle_sui

# residuenmatrix
residual_matr <- xtabs(residuals(mod.sui, type="pearson") ~
                         cause + sex_age, combined.df)
residual_matr

# SVD (singular value decomposition)
# fuer erste 2 komponenten
svd_res <- svd(residual_matr,2,2)
svd_res

# zerlege in u- und v-komponenten
sv1 <- svd_res$u %*% diag(sqrt(svd_res$d[1:2]))
sv2 <- svd_res$v %*% diag(sqrt(svd_res$d[1:2]))

# inertia 
inertia <- svd_res$d[1]^2+svd_res$d[2]^2
inertia

# vgl. der inertia mit chi-quadrat
# (vorauss.: erwartetete haeuf. groesser 5) >> hier erfuellt
res <- matrix(numeric(ncol(tabelle_sui)*nrow(tabelle_sui)),
              ncol=ncol(tabelle_sui))
for(i in 1:nrow(tabelle_sui)){
  for(j in 1:ncol(tabelle_sui)){
    res[i,j] <- rowSums(tabelle_sui)[i] * 
      colSums(tabelle_sui)[j]/sum(tabelle_sui)
  }
}
all(res>5)

chisq.test(tabelle_sui)
# ergebnis: chi-quadr.-test signifikant 
# >> zshg. zwischen cause und sex_age

gesamt <- t(svd_res$d) %*% svd_res$d
gesamt # ergibt ungefaehr selbiges wie chisq teststat. 
       # bzw. inertia
# anzeichen fuer korrekte skalierung 

# korrespondenzanalyse plot
aa <-1.1 * max(abs(sv1),abs(sv2)) 
par(mfrow=c(1,1))
plot(rbind(sv1,sv2), asp=1, 
     xlim=c(-aa,aa),ylim=c(-aa,aa), 
     xlab="SV1",ylab="SV2",type="n") 
abline(h=0,v=0) 
text(sv2,c("mm","mo","my","fm","fo","fy"), cex=0.9, col="mediumslateblue") 
text(sv1,c("drug","gas","gun","hang","jump","other"), cex=0.9, col="lightsalmon3")
legend("topright", cex=0.8, legend=c("sv2","sv1"), 
       col=c("mediumslateblue","lightsalmon3"), 
       lwd=2, lty=1)

tt <- rbind(sv2,sv1)
tt


#***********************************************************************************************
# AUFGABE 6 (Dreidimensionale Tabellen)
#***********************************************************************************************
# Die folgenden beiden Tabellen stellen die Überlebenden und Toten beim Untergang 
# der Titanic gegliedert nach Geschlecht, und Passagierklasse dar.
# Tabelle: siehe Angabe. 

# Man untersuche den Zusammenhang mit einem loglinearen Modell und vergleiche 
# die Ergebnisse mit einer Analyse mit einer logistischen Regression.

# daten anlegen und tabelle erstellen
y_titanic <- c(119,155,422,670, 62,25,88,192,
               5,14,106,3, 141,93,90,20)
ueberleben<-gl(2,4,16,labels=c("nein","ja")) 
geschlecht<-gl(2,8,16,labels=c("maennlich","weiblich")) 
passagierklasse<-gl(4,1,16,labels=c("klasse1","klasse2",
                                   "klasse3","mannschaft")) 
titanic<-data.frame(y,ueberleben,geschlecht,passagierklasse) 
titanic

tab.titanic <-xtabs(y_titanic ~ ueberleben + geschlecht + 
                                passagierklasse) 
tab.titanic

# MOSAICPLOT
par(mfrow=c(1,1))
mosaicplot(tab.titanic, color=c("dodgerblue","dodgerblue3",
                                "dodgerblue4","darkslateblue"),
           cex=0.7)

# LOG-LINEARES MODELL

# zshg. ueberleben - geschlecht 
mod1.titanic <- xtabs(y ~ ueberleben + geschlecht, titanic)
summary(mod1.titanic) # signifikant

prop.table(mod1.titanic,1)

# zshg. passagierklasse - geschlecht 
mod2.titanic <- xtabs(y ~ passagierklasse + geschlecht, titanic)
summary(mod2.titanic) # signifikant

prop.table(mod2.titanic,1)

# zshg. ueberleben - passagierklasse 
mod3.titanic <- xtabs(y ~ ueberleben + passagierklasse, titanic)
summary(mod3.titanic) # signifikant

prop.table(mod3.titanic,1)

# zshg. ueberleben - geschlecht - passagierklasse
mod4.titanic <- xtabs(y ~ ueberleben + geschlecht + 
                        passagierklasse, titanic)
summary(mod4.titanic) # signifikant

prop.table(mod4.titanic,1)

# analyse in den passagierklassen 
mod.pk1 <- xtabs(y ~ ueberleben + geschlecht, titanic,
                 subset=(passagierklasse=="klasse1")) 
mod.pk1
summary(mod.pk1) # signifikant


mod.pk2 <- xtabs(y ~ ueberleben + geschlecht, titanic,
                 subset=(passagierklasse=="klasse2")) 
mod.pk2
summary(mod.pk2) # signifikant


mod.pk3 <- xtabs(y ~ ueberleben + geschlecht, titanic,
                 subset=(passagierklasse=="klasse3")) 
mod.pk3
summary(mod.pk3) # signifikant


mod.pk4 <- xtabs(y ~ ueberleben + geschlecht, titanic,
                 subset=(passagierklasse=="klasse1")) 
mod.pk4
summary(mod.pk4) # signifikant


# modell der totalen unabhaengigkeit
mod.totaleunabh <- glm(y ~ ueberleben + geschlecht + 
                         passagierklasse, family=poisson)

summary(mod.totaleunabh)
c(deviance(mod.totaleunabh), df.residual(mod.totaleunabh))

qchisq(0.95, df.residual(mod.totaleunabh))

# modell der 2fach-interaktionen
mod.2fach <- glm(y ~ (ueberleben + geschlecht + 
                        passagierklasse)^2, family=poisson) 
summary(mod.2fach)

qchisq(0.95, df.residual(mod.2fach))

# modell vereinfachen
drop1(mod.2fach, test="Chi")

# modell der bedingten unabhaengigkeit
# von ueberleben und passagierklasse 
# geg. geschlecht
mod.bedingt <- glm(y ~ ueberleben*geschlecht + 
                       geschlecht*passagierklasse, 
                  family=poisson)
summary(mod.bedingt)

qchisq(0.95,df.residual(mod.bedingt))

par(mfrow=c(2,2))
plot(mod.bedingt)

# sonstige vereinfachung beurteilen
drop1(mod.bedingt, test="Chi")

mod.bedingt2 <- glm(y ~ ueberleben + 
                      geschlecht*passagierklasse, 
                    family=poisson) 
summary(mod.bedingt2)

qchisq(0.95, df.residual(mod.bedingt2))

# vgl. modelle der bedingten unabhaengigkeit
anova(mod.bedingt,mod.bedingt2)

# zsfg. hinsichtlich geschlecht

  # zsfg. ueberleben - passagierklasse
(zsfg_ue_pa<-xtabs(y~ueberleben+passagierklasse, titanic))
prop.table(zsfg_ue_pa,1)
summary(zsfg_ue_pa) # signifikant

# fuer geschlecht==maennlich
(zsfg_ue_pa_m <-xtabs(y~ueberleben+passagierklasse, titanic,
                      subset=(geschlecht=="maennlich")))
prop.table(zsfg_ue_pa_m,1)
summary(zsfg_ue_pa_m) # signifikant

# fuer geschlecht==weiblich
(zsfg_ue_pa_w <-xtabs(y~ueberleben+passagierklasse, titanic,
                      subset=(geschlecht=="weiblich")))
prop.table(zsfg_ue_pa_w,1)
summary(zsfg_ue_pa_w) # signifikant


# LOGISTISCHE REGRESSION
res.logit.full <- glm(ueberleben ~ geschlecht*passagierklasse, 
                 family=binomial(link=logit), data=titanic)
res.logit.full
summary(res.logit.full)

res.logit <- glm(ueberleben ~ geschlecht + passagierklasse, 
                 family=binomial(link=logit), data=titanic)
res.logit
summary(res.logit)

res.logit2 <- glm(ueberleben ~ passagierklasse, 
                 family=binomial(link=logit), data=titanic)
res.logit2
summary(res.logit2)

res.logit3 <- glm(ueberleben ~ geschlecht, 
                  family=binomial(link=logit), data=titanic)
res.logit3
summary(res.logit3)

 # vergleiche modelle
anova(res.logit, res.logit2, test ="Chisq")
anova(res.logit, res.logit3, test ="Chisq")




