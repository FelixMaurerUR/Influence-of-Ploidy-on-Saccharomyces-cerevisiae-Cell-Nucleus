#######################################################################################################################################
#Beta-Mischmodell für Chromatin Coverage (Biomni)

#Packages installieren
 #install.packages("readxl")
 #install.packages("tidyr")
 #install.packages("dplyr")
 #install.packages("ggplot2")
 #install.packages("glmmTMB")
 #install.packages("DHARMa")
 #install.packages("ggeffects")
 #install.packages("MuMIn")
 #install.packages("emmeans")
 #install.packages("showtext")
 #install.packages("scales")
 #install.packages("ggbeeswarm")
 #install.packages("performance")


#Packages laden
 library(readxl)
 library(tidyr)
 library(dplyr)
 library(ggplot2)
 library(glmmTMB)
 library(DHARMa)
 library(ggeffects)
 library(MuMIn)
 library(emmeans)
 library(showtext)
 library(scales)
 library(ggbeeswarm)
 library(performance)


#Schriftart Jost laden
 font_add(family = "jost", regular = "C:/Users/LocalAdmin/Jost/static/Jost-Regular.ttf", bold = "C:/Users/LocalAdmin/Jost/static/Jost-Bold.ttf")
 showtext_auto() 


#Excel Tabelle einlesen
 df <- read_excel("C:/Users/LocalAdmin/Desktop/Felix Maurer/Auswertung/Zellkernvolumen/Zellkernvolumen - Python/Auswertung - Projekt 2/Volumes_all_strains.xlsx",
                 sheet = "ratio_percent" #Falls nur ein Sheet ausgewählt werden soll
                                         #ratio_percent
)


#Daten vorbereiten
#Daten in langes Format bringen
 df_long <- df %>%
   pivot_longer(
     cols      = everything(),
     names_to  = "Stamm",
     values_to = "ChromatinCoverage_percent"
   ) %>%
   filter(!is.na(ChromatinCoverage_percent)) #%>%
 #  group_by(Stamm) %>%              #Nur 200 Zellen verwenden pro Stamm
 #  slice_head(n = 200) %>%          #Nur 200 Zellen verwenden pro Stamm
 #  ungroup()                        #Nur 200 Zellen verwenden pro Stamm


#Text-Fehler korrigieren, Whitespace entfernen "W21605 " (Code eigentlich irrelevant)
 df_long$Stamm <- trimws(df_long$Stamm)


#Ploidie-Zuordnung
 ploidiemap <- data.frame(
   Stamm = c(
     "W21475",
     "W21476",
     "W21477",
     "W21478",
     "W21620",
     "W21656",
     "W21657",
     "W21658",
     "W21691",
     "W21665",
     "W21695",
     "W21692"
     ),
   Ploidie = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4)
)

 df_long <- df_long %>%
   left_join(ploidiemap, by = "Stamm")


#Werte außerhalb 0-100 % ausschließen (biologisch unmöglich) und Ausschluss dokumentieren
 n_total <- nrow(df_long)

 excluded <- df_long %>%
   filter(ChromatinCoverage_percent < 0 | ChromatinCoverage_percent > 100)

 n_excluded <- nrow(excluded)


#Ausgabe in Konsole wie viele Werte ausgeschlossen werden
 cat(sprintf(
  "Ausgeschlossene Werte außerhalb 0-100 %%: %d von %d (%.2f %%)\n",
   n_excluded,
   n_total,
   100 * n_excluded / n_total
))


#Ausschluss pro Ploidiestufe ausgeben
 if (n_excluded > 0) {
   print(
     excluded %>%
       count(Ploidie, name = "n_excluded")
   )
}


#Messfehler endgültig ausschließen
 df_long <- df_long %>%
   filter(
     ChromatinCoverage_percent >= 0,    #Nur biologisch mögliche Werte zwischen 0 und 100 % behalten
     ChromatinCoverage_percent <= 100   #Nur biologisch mögliche Werte zwischen 0 und 100 % behalten
   ) %>%
   mutate(ChromatinCoverage = ChromatinCoverage_percent / 100)   #Prozentwerte auf Bereich 0-1 umrechnen


#Separate Faktorvariable erstellen
 df_long <- df_long %>%
   mutate(Ploidie_factor = factor(Ploidie))

 
#Zusammenfassung
 print(summary(df_long$ChromatinCoverage_percent))
 print(range(df_long$ChromatinCoverage_percent))   #Wertebereich in Prozent
 print(range(df_long$ChromatinCoverage))           #Wertebereich nach Umrechnung auf 0-1


#Prüfen, ob exakt 0% und/oder 100% vorkommt
 n_zero <- sum(df_long$ChromatinCoverage == 0)     #Prüfen, ob exakt 0 % vorkommen
 n_one  <- sum(df_long$ChromatinCoverage == 1)     #Prüfen, ob exakt 100 % vorkommen
 cat("Exact 0% values:", n_zero, "\n")
 cat("Exact 100% values:", n_one, "\n")


#Standard-Beta-Regression ist nur für Werte strikt zwischen 0 und 1 definiert
#Smithson-Verkuilen-Transformation drückt exakte 0/1-Werte leicht ins Innere (falls Werte genau bei 0% oder 100% liegen)
#if (n_zero > 0 | n_one > 0) {
#  n <- nrow(df_long)
#  df_long <- df_long %>%
#    mutate(ChromatinCoverage = (ChromatinCoverage * (n - 1) + 0.5) / n)
#  cat("Exakte 0%/100%-Werte vorhanden -> Smithson-Verkuilen-Transformation angewendet.\n")
#}

 
####################
#Beta-Modelle fitten
####################
 
#Beta-Modell mit linearem Ploidie-Effekt
 model_beta_lin <- glmmTMB(ChromatinCoverage ~ Ploidie + (1 | Stamm), data = df_long, family = beta_family(link = "logit"))

#Beta-Modell mit Ploidie als Faktor
 model_beta_fac <- glmmTMB(ChromatinCoverage ~ Ploidie_factor + (1 | Stamm), data = df_long, family = beta_family(link = "logit"))

#Beta-Modell lineares Modell inklusive stammspezifischer Dispersion
 model_beta_lin_disp <- glmmTMB(ChromatinCoverage ~ Ploidie + (1 | Stamm), data = df_long, family = beta_family(link = "logit"), dispformula = ~ Stamm)
 
#Beta-Modell Faktor-Modell inklusive stammspezifischer Dispersion
 model_beta_fac_disp <- glmmTMB(ChromatinCoverage ~ Ploidie_factor + (1 | Stamm), data = df_long, family = beta_family(link = "logit"), dispformula = ~ Stamm)
 
 
 print(summary(model_beta_lin))      #Modell Zusammenfassung Lineares Modell
 print(summary(model_beta_fac))      #Modell Zusammenfassung Faktor Modell
 print(summary(model_beta_lin_disp)) #Modell Zusammenfassung Lineares Modell mit stammspezifischer Dispersion
 print(summary(model_beta_fac_disp)) #Modell Zusammenfassung Faktor Modell mit stammspezifischer Dispersion

 
#Likelihood-Ratio-Test (Modellvergleich)
#Linearer vs. Faktor-Ploidieeffekt
lrt_beta <- anova(model_beta_lin, model_beta_fac)
print(lrt_beta)

#Linearer vs. Faktor-Plodieeffekt mit stammspezifischer Dispersion
lrt_beta_disp <- anova(model_beta_lin_disp, model_beta_fac_disp)
print(lrt_beta_disp)

#Verbessert stammspezfische Dispersion das lineare Modell?
lrt_disp_lin <- anova(model_beta_lin, model_beta_lin_disp)
print(lrt_disp_lin)

#Verbessert stammspezfische Dispersion das Faktor-Modell?
lrt_disp_fac <- anova(model_beta_fac, model_beta_fac_disp)
print(lrt_disp_fac)


#AIC aller vier Kandidatenmodelle
print(AIC(model_beta_lin, model_beta_fac, model_beta_lin_disp, model_beta_fac_disp))


#p-Werte der Modellvergleiche
p_lrt_beta <- lrt_beta$`Pr(>Chisq)`[2]

p_lrt_beta_disp <- lrt_beta_disp$`Pr(>Chisq)`[2]

p_lrt_disp_lin <- lrt_disp_lin$`Pr(>Chisq)`[2]

p_lrt_disp_fac <- lrt_disp_fac$`Pr(>Chisq)`[2]


#############
#Diagnoselots 
#############

set.seed(123)

#Simulierte Residuen für Modell 1
diagnostics_beta_fac_disp <- simulateResiduals(fittedModel = model_beta_fac_disp, n = 1000)       #An Modelauswahl anpassen

#Simulierte Residuen für Modell 2
diagnostics_beta_lin_disp <- simulateResiduals(fittedModel = model_beta_lin_disp, n = 1000)  #An Modelauswahl anpassen

#Skalierte DHARMa-Residuen extrahieren
#(Werte zwischen 0 und 1, bei korrekt spezifiziertem Modell uniform verteilt)
res_lin <- diagnostics_beta_fac_disp$scaledResiduals                                       #An Modelauswahl anpassen
res_fac <- diagnostics_beta_lin_disp$scaledResiduals                                  #An Modelauswahl anpassen

#Vier Diagnoseplots
par(mfrow = c(2, 2))
par(
  family = "jost",
  cex.lab = 1.1,
  cex.axis = 1.1,
  cex.main = 1.1,
  font.lab = 1,
  font.main = 1,
  lwd = 1.2,
  mar = c(3.5, 4.2, 1.5, 0.8),
  oma = c(0, 0, 3, 0),             #Zusätzlicher Rand oben
  mgp = c(2.0, 0.5, 0),
  bty = "o"
)

#Modell 1 Residuen vs. Fitted
plot(
  fitted(model_beta_fac_disp),                              #An Modelauswahl anpassen     
  res_lin,
  xlab = "",                                            #Leere "" unterdrückt Achsenbeschriftung
  ylab = "",                                            #Leere "" unterdrückt Achsenbeschriftung
  main = "Residuals vs. Fitted (factor beta dispersion model)",    #An Modelauswahl anpassen
  pch = 16,
  col = rgb(0.15, 0.3, 0.6, 0.45),
  cex = 1.1,
  las = 1
)

title(xlab = "Fitted values", mgp = c(1.8, 0.5, 0))   #X-Achsenbeschriftung weg schieben (1.8 Default)
title(ylab = "Residuals", mgp = c(2.5, 0.5, 0))       #Y-Achsenbeschriftung wegschieben (2.0 Default)

abline(h = 0.5, lwd = 2, col = "darkred", lty = 2)            #Erwarteter Median der skalierten Residuen
abline(h = c(0.25, 0.75), lwd = 1, col = "darkred", lty = 3)  #Quartilslinien (optional, bei Bedarf auskommentieren)


#Modell 1 Q-Q Plot (gegen Uniformverteilung)
qqplot(
  qunif(ppoints(length(res_lin))),
  res_lin,
  xlab = "",                                                    #Leere "" unterdrückt Achsenbeschriftung
  ylab = "",                                                    #Leere "" unterdrückt Achsenbeschriftung
  main = "Uniform Q-Q Plot (factor beta dispersion model)",     #An Modelauswahl anpassen
  pch = 16,
  cex = 1.1,
  col = rgb(0.15, 0.3, 0.6, 0.35),
  las = 1
)

title(xlab = "Theoretical Quantiles", mgp = c(1.8, 0.5, 0))   #X-Achsenbeschriftung weg schieben (1.8 Default)
title(ylab = "Sample Quantiles", mgp = c(2.5, 0.5, 0))        #Y-Achsenbeschriftung wegschieben (2.0 Default)

abline(a = 0, b = 1, lwd = 2, col = "darkred")   #1:1-Linie (Analogon zu qqline bei uniformen Residuen)


#Modell 2 Residuen vs. Fitted
plot(
  fitted(model_beta_lin_disp),                        #An Modelauswahl anpassen
  res_fac,
  xlab = "",
  ylab = "",
  main = "Residuals vs. Fitted (linear beta dispersion model)",  #An Modelauswahl anpassen
  pch = 16,
  cex = 1.1,
  col = rgb(0.15, 0.3, 0.6, 0.45),
  las = 1
)

title(xlab = "Fitted values", mgp = c(1.8, 0.5, 0))
title(ylab = "Residuals", mgp = c(2.5, 0.5, 0))

abline(h = 0.5, lwd = 2, col = "darkred", lty = 2)
abline(h = c(0.25, 0.75), lwd = 1, col = "darkred", lty = 3)


#Modell 2 Q-Q Plot (gegen Uniformverteilung)
qqplot(
  qunif(ppoints(length(res_fac))),
  res_fac,
  xlab = "",
  ylab = "",
  main = "Uniform Q-Q Plot (linear beta dispersion model)",     #An Modelauswahl anpassen
  pch = 16,
  cex = 1.1,
  col = rgb(0.15, 0.3, 0.6, 0.35),
  las = 1
)

title(xlab = "Theoretical Quantiles", mgp = c(1.8, 0.5, 0))
title(ylab = "Sample Quantiles", mgp = c(2.5, 0.5, 0))

abline(a = 0, b = 1, lwd = 2, col = "darkred")

#Gemeinsamer Titel
mtext(
  "Model diagnostics for chromatin coverage analysis - Series 1",
  outer = TRUE,
  side = 3,
  line = 1,
  adj = 0.5,
  cex = 1.5,
  font = 1,
  family = "jost"
)


################################################################
#Qualitätstests (plot = FALSE verhindert zusätzliche Test-Plots)
################################################################
#Lineares Modell
print(testUniformity(diagnostics_beta_lin_disp, plot = FALSE))        #An Modelauswahl anpassen
print(testDispersion(diagnostics_beta_lin_disp, plot = FALSE))        #An Modelauswahl anpassen

#Faktormodell
print(testUniformity(diagnostics_beta_fac_disp, plot = FALSE))        #An Modelauswahl anpassen
print(testDispersion(diagnostics_beta_fac_disp, plot = FALSE))        #An Modelauswahl anpassen


#####################################################################################################################################
#Estimated Marginal Means (emmeans) berechnen (modellbasierte geschätzte mittlere Chromatin Coverage für die einzelnen Ploidiestufen)
#####################################################################################################################################

#Beobachtete Ploidiestufen
ploidie_werte <- sort(unique(df_long$Ploidie))

#Emmeans des Faktormodells ohne Dispersion 
emm_beta_fac <- emmeans(model_beta_fac, specs = "Ploidie_factor")     
print(emm_beta_fac)                                                #Werte auf Logit-Skala
print(emm_beta_fac, type = "response")                             #Werte rücktransformiert auf 0-1-Skala (x 100 = Prozent)

#Emmeans des Faktormodells mit Dispersion 
emm_beta_fac_disp <- emmeans(model_beta_fac_disp, specs = "Ploidie_factor")
print(emm_beta_fac_disp)                                           #Werte auf Logit-Skala
print(emm_beta_fac_disp, type = "response")                        #Werte rücktransformiert auf 0-1-Skala (x 100 = Prozent)
 
#Emmeans des linearen Modells ohne Dispersion 
emm_beta_lin <- emmeans(model_beta_lin, specs = "Ploidie", at = list(Ploidie = ploidie_werte))
print(emm_beta_lin)                                                #Werte auf Logit-Skala
print(emm_beta_lin, type = "response")                             #Werte rücktransformiert auf 0-1-Skala (x 100 = Prozent)

#Emmeans des linearen Modells mit Dispersion 
emm_beta_lin_disp <- emmeans(model_beta_lin_disp, specs = "Ploidie", at = list(Ploidie = ploidie_werte))
print(emm_beta_lin_disp)                                           #Werte auf Logit-Skala
print(emm_beta_lin_disp, type = "response")                        #Werte rücktransformiert auf 0-1-Skala (x 100 = Prozent)


#################################################
#Post-hoc paarweise Vergleiche (Tukey-korrigiert)
#################################################
#Fatormodell ohne stammspezifische Dispersion
pairs_beta_fac <- pairs(emm_beta_fac, adjust = "tukey")
print(pairs_beta_fac)                                                  #Vergleiche auf Logit-Skala
print(pairs_beta_fac, type = "response")                               #Rücktransformierte Vergleiche als Odds Ratios

#Faktormodell mit stammspezifischer Dispersion
pairs_beta_fac_disp <- pairs(emm_beta_fac_disp, adjust = "tukey")
print(pairs_beta_fac_disp)                                             #Vergleiche auf Logit-Skala
print(pairs_beta_fac_disp, type = "response")                          #Rücktransformierte Vergleiche als Odds Ratios
 

######################
#Vorhersagen für Plots
######################
#Falls Faktor-Beta-Modell ohne stammspezifischer Dispersion ausgewählt wurde
pred_beta_fac <- ggpredict(model_beta_fac, terms = "Ploidie_factor", type ="fixed", bias_correction = FALSE)

#Falls Faktor-Beta-Modell mit stammspezifischer Dispersion ausgewählt wurde
pred_beta_fac_disp <- ggpredict(model_beta_fac_disp, terms = "Ploidie_factor", type = "fixed", bias_correction = FALSE)

#Falls lineares Beta-Modell ohne stammspezifischer Dispersion ausgewählt wurde
pred_beta_lin <- ggpredict(model_beta_lin, terms = "Ploidie [all]", type = "fixed", bias_correction = FALSE)
ploidy_trend_lin <- emtrends(model_beta_lin, specs = ~1, var = "Ploidie")
ploidy_trend_summary_lin <- as.data.frame(summary(ploidy_trend_lin, infer = c(TRUE, TRUE)))   #Steigung, 95%-KI und p-Wert
print(ploidy_trend_summary_lin)                                                               #Steigung, 95%-KI und p-Wert

#Falls lineares Beta-Modell mit stammspezifischer Dispersion ausgewählt wurde
pred_beta_lin_disp <- ggpredict(model_beta_lin_disp, terms = "Ploidie [all]", type = "fixed", bias_correction = FALSE)
ploidy_trend_lin_disp <- emtrends(model_beta_lin_disp, specs = ~ 1, var = "Ploidie")
ploidy_trend_summary_lin_disp <- as.data.frame(summary(ploidy_trend_lin_disp, infer = c(TRUE, TRUE)))   #Steigung, 95%-KI und p-Wert
print(ploidy_trend_summary_lin_disp)                                                                    #Steigung, 95%-KI und p-Wert

#Ploidie-Beschriftungen
ploidie_labels <- paste0(ploidie_werte, "N")


##########################
#Vorarbeit für beide Plots
##########################

#Median der einzelnen Stämme berechnen
strain_medians <- df_long %>%
  group_by(Stamm, Ploidie) %>%
  summarise(
    Median = median(ChromatinCoverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
#Stämme entsprechend der Ploidie-Map sortieren
  mutate(
    Stamm = factor(
      Stamm,
      levels = ploidiemap$Stamm
    )
  ) %>%
  
  arrange(Ploidie, Stamm) %>%
  group_by(Ploidie) %>%
  
#Symmetrische horizontale Position innerhalb jeder Ploidiegruppe
  mutate(
    n_straens = n(),
    
    x_position = ifelse(
      n_straens == 1,
      Ploidie,
      Ploidie +
        ((row_number() - 1) / (n_straens - 1) - 0.5) * 0.28
    )
  ) %>%
  
  ungroup()


#Vorhandene Stämme in der Reihenfolge der Ploidie-Map
strain_order <- ploidiemap$Stamm[
  ploidiemap$Stamm %in% as.character(strain_medians$Stamm)
]


#Pool verschiedener ggplot-Symbole
shape_pool <- c(
  15,       # Quadrat, gefüllt
  16,       # Kreis, gefüllt
  17,       # Dreieck Spitze oben, gefüllt
  25,       # Dreieck Spitze unten, gefüllt
  18,       # Raute, gefüllt
  5,        # Raute, ohne Füllung
  0,        # Quadrat, ohne Füllung
  1,        # Kreis, ohne Füllung
  2,        # Dreieck Spitze oben, ohne Füllung
  6,        # Dreieck Spitze unten, ohne Füllung
  4,        # Multiplikationszeichen
  3,        # Pluszeichen
  10,       # Quadrat mit Pluszeichen
  12,       # Kreis mit Pluszeichen
  9,        # Raute mit Pluszeichen
  13        # Kreis mit Kreuz
)


#Prüfen, ob genügend Symbole vorhanden sind
if (length(strain_order) > length(shape_pool)) {
  stop("Not enough unique symbols available for all strains.")
}


#Jedem Stamm automatisch ein Symbol zuweisen
strain_shapes <- setNames(
  shape_pool[seq_along(strain_order)],
  strain_order
)


###############################################
#Plot - Beta-Modell mit linearem Ploidie-Effekt
###############################################

#Geschätzte Estimated Marginal Means 
pred_beta_points_lin_disp <- as.data.frame(summary(emm_beta_lin_disp, type = "response", infer = c(TRUE, TRUE)))
print(pred_beta_points_lin_disp)
names(pred_beta_points_lin_disp)

#Feines Raster zwischen kleinster und größter Ploidiestufe
ploidie_grid <- seq(
  min(df_long$Ploidie),
  max(df_long$Ploidie),
  length.out = 200
)

#Modellvorhersagen entlang des linearen Ploidietrends
emm_beta_curve_lin_disp <- emmeans(
  model_beta_lin_disp,
  specs = "Ploidie",
  at = list(Ploidie = ploidie_grid)
)

#Rücktransformation auf die ursprüngliche 0-1-Skala
pred_beta_curve_lin_disp <- as.data.frame(
  summary(
    emm_beta_curve_lin_disp,
    type = "response",
    infer = c(TRUE, TRUE)
  )
)

#CI-Spalten vereinheitlichen
if ("asymp.LCL" %in% names(pred_beta_curve_lin_disp)) {
  pred_beta_curve_lin_disp$lower <- pred_beta_curve_lin_disp$asymp.LCL
  pred_beta_curve_lin_disp$upper <- pred_beta_curve_lin_disp$asymp.UCL
}

if ("lower.CL" %in% names(pred_beta_curve_lin_disp)) {
  pred_beta_curve_lin_disp$lower <- pred_beta_curve_lin_disp$lower.CL
  pred_beta_curve_lin_disp$upper <- pred_beta_curve_lin_disp$upper.CL
}


#Plot
p_beta_lin_disp <- ggplot() +                     #An Modellauswahl anpassen!
  
 geom_point(
    data = df_long,
    aes(
      x = Ploidie,
      y = ChromatinCoverage,
      color = factor(Ploidie),
      group = factor(Ploidie)
    ),
    shape = 1,
    alpha = 1.0,     #Werte von 0-1 (1 voll deckend, 0 voll transparent)
    size = 3.0,
    position = ggbeeswarm::position_quasirandom(
    width = 0.3,
    method = "quasirandom"
    )
  ) +
 
  #95%-Konfidenzband der modellgeschätzten Chromatin Coverage
  geom_ribbon(
    data = pred_beta_curve_lin_disp,
    aes(
      x = Ploidie,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.5            #Transparenz des Konfidenzbandes
  ) +
  
  #Modellierter Ploidietrend
  geom_line(
    data = pred_beta_curve_lin_disp,
    aes(
      x = Ploidie,
      y = response
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.2,
    group = 1,
    lineend = "round",
    linejoin = "round"
  ) +
  
  #Mediane pro Stamm
# geom_point(
#    data = strain_medians,
#    
#    aes(
#      x = x_position,
#      y = Median,
#      shape = Stamm,
#      color = factor(Ploidie),
#      fill = factor(Ploidie)
#    ),
#    
#    size = 2.5,
#    stroke = 2.0,
#    colour = "black",    #Farbe der Symbole ändern
#    fill = "black"       #Fill der Symbole ändern
#  ) +
  
  #95%-Konfidenzintervalle der modellgeschätzten Werte
  geom_errorbar(
    data = pred_beta_points_lin_disp,
    aes(
      x = Ploidie,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 1.0,
    colour = "black"
  ) +
  
  #Modellgeschätzte mittlere Chromatin Coverage
  geom_point(
    data = pred_beta_points_lin_disp,
    aes(
      x = Ploidie,
      y = response
    ),
    inherit.aes = FALSE,
    shape = 16,
    size = 3,               #Größe des Kreis geschätzten Mittelwerts (Default 3)
    stroke = 0.8,           
    colour = "black"
  ) +
  
 scale_color_manual(
   values = c(
     "1" = "#B7D7F0",   # Haploid
     "2" = "#F2B6B6",   # Diploid
     "3" = "#C5DEC5",   # Triploid
     "4" = "#F2D4B6"    # Tetraploid
   )
) +
  
 scale_fill_manual(
   values = c(
     "1" = "#B7D7F0",   # Haploid
     "2" = "#F2B6B6",   # Diploid
     "3" = "#C5DEC5",   # Triploid
     "4" = "#F2D4B6"    # Tetraploid
  )
) +

#Irgendwas  
  scale_shape_manual(values = strain_shapes) +
  
#X-Achse  
 scale_x_continuous(breaks = ploidie_werte, labels = ploidie_labels) +

#Y-Achse    
 scale_y_continuous(
   breaks = seq(0, 1, by = 0.2),
   labels = scales::label_percent(accuracy = 1),
   limits = c(0, 1),
   expand = expansion(mult = c(0, 0.05))
) +

 labs(
   x = "Ploidy level",
   y = "Chromatin coverage (%)"
) +

 theme_bw() +
  
  theme(
    panel.grid = element_blank(),                         #Verhindert die Doppelung der X-Achse
    panel.background = element_blank(),                   #Verhindert die Doppelung der X-Achse
    text = element_text(family = "jost"),                 #Schriftart auf Jost ändern
    axis.title = element_text(family = "jost"),           #Schriftart auf Jost ändern
    axis.text = element_text(family = "jost"),            #Schriftart auf Jost ändern
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.title.x = element_text(face = "plain",
      size = 20,
      margin = margin (t = 12)
    ),
    axis.title.y = element_text(face = "plain",
      size = 20,
      margin = margin(r = 10)
    ),
    
    axis.text.x = element_text(
      size = 19,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 19,
      colour = "black"
    ),
    legend.position = "none"
  )

print(p_beta_lin_disp)


##########################################
#Plot - Beta-Modell mit Ploidie als Faktor
##########################################

#Ploidie-Werte der Predictions für Plot in numerische Werte umwandeln
pred_beta_fac_disp_plot <- as.data.frame(pred_beta_fac_disp) %>%
  mutate(
    Ploidie = as.numeric(as.character(x))
  )
 
p_beta_fac_disp <- ggplot() +                                                   #An Modellauswahl anpassen!
  
  
  #Einzelne Zellmessungen
  geom_point(
    data = df_long,
    aes(
      x = Ploidie,
      y = ChromatinCoverage,
      color = factor(Ploidie),
      group = factor(Ploidie)
    ),
    
    shape = 1,
    alpha = 1.0,                     #Zahl zwischen 0-1 (1 entspricht voll deckend, 0 entspricht voll transparent)
    size = 3.0,
    
    position = ggbeeswarm::position_quasirandom(
      width = 0.3,
      method = "quasirandom"
    )
  ) +
  
  
  #Mediane pro Stamm
  geom_point(
    data = strain_medians,
    
    aes(
      x = x_position,
      y = Median,
      shape = Stamm,
      color = factor(Ploidie),
      fill = factor(Ploidie)
    ),
    
    size = 2.5
    #stroke = 2.0,
    #colour = "black",
    #fill = "black"
  ) +
  
  
  #Geschätztes 95%-Konfidenzintervall
  geom_errorbar(
    data = pred_beta_fac_disp_plot,              #An Modellauswahl anpassen!
    aes(
      x = Ploidie,
      ymin = conf.low,
      ymax = conf.high
    ),
    
    width = 0.08,       #Breite der Enden des Konfidenzintervalls
    linewidth = 1.0,    #Dicke der Linie
    color = "black"
  ) +
  
  
  #Geschätzter Mittelwert des Faktor-Modells
  geom_point(
    data = pred_beta_fac_disp_plot,             #An Modellauswahl anpassen!
    aes(
      x = Ploidie,
      y = predicted
    ),
    
    shape = 21,
    size = 3.0,
    stroke = 0.5,
    
    color = "black",
    fill = "black"
  ) +
  
  
  #Farben der einzelnen Ploidiestufen
  scale_color_manual(
    values = c(
      "1" = "#B7D7F0",   #Haploid
      "2" = "#F2B6B6",   #Diploid
      "3" = "#C5DEC5",   #Triploid
      "4" = "#F2D4B6"    #Tetraploid
    )
  ) +
  
  
  #Füllfarben
  scale_fill_manual(
    values = c(
      "1" = "#B7D7F0",   #Haploid
      "2" = "#F2B6B6",   #Diploid
      "3" = "#C5DEC5",   #Triploid
      "4" = "#F2D4B6"    #Tetraploid
    )
  ) +
  
  
  #Symbole der einzelnen Stämme
  scale_shape_manual(
    values = strain_shapes
  ) +
  
  
  #X-Achse
  scale_x_continuous(
    breaks = ploidie_werte,
    labels = ploidie_labels
  ) +
  
  
  #Y-Achse
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  
  #Achsentitel
  labs(
    x = "Ploidy level",
    y = "Chromatin coverage (%)"
  ) +
  
  
  #Theme
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    text = element_text(
      family = "jost"
    ),
    
    axis.title = element_text(
      family = "jost"
    ),
    
    axis.text = element_text(
      family = "jost"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 20,
      margin = margin(t = 12)
    ),
    
    axis.title.y = element_text(
      face = "plain",
      size = 20,
      margin = margin(r = 10)
    ),
    
    axis.text.x = element_text(
      size = 19,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 19,
      colour = "black"
    ),
    
    legend.position = "none"
  )

print(p_beta_fac_disp)


#####################################################
#Säulendiagramm (falls lineares Modell gewählt wurde)
#####################################################

#Emmeans des Linearen-Modells auf die 0-1-Skala bringen
emm_resp <- regrid(emm_beta_lin_disp, transform = "response") #An Modellauswahl anpassen!

#Kontraste gegenüber 1N (Differenzen der geschätzten Mittelwerte, mvt-korrigiert)
diff_contrasts <- contrast(emm_resp, method = "trt.vs.ctrl", ref = 1)
diff_summary <- summary(diff_contrasts, infer = c(TRUE, TRUE), adjust = "mvt")

#Vorhandene Ploidiestufen außer Referenz 1N bestimmen
comparison_ploidies <- sort(setdiff(unique(df_long$Ploidie), 1))

#Differenzen (in Prozentpunkten) und p-Werte
diff_table <- as.data.frame(diff_summary)

#CI-Spalten vereinheitlichen
if ("asymp.LCL" %in% names(diff_table)) {
  diff_table$lower <- diff_table$asymp.LCL
  diff_table$upper <- diff_table$asymp.UCL
}

if ("lower.CL" %in% names(diff_table)) {
  diff_table$lower <- diff_table$lower.CL
  diff_table$upper <- diff_table$upper.CL
}

diff_table <- diff_table %>%
  mutate(
    Ploidie = comparison_ploidies,
    Diff_pp = estimate * 100,
    Lower_pp = lower * 100,
    Upper_pp = upper * 100
  ) %>%
  select(
    Ploidie,
    Diff_pp,
    Lower_pp,
    Upper_pp,
    p_value = p.value
  )

#Referenzgruppe 1N ergänzen (Differenz = 0)
diff_final <- data.frame(
  Ploidie = 1,
  Diff_pp = 0,
  Lower_pp = NA_real_,
  Upper_pp = NA_real_,
  p_value = NA_real_
) %>%
  bind_rows(diff_table) %>%
  mutate(
    Ploidie_Label = paste0(Ploidie, "N")
  )

#Signifikanzsterne
diff_final <- diff_final %>%
  mutate(
    significance = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(diff_final)

#Odds Ratios der gleichen Kontraste als Konsolen-Ausgabe
cat("Odds Ratios vs. 1N (Logit-Skala des Faktormodells):\n")
print(
  contrast(emm_beta_lin_disp, method = "trt.vs.ctrl", ref = 1), #An Modellauswahl anpassen
  type = "response"
)

#Anzahl der Vergleiche
n_comparisons <- length(comparison_ploidies)

#Wertebereich inklusive 95%-Konfidenzintervalle bestimmen
y_max <- max(
  0,
  diff_final$Diff_pp,
  diff_final$Upper_pp,
  na.rm = TRUE
)

y_min <- min(
  0,
  diff_final$Diff_pp,
  diff_final$Lower_pp,
  na.rm = TRUE
)

y_span <- y_max - y_min

if (y_span == 0) {
  y_span <- 1
}

bracket_y <- seq(
  from = y_max + 0.20 * y_span,
  by = 0.15 * y_span,
  length.out = n_comparisons
)

#Beschriftungspositionen am Anfang der Säulen nahe der Nulllinie
diff_final <- diff_final %>%
  mutate(
    label_y = case_when(
      Diff_pp > 0  ~  0.055 * y_span,
      Diff_pp < 0  ~ -0.055 * y_span,
      TRUE         ~  0.055 * y_span
    )
  )

#Tabelle für Signifikanzklammern
brackets <- data.frame(
  x_start = rep(1, n_comparisons),
  x_end = comparison_ploidies,
  y = bracket_y,
  significance = diff_final %>%
    filter(Ploidie != 1) %>%
    arrange(Ploidie) %>%
    pull(significance)
)

print(brackets)

#Plot
p_diff <- ggplot() +

#Nulllinie (Referenz)
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.6
  ) +

#Säulen
  geom_col(
    data = diff_final,
    aes(
      x = Ploidie,
      y = Diff_pp,
      fill = Ploidie_Label
    ),
    color = "black",       #Rahmenlinie um jede Säule
    width = 0.6,
    linewidth = 0.6
  ) +
  
  geom_errorbar(
    data = diff_final %>% filter(Ploidie != 1),
    aes(
      x = Ploidie,
      ymin = Lower_pp,
      ymax = Upper_pp
    ),
    width = 0.08,
    linewidth = 1.0,
    colour = "black"
  ) +

#Differenz-Wert an jeder Säule (oberhalb positiver, unterhalb negativer Säulen)
  
  geom_text(
    data = diff_final,
    aes(
      x = Ploidie,
      y = label_y,
      label = sprintf("%.2f", Diff_pp)
    ),
    family = "jost",
    size = 7.0,
    color = "black",
    vjust = 0.5
  ) +

#Horizontale Linien der Signifikanzklammern
  geom_segment(
    data = brackets,
    aes(
      x = x_start,
      xend = x_end,
      y = y,
      yend = y
    ),
    color = "black",
    linewidth = 0.6
  ) +

#Linkes Ende der Signifikanzklammern
  geom_segment(
    data = brackets,
    aes(
      x = x_start,
      xend = x_start,
      y = y,
      yend = y - 0.03 * y_span
    ),
    color = "black",
    linewidth = 0.6
  ) +

#Rechtes Ende der Signifikanzklammern
  geom_segment(
    data = brackets,
    aes(
      x = x_end,
      xend = x_end,
      y = y,
      yend = y - 0.03 * y_span
    ),
    color = "black",
    linewidth = 0.6
  ) +

#Signifikanzangaben über den Klammern
  geom_text(
    data = brackets,
    aes(
      x = (x_start + x_end) / 2,
      y = y + 0.02 * y_span,
      label = significance
    ),
    family = "jost",
    size = 5.0,
    color = "black",
    vjust = 0
  ) +

#X-Achse
  scale_x_continuous(breaks = diff_final$Ploidie, labels = diff_final$Ploidie_Label) +

#Y-Achse (unterer Rand mit Platz für Beschriftung negativer Säulen)
  scale_y_continuous(limits = c(y_min - 0.18 * y_span, max(bracket_y) + 0.15 * y_span),
    expand = expansion(mult = c(0, 0.02))
  ) +

#Farben der Ploidiestufen
  scale_fill_manual(
    values = c(
      "1N" = "#B7D7F0",
      "2N" = "#F2B6B6",
      "3N" = "#C5DEC5",
      "4N" = "#F2D4B6"
    )
  ) +

  labs(
    x = "Ploidy level",
    y = "Difference in chromatin coverage \n(percentage points)"
  ) +

  theme_bw() +                                      #Löst Problem mit doppelten Rahmenlinien

  theme(
    text = element_text(family = "jost"),

    axis.title = element_text(family = "jost"),

    axis.text = element_text(family = "jost"),

    #Keine Gitternetzlinien
    panel.grid = element_blank(),                    #Löst Problem mit doppelten Rahmenlinien

    #Panel-Hintergrund entfernen
    panel.background = element_blank(),              #Löst Problem mit doppelten Rahmenlinien

    #Keine zusätzlichen Achsenlinien zeichnen 
    axis.line = element_blank(),                     #Löst Problem mit doppelten Rahmenlinien

    #X-Achsenwerte
    axis.text.x = element_text(size = 19, colour = "black"),

    #Y-Achsenwerte
    axis.text.y = element_text(size = 19, colour = "black"),

    #X-Achsentitel
    axis.title.x = element_text(size = 22, margin = margin(t = 12), face = "plain"),

    #Y-Achsentitel
    axis.title.y = element_text(size = 22, margin = margin(r = 10), face = "plain"),

    #Nur ein geschlossener Rahmen
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),

    #Keine Legende
    legend.position = "none"
  )

print(p_diff)


###################################################
#Säulendiagramm (falls Faktor Modell gewählt wurde)
###################################################

#Das für die Darstellung verwendete Modell auswählen
emm_for_diff <- emm_beta_fac_disp                       #An Modellauswahl anpassen!

#EMMs auf die ursprüngliche 0-1-Skala bringen
emm_resp <- regrid(
  emm_for_diff,
  transform = "response"
)

#Kontraste gegenüber 1N
diff_contrasts <- contrast(
  emm_resp,
  method = "trt.vs.ctrl",
  ref = 1
)

#Zusammenfassung mit multipler Testkorrektur
diff_summary <- summary(diff_contrasts, infer = c(TRUE, TRUE), adjust = "mvt")

#Vorhandene Ploidiestufen außer Referenz 1N
comparison_ploidies <- sort(setdiff(unique(df_long$Ploidie), 1))

#Differenzen in Prozentpunkten und p-Werte
diff_table <- as.data.frame(diff_summary) %>%
  mutate(
    Ploidie = comparison_ploidies,
    Diff_pp = estimate * 100
  ) %>%
  select(
    Ploidie,
    Diff_pp,
    p_value = p.value
  )

#Referenzgruppe 1N ergänzen
diff_final <- data.frame(
  Ploidie = 1,
  Diff_pp = 0,
  p_value = NA_real_
) %>%
  bind_rows(diff_table) %>%
  mutate(
    Ploidie_Label = paste0(Ploidie, "N")
  )

#Signifikanzsterne
diff_final <- diff_final %>%
  mutate(
    significance = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(diff_final)

y_max <- max(diff_final$Diff_pp)
y_min <- min(0, min(diff_final$Diff_pp))

y_span <- y_max - y_min

if (y_span == 0) {
  y_span <- 1
}


diff_final <- diff_final %>%
  mutate(
    label_y = ifelse(
      Diff_pp >= 0,
      Diff_pp + 0.025 * y_span,
      Diff_pp - 0.025 * y_span
    ),
    
    label_vjust = ifelse(
      Diff_pp >= 0,
      0,
      1
    )
  )

n_comparisons <- length(comparison_ploidies)


bracket_y <- seq(
  from = y_max + 0.20 * y_span,
  by = 0.16 * y_span,
  length.out = n_comparisons
)


brackets <- data.frame(
  x_start = rep(1, n_comparisons),
  x_end = comparison_ploidies,
  y = bracket_y,
  
  significance = diff_final %>%
    filter(Ploidie != 1) %>%
    arrange(Ploidie) %>%
    pull(significance)
)


print(brackets)


#Plot

p_diff <- ggplot() +
  
  #Nulllinie
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.6
  ) +
  
  #Säulen
  geom_col(
    data = diff_final,
    aes(
      x = Ploidie,
      y = Diff_pp,
      fill = Ploidie_Label
    ),
    color = "black",
    width = 0.6,
    linewidth = 0.6
  ) +
  
  #Werte an den Säulen
  geom_text(
    data = diff_final,
    aes(
      x = Ploidie,
      y = label_y,
      vjust = label_vjust,
      label = sprintf("%.2f", Diff_pp)
    ),
    family = "jost",
    size = 7.0,
    color = "black"
  ) +
  
  #Horizontale Linien der Signifikanzklammern
  geom_segment(
    data = brackets,
    aes(
      x = x_start,
      xend = x_end,
      y = y,
      yend = y
    ),
    color = "black",
    linewidth = 0.6
  ) +
  
  #Linkes Ende der Signifikanzklammern
  geom_segment(
    data = brackets,
    aes(
      x = x_start,
      xend = x_start,
      y = y,
      yend = y - 0.03 * y_span
    ),
    color = "black",
    linewidth = 0.6
  ) +
  
  #Rechtes Ende der Signifikanzklammern
  geom_segment(
    data = brackets,
    aes(
      x = x_end,
      xend = x_end,
      y = y,
      yend = y - 0.03 * y_span
    ),
    color = "black",
    linewidth = 0.6
  ) +
  
  #Signifikanzangaben
  geom_text(
    data = brackets,
    aes(
      x = (x_start + x_end) / 2,
      y = y + 0.02 * y_span,
      label = significance
    ),
    family = "jost",
    size = 5.0,
    color = "black",
    vjust = 0
  ) +
  
  #X-Achse
  scale_x_continuous(
    breaks = diff_final$Ploidie,
    labels = diff_final$Ploidie_Label
  ) +
  
  #Y-Achse
  scale_y_continuous(
    limits = c(
      y_min - 0.18 * y_span,
      max(bracket_y) + 0.15 * y_span
    ),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  #Farben
  scale_fill_manual(
    values = c(
      "1N" = "#B7D7F0",
      "2N" = "#F2B6B6",
      "3N" = "#C5DEC5",
      "4N" = "#F2D4B6"
    )
  ) +
  
  labs(
    x = "Ploidy level",
    y = "Difference in chromatin coverage\n(percentage points)"
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(family = "jost"),
    axis.title = element_text(family = "jost"),
    axis.text = element_text(family = "jost"),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    
    axis.text.x = element_text(
      size = 19,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 19,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 22,
      margin = margin(t = 12),
      face = "plain"
    ),
    
    axis.title.y = element_text(
      size = 22,
      margin = margin(r = 10),
      face = "plain"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.position = "none"
  )


print(p_diff)

##############################################################################################################################