############################################################################################################################
##Boxplots

#Packages installieren
#install.packages("readxl")
#install.packages("ggplot2")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("ggh4x")
#install.packages("showtext")

#Packages laden
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
library(ggh4x)
library(showtext)

#Schriftart Jost laden
font_add(family = "jost", regular = "C:/Users/LocalAdmin/Jost/static/Jost-Regular.ttf", bold = "C:/Users/LocalAdmin/Jost/static/Jost-Bold.ttf")
showtext_auto() 

#Excel Tabelle einlesen
df <- read_excel("C:/Users/LocalAdmin/Desktop/Felix Maurer/Auswertung/Septin-Ring Durchmesser/Auswertung - Projekt 6 cdc11/Diploide/ring_diameter_per_strain_final.xlsx",
                 #sheet = "histone_volume_um3" #Falls nur ein Sheet ausgewählt werden soll
                                               #nuclear_volume_um3 
                                               #histone_volume_um3
                                               #ratio_percent
                 )

#Stämme auswählen
auswahl <- c(
  "W21503",
  "W21504",
  "W21584",
  "W21585",
  "W21586",
  "W21587"
    )

df <- df %>%
  select(all_of(auswahl))

#Daten in langes Format bringen
df_long <- df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Stamm",
    values_to = "Kernvolumen"
  ) %>%
  filter(!is.na(Kernvolumen)) %>%
#  filter(   
#    Kernvolumen >= 0,       #Nur biologisch mögliche Werte zwischen 0 und 100 % behalten (nur Chromatin Coverage)
#    Kernvolumen <= 100      #Nur biologisch mögliche Werte zwischen 0 und 100 % behalten (nur Chromatin Coverage)
#  ) %>%
  group_by(Stamm) %>%              #Nur 200 Zellen verwenden pro Stamm
  slice_head(n = 200) %>%          #Nur 200 Zellen verwenden pro Stamm
  ungroup()                        #Nur 200 Zellen verwenden pro Stamm
 

#Boxplots verschieben (falls notwendig!!!)
#df_long$Stamm <- factor(
#  df_long$Stamm,
#  levels = c(
#       "W21475 7.5%",
#       "W21476 7.5%",
#       "W21477 7.5%",
#       "W21478 7.5%",
#       "W21475 10%",
#       "W21476 10%",
#       "W21477 10%",
#       "W21478 10%",
#       "W21475 20%",
#       "W21476 20%",
#       "W21477 20%",
#       "W21478 20%"
# )
#)

#Median, Mittelwert und Stichprobengröße berechnen
summary.df <- df_long %>%
  group_by(Stamm) %>%
  summarise(
    Mean = mean(Kernvolumen),
    Median = median(Kernvolumen),
    SD = sd(Kernvolumen),
    n = n(),
    .groups = "drop"
  )


#X-Achsenbeschriftung mit Stichprobengröße erstellen
labels_n <- paste0(summary.df$Stamm, "\n(n=", summary.df$n, ")")

#Falls die x-Achsenbeschriftung gedreht werden muss
#labels_n <- paste0(
#  summary.df$Stamm,
#  " (n=",
#    summary.df$n,
#  ")"
#)


#Boxplot
p <- ggplot(df_long, aes(x = Stamm, y = Kernvolumen)) +

#Datenpunkte im Hintergrund verteilen (Jitter)
#  geom_jitter(
#    color = "#F2D4B6",       #Farbe der Punkte grey85 als Default, #B7D7F0 für Haploide, #F2B6B6 für Diploide, #C5DEC5 für Triploide, #F2D4B6 für Tetraploide
#    alpha = 0.8,             # Transparenz der Punkte (0 = unsichtbar, 1 = deckend)
#    width = 0.2,             # Streubreite der Punkte nach links/rechts
#    height = 0               # Keine vertikale Verzerrung der echten Werte
#  ) +
  
#Schnurrhaare anzeigen lassen    
  stat_boxplot(
    geom = "errorbar",
    width = 0.15,         #Größe der Schnurhaare verändern (Default bei 0.15, Default, falls nur ein Tetraploider 0.1)
    colour = "black"
   ) +
  
  geom_boxplot(
    width = 0.7,          #Wert erhöhen, falls Box zu eng (0.7 als Default, Default, falls nur ein Tetraploider 0.4)
    fill = "#F2B6B6",     #grey85 als Default, #B7D7F0 für Haploide, #F2B6B6 für Diploide, #C5DEC5 für Triploide, #F2D4B6 für Tetraploide
    colour = "black",     #black als Default
#    outlier.shape = NA,    #Falls Datenpunkte angezeigt werden sollen
    outlier.shape = 1,   #Falls Datenpunkte angezeigt werden ausgrauen
    outlier.size = 2.5,  #Falls Datenpunkte angezeigt werden ausgrauen
    fatten = 1.0
  ) +
  
  geom_text(
    data = summary.df,
    aes(
      y = Median,
      label = round(Median, 2)
    ),
    vjust = -0.15, #Median in der Box verschieben, Default: -0.15
    size = 7.0,    #Schriftgröße des Medians erhöhen, Default: 7.0
    family = "jost"
  ) +
  
  scale_x_discrete(labels = labels_n        #X-Achsenbeschriftung mit Stichprobengröße
  ) +

#  scale_x_discrete() +  #X-Achsenbeschriftung ohne Stichprobengröße
  
#  scale_x_discrete(      #X-Achsenebschriftung manuell setzen
#    labels = c(
#      "W21475",
#      "W21476",
#      "W21477",
#      "W21478",
#      "W21475",
#      "W21476",
#      "W21477",
#      "W21478",
#      "W21475",
#      "W21476",
#      "W21477",
#      "W21478"
#    )
#  ) +
  
#Y-Achse für Zellkern- und Histonvolumen    
#  scale_y_continuous(
#    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60),
#    minor_breaks = seq(0, 10, by = 1),
#    limits = c(0, 4) #Zweite Zahl variieren, je nachdem welche Ploidie Stufe untersucht wird
#  ) +

#Y-Achse für Septin Ring Durchmesser
  scale_y_continuous(
    breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4),
    minor_breaks = seq(0, 10, by = 1),
    limits = c(0, 3.5) #Zweite Zahl variieren, je nachdem welche Ploidie Stufe untersucht wird
  ) +
  
#Y-Achse Chromatin Coverage (%)
# scale_y_continuous(
#   breaks = c(0, 20, 40, 60, 80, 100),
#   limits = c(0, 100)
# ) +
      
  labs(
#    title = "ImageJ",        #Titel für Graph hinzufügen
     x = "Strains",                 #Titel X-Achse: Für Haploide und Diploide = NULL für Tri- und Tetraploide = "Strains"
#     y = "Nuclear volume (µm³)"      #Zellkernvolumina
#    y = "Histone volume (µm³) "    #Histonvolumina
     y = "septin ring diameter (µm)"
#     y = "Chromatin Coverage (%)"  #Histon-to-Volume Ratio (%)
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(family = "jost"),       #Schriftart auf Jost ändern
    axis.title = element_text(family = "jost"), #Schriftart auf Jost ändern
    axis.text = element_text(family = "jost"),  #Schriftart auf Jost ändern
    
    plot.title = element_text(family = "jost", size = 28, face = "plain", margin = margin(b = 3), hjust = 0.5), #Titel für Graph hinzufügen
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_text (size = 22, margin = margin(t = 12), face = "plain"),  
    #Schriftgröße Strains verändern und Titel von x-Achse wegschieben, Schriftgröße Default: 15, Wegschiebung Default: 8
    axis.title.y = element_text (size = 22, margin = margin(r = 10), face = "plain"),  
    #Schriftgröße Nuclear volume verändern und Titel von y-Achse wegschieben, Schriftgröße Default: 15, Wegschiebung Default: 8
    axis.text.x = element_text(size = 19, colour = "black"),    #Schriftgröße x-Achse verändern, Default: 14, ausgrauen, falls X-Achse gedreht werden muss
    axis.text.y = element_text(size = 19, colour = "black"),    #Schriftgröße y-Achse verändern, Default: 15
    panel.border = element_rect( 
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
#Falls x-Achsen Beschriftung gedreht werden muss
#   axis.text.x = element_text(
#     angle = 50,
#      hjust = 1,
#      vjust = 1,
#      size  = 14,          #Default 14
#     colour = "black"
#   )
  )    

print(p)
print(summary.df)

#######################################################################################################################################
#Linear Mixed Model für Zellkern- und Histonvolumen (Biomni)
#Packages installieren
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("ggeffects")
#install.packages("MuMIn")
#install.packages("emmeans")
#install.packages("ggbeeswarm")
#install.packages("patchwork")


#Packages laden
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ggeffects)
library(MuMIn)
library(emmeans)
library(showtext)
library(ggbeeswarm)
library(patchwork)


#Schriftart Jost laden
font_add(family = "jost", regular = "C:/Users/LocalAdmin/Jost/static/Jost-Regular.ttf", bold = "C:/Users/LocalAdmin/Jost/static/Jost-Bold.ttf")
showtext_auto() 


#Daten laden und aufbereiten
df <- read_excel("C:/Users/LocalAdmin/Desktop/Felix Maurer/Auswertung/Septin-Ring Durchmesser/Auswertung - Projekt 6 cdc11/Septin_ring_diameter_all_strains.xlsx",
                # sheet = "histone_volume_um3" #Falls nur ein Sheet ausgewählt werden soll
                                         #nuclear_volume_um3
                                         #histone_volume_um3
                 )


#Daten in langes Format bringen
df_long <- df %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "Stamm",
    values_to = "Kernvolumen"
  ) %>%
  filter(!is.na(Kernvolumen)) #%>%
#  group_by(Stamm) %>%              #Nur 200 Zellen verwenden pro Stamm
#  slice_head(n = 200) %>%          #Nur 200 Zellen verwenden pro Stamm
#  ungroup()                        #Nur 200 Zellen verwenden pro Stamm


df_long$Stamm <- trimws(df_long$Stamm)


#Ploidie-Zuordnung
ploidiemap <- data.frame(
  Stamm = c(
    "W21447",
    "W21448",
    "W21449",
    "W21450",
    "W21503",
    "W21504",
    "W21584",
    "W21585",
    "W21586",
    "W21587",
    "W21589",
    "W21592",
    "W21588"
    ),
  Ploidie = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 4)
)

df_long <- df_long %>%
  left_join(ploidiemap, by = "Stamm")


####################
#Vier Modelle fitten
####################

#Roh-Skala
model_raw_lin <- lmer(Kernvolumen ~ Ploidie + (1 | Stamm),data = df_long,REML = TRUE)

model_raw_fac <- lmer(Kernvolumen ~ factor(Ploidie) + (1 | Stamm), data = df_long, REML = TRUE)


#Log-Skala
model_log_lin <- lmer(log(Kernvolumen) ~ Ploidie + (1 | Stamm), data = df_long, REML = TRUE)

model_log_fac <- lmer(log(Kernvolumen) ~ factor(Ploidie) + (1 | Stamm), data = df_long, REML = TRUE)


#Likelihood-Ratio-Test
lrt_raw <- anova(model_raw_lin, model_raw_fac, refit = TRUE)
print(lrt_raw)
p_lrt_raw <- lrt_raw$`Pr(>Chisq)`[2]

lrt_log <- anova(model_log_lin, model_log_fac, refit = TRUE)
print(lrt_log)
p_lrt_log <- lrt_log$`Pr(>Chisq)`[2]


###############
#Diagnose-Plots
###############

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

###########
#Raw-Modell
###########

#Plot 1 - Raw-Modell Residuen vs. Fitted

plot(fitted(model_raw_fac),                               #An Modellwahl anpassen!
     resid(model_raw_fac),                                #An Modellwahl anpassen!
     xlab = "",                                           #Leere "" unterdrückt Achsenbeschriftung
     ylab = "",                                           #Leere "" unterdrückt Achsenbeschriftung
     main = "Residuals vs. Fitted (raw_fac_model)",       #An Modellwahl anpassen!
     pch = 16, col = rgb(0.15, 0.3, 0.6, 0.45), cex = 1.1, las = 1)

title(xlab = "Fitted values", mgp = c(1.8, 0.5, 0))   #X-Achsenbeschriftung weg schieben (1.8 Default)
title(ylab = "Residuals", mgp = c(2.5, 0.5, 0))       #Y-Achsenbeschriftung wegschieben  (2.0 Default)

abline(h = 0, lwd = 2, col = "darkred", lty = 2)


#Plot 2 - Raw-Modell QQ-Plot
qqnorm(
  resid(model_raw_fac),                            #An Modellwahl anpassen!
  pch = 16,
  cex = 1.1,
  col = rgb(0.15,0.3,0.6,0.35),
  main = "Normal Q-Q Plot (raw_fac_model)",        #An Modellwahl anpassen!
  xlab = "",                                       #Leere "" unterdrückt Achsenbeschriftung
  ylab = "",                                       #Leere "" unterdrückt Achsenbeschriftung
  las = 1
)

title(xlab = "Theoretical Quantiles", mgp = c(1.8, 0.5, 0))   #X-Achsenbeschriftung wegschieben (1.8 Default)
title(ylab = "Sample Quantiles", mgp = c(2.5, 0.5, 0))       #Y-Achsenbeschriftung wegschieben (2.0 Default)

qqline(
  resid(model_raw_fac),                           #An Modellwahl anpassen! 
  lwd = 2,
  col = "darkred"
)


###########
#Log-Modell
###########

#Plot 1 - Log-Modell Residuen vs. Fitted

plot(
  fitted(model_log_fac),                          #An Modellwahl anpassen!
  resid(model_log_fac),                           #An Modellwahl anpassen!
  xlab = "",                                      #Leere "" unterdrückt Achsenbeschriftung
  ylab = "",                                      #Leere "" unterdrückt Achsenbeschriftung 
  main = "Residuals vs. Fitted (log_fac_model)",  #An Modellwahl anpassen!
  pch = 16,
  cex = 1.1,
  col = rgb(0.15, 0.3, 0.6, 0.45),
  las = 1
)

title(xlab = "Fitted values (log)", mgp = c(1.8, 0.5, 0))   #X-Achsenbeschriftung weg schieben (1.8 Default)
title(ylab = "Residuals (log)", mgp = c(2.5, 0.5, 0))       #Y-Achsenbeschriftung wegschieben  (2.0 Default)

abline(h = 0, lwd = 2, col = "darkred", lty = 2)


#Plot 2 - Log-Modell QQ-Plot

qqnorm(
  resid(model_log_fac),                                     #An Modellwahl anpassen!
  pch = 16,
  cex = 1.1,
  col = rgb(0.15,0.3,0.6,0.35),
  main = "Normal Q-Q Plot (log_fac_model)",  #An Modellwahl anpassen!
  xlab = "",                                 #Leere "" unterdrückt Achsenbeschriftung
  ylab = "",                                 #Leere "" unterdrückt Achsenbeschriftung
  las = 1
)

title(xlab = "Theoretical Quantiles (log)", mgp = c(1.8, 0.5, 0))   #X-Achsenbeschriftung weg schieben (1.8 Default)
title(ylab = "Sample Quantiles (log)", mgp = c(2.5, 0.5, 0))        #Y-Achsenbeschriftung wegschieben  (2.0 Default)

qqline(
  resid(model_log_fac), #An Modellauswahl anpassen!
  lwd = 2,
  col = "darkred"
)

mtext(
#  "Model diagnostics for nuclear volume analysis - Series 6",        #Zellkernvolumen, Nummer des Projekts ändern!
#  "Model diagnostics for histone volume analysis - Series 2",        #Histonvolumen, Nummer des Projekts ändern!
  "Model diagnostics for septin ring diameter analysis - Series 6",  #Septin-Ring Durchmesser, Nummer des Projekts ändern!
  outer = TRUE,
  side = 3,
  line = 1,
  adj = 0.5,
  cex = 1.5,
  font = 1,
  family = "jost"
)


#########################
#Weitere Diagnosemethoden
#########################

model_final <- model_log_fac   #An Modellauswahl anpassen!

#Test 1 - Test auf Singularität

lme4::isSingular(model_final, tol = 1e-4)
VarCorr(model_final)
print(VarCorr(model_final), comp = c("Variance", "Std.Dev."))


#Test 2 - Normalverteilung der Random Effects

#Random Intercepts der Stämme extrahieren
random_effects <- ranef(model_final)$Stamm[, 1]

#Dataframe für ggplot erstellen
random_effects_df <- data.frame(Random_effect = random_effects)

#Random Effects extrahieren
random_effects <- ranef(model_final)$Stamm[, 1]

#Diagnosedatensatz erstellen
diagnostic_df <- df_long %>%
  mutate(
    Fitted = fitted(model_final),
    Residual = resid(model_final)
  )

#Reihenfolge der Stämme festlegen
strain_levels <- ploidiemap$Stamm[
  ploidiemap$Stamm %in% unique(diagnostic_df$Stamm)
]

diagnostic_df$Stamm <- factor(
  diagnostic_df$Stamm,
  levels = strain_levels
)


##########################
#Zusätzliche Diagnoseplots
##########################

#Einheitliches Theme für alle drei Diagnoseplots
theme_diagnostics <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    text = element_text(
      family = "jost"),
    
    #Titel der einzelnen Panels
    plot.title = element_text(
      family = "jost",
      size = 11,
      face = "plain",
      hjust = 0.5,
      margin = margin(b = 3)
    ),
    
    #X-Achsenwerte
    axis.text.x = element_text(
      family = "jost",
      size = 11,
      colour = "black"
    ),
    
    #Y-Achsenwerte
    axis.text.y = element_text(
      family = "jost",
      size = 11,
      colour = "black"
    ),
    
    #X-Achsentitel
    axis.title.x = element_text(
      family = "jost",
      size = 11,
      face = "plain",
      margin = margin(t = 8)
    ),
    
    #Y-Achsentitel
    axis.title.y = element_text(
      family = "jost",
      size = 11,
      face = "plain",
      margin = margin(r = 10)
    ),
    
    #Geschlossener Rahmen
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )


############################
#Plot 1 - Random Effects Q-Q
############################

p_random_effects <- ggplot(
  random_effects_df,
  aes(sample = Random_effect)
) +
  
  stat_qq(
    shape = 16,
    size = 2.5,
    color = rgb(0.15, 0.3, 0.6, 0.35)
  ) +
  
  stat_qq_line(
    color = "darkred",
    linewidth = 1
  ) +
  
  labs(
    title = "Random effects",
    x = "Theoretical Quantiles",
    y = "Random intercepts"
  ) +
  
  theme_diagnostics


###############################
#Plot 2 - Residuen nach Ploidie
###############################

p_res_ploidy <- ggplot(
  diagnostic_df,
  aes(
    x = factor(Ploidie),
    y = Residual
  )
) +
  
  geom_boxplot(
    width = 0.6,
    colour = "black",
    fill = rgb(0.15, 0.3, 0.6, 0.35),
    staplewidth = 0.3,                          #Whisker einfügen
    outlier.shape = 1,
    outlier.colour = rgb(0.15, 0.3, 0.6, 0.35)
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    color = "darkred"
  ) +
  
  labs(
    title = "Residuals by ploidy",
    x = "Ploidy level",
    y = "Residuals"
  ) +
  
  scale_x_discrete(
    labels = function(x) paste0(x, "N")
  ) +
  
  theme_diagnostics


#############################
#Plot 3 - Residuen nach Stamm
#############################

 p_res_strain <- ggplot(
   diagnostic_df,
   aes(
     x = Stamm,
     y = Residual
   )
 ) +
  
  geom_boxplot(
    width = 0.6,
    colour = "black",
    fill = rgb(0.15, 0.3, 0.6, 0.35),
    staplewidth = 0.3,                   #Whisker einfügen
    outlier.shape = 1,
    outlier.colour = rgb(0.15, 0.3, 0.6, 0.35),
    
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    color = "darkred"
  ) +
  
  labs(
    title = "Residuals by strain",
    x = "Strains",
    y = "Residuals"
  ) +
  
  theme_diagnostics +
  
  
#Nur für die Stammnamen: 40° drehen und etwas kleiner darstellen
  theme(
    axis.text.x = element_text(
      family = "jost",
      angle = 40,
      hjust = 1,
      vjust = 1,
      size = 10,
      colour = "black"
    ),
    
    
#Abstand und Größe X-Achsentitel
    axis.title.x = element_text(
      family = "jost",
      size = 11,
      face = "plain",
      margin = margin(t = 8)
    ),

#Abstand und Größe Y-Achsentitel
   axis.title.y.left = element_text(
     family = "jost",
     size = 11,
     face = "plain",
     margin = margin(r = 10)
)

  )


diagnostic_combined <- (
  (p_random_effects | p_res_ploidy) /
    p_res_strain
) +
  
  plot_layout(
    heights = c(1, 1.2)
  ) +
  
  plot_annotation(
#    title = "Additional model diagnostics for nuclear volume analysis - Series 1",        #Zellkernvolumen, Nummer des Projekts ändern!
#    title = "Additional model diagnostics for histon volume analysis - Series 2",         #Histonvolumen, Nummer des Projekts ändern!
     title = "Additional model diagnostics for septin ring diameter analysis - Series 6",  #Septin Ring Durchmesser, Nummer des Projekts ändern!
     
    theme = theme(
      plot.title = element_text(
        family = "jost",
        size = 18,
        face = "plain",
        hjust = 0.5,
        margin = margin(
          t = 10,
          b = 8
        )
      )
    )
  )

print(diagnostic_combined)


#R-squared ausgeben (nur den Wert für das gewählte Modell ausgeben lassen)
print(r.squaredGLMM(model_raw_lin))
print(r.squaredGLMM(model_raw_fac))
print(r.squaredGLMM(model_log_lin))
print(r.squaredGLMM(model_log_fac))


#Estimated Marginal Means (emmeans) berechnen (Werte entsprechen dem geschätzten Mittelwert über alle Stämme hinweg)
ploidie_werte <- sort(unique(df_long$Ploidie))

#Estimated Marginal Mean für Raw_Lin Modell
emm_raw_lin <- emmeans(model_raw_lin, specs = "Ploidie", at = list(Ploidie = ploidie_werte))        
print(emm_raw_lin)

#Estimated Marginal Mean für Raw_Fac Modell
emm_raw_fac <- emmeans(model_raw_fac, specs = "Ploidie")        
print(emm_raw_fac)               

#Estimated Marginal Mean für Log_Lin Modell
emm_log_lin <- emmeans(model_log_lin, specs = "Ploidie", at = list(Ploidie = ploidie_werte))        
print(emm_log_lin)                                              
print(emm_log_lin, type = "response")

#Estimated Marginal Mean für Log_Fac Modell
emm_log_fac <- emmeans(model_log_fac, specs = "Ploidie")        
print(emm_log_fac)                                              
print(emm_log_fac, type = "response")                           


#Post-hoc paarweise Vergleiche (Tukey-korrigiert)

#Raw_Lin Modell
#pairs_raw_lin <- pairs(emm_raw_lin, adjust = "tukey", reverse = TRUE)         #Paarweise Vergleiche bei linearen Modell nicht sinnvoll (eine Steigung!)
#print(pairs_raw_lin)

#Raw_Fac Modell
pairs_raw_fac <- pairs(emm_raw_fac, adjust = "tukey", reverse = TRUE)
print(pairs_raw_fac)

#Log_Lin Modell     
#pairs_log_lin <- pairs(emm_log_lin, adjust = "tukey", reverse = TRUE)         #Paarweise Vergleiche bei linearen Modell nicht sinnvoll (eine Steigung!)
#print(pairs_log_lin)
#print(summary(pairs_log_lin, type = "response"))

#Log_Fac Modell
pairs_log_fac <- pairs(emm_log_fac, adjust = "tukey", reverse = TRUE)
print(pairs_log_fac)
print(summary(pairs_log_fac, type = "response"))


###################################################
#Plots (je nachdem welches Modell ausgewählt wurde)
###################################################

pred_raw_lin <- ggpredict(model_raw_lin, terms = "Ploidie [all]")                          #Falls Modell Raw Linear ausgewählt wurde
pred_raw_fac <- ggpredict(model_raw_fac, terms = "Ploidie")                                #Falls Modell Raw Factor ausgewählt wurde
pred_log_lin <- ggpredict(model_log_lin, terms = "Ploidie [all]", back_transform = TRUE)   #Falls Modell Log Linear ausgewählt wurde
pred_log_fac <- ggpredict(model_log_fac, terms = "Ploidie", back_transform = TRUE)         #Falls Modell Log Factor ausgewählt wurde  

ploidie_labels <- paste0(ploidie_werte, "N")


#Vorarbeit für alle Plots
#Median der einzelnen Stämme berechnen
strain_medians <- df_long %>%
  group_by(Stamm, Ploidie) %>%
  summarise(
    Median = median(Kernvolumen, na.rm = TRUE),
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
    n_strains = n(),
    
    x_position = ifelse(
      n_strains == 1,
      Ploidie,
      Ploidie +
        ((row_number() - 1) / (n_strains - 1) - 0.5) * 0.28
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
  10,       # Kreis mit Pluszeichen
  12,       # Quadrat mit Pluszeichen
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


#########################################
#Falls Modell Raw Linear ausgewählt wurde
#########################################
#Linear skaliert (Y = Volumen in µm³)

p1 <- ggplot() +
  
  geom_point(
    data = df_long,
    aes(
      x = Ploidie,
      y = Kernvolumen,
      color = factor(Ploidie),
      group = factor(Ploidie)
    ),
    shape = 1,
    alpha = 1.0,                                       #Transparenz der Datenpunkte verändern (1.0 vollständig deckend)
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
    
    size = 2.5,
    #    stroke = 2.0,
    #    colour = "black",    #Farbe der Symbole ändern
    #    fill = "black"       #Fill der Symbole ändern
  ) +
  
  
  #95%-Konfidenzintervall des linearen Modells
  geom_ribbon(
    data = pred_raw_lin,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    fill = "grey60",
    alpha = 0.20
  ) +
  
  #Geschätzter linearer Zusammenhang
  geom_line(
    data = pred_raw_lin,
    aes(
      x = x,
      y = predicted
    ),
    linewidth = 1.2,
    color = "black",
    lineend = "round"
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
  
  #Symbole der einzelnen Stämme  
  scale_shape_manual(
    values = strain_shapes
  ) +
  
  #X-Achse  
  scale_x_continuous(
    breaks = ploidie_werte,
    labels = ploidie_labels
  ) +
  
  # Y-Achse 
  scale_y_continuous(
    breaks = seq(0, 100, by = 5),           
    limits = c(0, 36),                     #Zweiten Wert verändern, je nachdem wie hoch die y-Achse sein muss
    expand = expansion(mult = c(0, 0.05))  
  ) +
  
  labs(
    x = "Ploidy level",
    y = bquote("Nuclear volume (" * mu * "m"^3 * ")")
#    y = bquote("Histone volume (" * mu * "m"^3 * ")")
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),                #Verhindert die Doppelung der X-Achse
    panel.background = element_blank(),          #Verhindert die Doppelung der X-Achse
    text = element_text(family = "jost"),        #Schriftart auf Jost ändern
    axis.title = element_text(family = "jost"),  #Schriftart auf Jost ändern
    axis.text = element_text(family = "jost"),   #Schriftart auf Jost ändern
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 20,
      margin = margin (t = 12)
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

print(p1)


#########################################
#Falls Modell Raw Factor ausgewählt wurde
#########################################
#Linear skaliert (Y = Volumen in µm³)

p2 <- ggplot() +
  
  geom_point(
    data = df_long,
    aes(
      x = Ploidie,
      y = Kernvolumen,
      color = factor(Ploidie),
      group = factor(Ploidie)
    ),
    shape = 1,
    alpha = 1.0,                            #Transparenz der Datenpunkte verändern (1.0 vollständig deckend)
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
    
    size = 2.5,
    #    stroke = 2.0,
    #    colour = "black",    #Farbe der Symbole ändern
    #    fill = "black"       #Fill der Symbole ändern
  ) +
  
  
  #Geschätztes 95% Konfidenzintervall  
  geom_errorbar(
    data = pred_raw_fac,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    
    width = 0.08,       #Größe ändern
    linewidth = 1.0,    #Größe ändern
    color = "black"
  ) +
  
  #Geschätzter Mittelwert des Modells als Punkt
  #  geom_point(
  #    data = pred_raw_fac,
  #    aes(
  #      x = x,
  #      y = predicted
  #    ),
  #    shape = 21,         #Shape = 21 Kreis, Shape = 24 Dreieck
  #    size = 3.0,         #Größe ändern
  #    stroke = 0.5,       #Größe ändern
  #    color = "black",
  #    fill = "black"
  #  ) +
  
  #Geschätzter Mittelwert des Modells als horizontale Linie
  geom_segment(
    data = pred_raw_fac,
    aes(
      x = x - 0.10,
      xend = x + 0.10,
      y = predicted,
      yend = predicted
    ),
    
    linewidth = 1.0,
    color = "black",
    lineend = "butt"
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
  scale_shape_manual(
    values = strain_shapes
  ) +
  
  #X-Achse  
  scale_x_continuous(
    breaks = ploidie_werte,
    labels = ploidie_labels
  ) +
  
  # Y-Achse 
  scale_y_continuous(
    breaks = seq(0, 100, by = 5),           
    limits = c(0, 36),                     #Zweiten Wert verändern, je nachdem wie hoch die y-Achse sein muss
    expand = expansion(mult = c(0, 0.05))  
  ) +
  
  labs(
    x = "Ploidy level",
    y = bquote("Nuclear volume (" * mu * "m"^3 * ")")
#    y = bquote("Histone volume (" * mu * "m"^3 * ")")
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),                #Verhindert die Doppelung der X-Achse
    panel.background = element_blank(),          #Verhindert die Doppelung der X-Achse
    text = element_text(family = "jost"),        #Schriftart auf Jost ändern
    axis.title = element_text(family = "jost"),  #Schriftart auf Jost ändern
    axis.text = element_text(family = "jost"),   #Schriftart auf Jost ändern
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 20,
      margin = margin (t = 12)
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

print(p2)


#########################################
#Falls Log-Linear-Modell ausgewählt wurde
#########################################
#Vorhersagen rücktransformiert, Darstellung auf ursprünglicher Volumenskala (µm³)

p3 <- ggplot() +
  
  geom_point(
    data = df_long,
    aes(
      x = Ploidie,
      y = Kernvolumen,
      color = factor(Ploidie),
      group = factor(Ploidie)
    ),
    shape = 1,
    alpha = 1.0,                                       #Transparenz der Datenpunkte verändern (1.0 vollständig deckend)
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
    
    size = 2.5,
    #    stroke = 2.0,
    #    colour = "black",    #Farbe der Symbole ändern
    #    fill = "black"       #Fill der Symbole ändern
  ) +
  
  
  #95%-Konfidenzintervall des Log-Linear-Modells
  geom_ribbon(
    data = pred_log_lin,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    fill = "grey60",
    alpha = 0.20
  ) +
  
  #Geschätzter Zusammenhang
  geom_line(
    data = pred_log_lin,
    aes(
      x = x,
      y = predicted
    ),
    linewidth = 1.2,
    color = "black",
    lineend = "round"
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
  scale_shape_manual(
    values = strain_shapes
  ) +
  
  #X-Achse  
  scale_x_continuous(
    breaks = ploidie_werte,
    labels = ploidie_labels
  ) +
  
  # Y-Achse 
  scale_y_continuous(
    breaks = seq(0, 100, by = 5),           
    limits = c(0, 36),                     #Zweiten Wert verändern, je nachdem wie hoch die y-Achse sein muss
    expand = expansion(mult = c(0, 0.05))  
  ) +
  
  labs(
    x = "Ploidy level",
    y = bquote("Nuclear volume (" * mu * "m"^3 * ")")
    #    y = bquote("Histone volume (" * mu * "m"^3 * ")")
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),                #Verhindert die Doppelung der X-Achse
    panel.background = element_blank(),          #Verhindert die Doppelung der X-Achse
    text = element_text(family = "jost"),        #Schriftart auf Jost ändern
    axis.title = element_text(family = "jost"),  #Schriftart auf Jost ändern
    axis.text = element_text(family = "jost"),   #Schriftart auf Jost ändern
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 20,
      margin = margin (t = 12)
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

print(p3)

#########################################
#Falls Modell Log Factor ausgewählt wurde
#########################################
#Vorhersagen rücktransformiert, Darstellung auf ursprünglicher Volumenskala (µm³)

p4 <- ggplot() +
  
  geom_point(
    data = df_long,
    aes(
      x = Ploidie,
      y = Kernvolumen,
      color = factor(Ploidie),
      group = factor(Ploidie)
    ),
    shape = 1,
    alpha = 1.0,                                       #Transparenz der Datenpunkte verändern (1.0 vollständig deckend)
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
    
    size = 2.5,
    #    stroke = 2.0,
    #    colour = "black",    #Farbe der Symbole ändern
    #    fill = "black"       #Fill der Symbole ändern
  ) +
  
  
  #Geschätztes 95% Konfidenzintervall  
  geom_errorbar(
    data = pred_log_fac,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    
    width = 0.08,       #Größe ändern
    linewidth = 1.0,    #Größe ändern
    color = "black"
  ) +
  
  #Rücktransformierter Modellschätzer als Punkt
  #  geom_point(
  #    data = pred_log_fac,
  #    aes(
  #      x = x,
  #      y = predicted
  #    ),
  #    shape = 21,         #Shape = 21 Kreis, Shape = 24 Dreieck
  #    size = 3.0,         #Größe ändern
  #    stroke = 0.5,       #Größe ändern
  #    color = "black",
  #    fill = "black"
  #  ) +
  
  #Rücktransformierter Modellschätzer als Punkt
  geom_segment(
    data = pred_log_fac,
    aes(
      x = x - 0.10,
      xend = x + 0.10,
      y = predicted,
      yend = predicted
    ),
    
    linewidth = 1.0,
    color = "black",
    lineend = "butt"
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
  scale_shape_manual(
    values = strain_shapes
  ) +
  
  #X-Achse  
  scale_x_continuous(
    breaks = ploidie_werte,
    labels = ploidie_labels
  ) +
  
  # Y-Achse für Histon und Zellkernvolumen
#  scale_y_continuous(
#    breaks = seq(0, 100, by = 5),           
#    limits = c(0, 18),                     #Zweiten Wert verändern, je nachdem wie hoch die y-Achse sein muss!
#    expand = expansion(mult = c(0, 0.05))  
#  ) +
  
  # Y-Achse für Septin-Ring Durchmesser
  scale_y_continuous(
    breaks = seq(0, 4.5, by = 0.5), 
    limits = c(0, 3.5),                    #Zweiten Wert verändern, je nachdem wie hoch die y-Achse sein muss!
    expand = expansion(mult = c(0,0.05))
  ) +
  
  labs(
    x = "Ploidy level",
#    y = bquote("Nuclear volume (" * mu * "m"^3 * ")")              #Zellkernvolumen! 
#    y = "Histone volume (µm³)"                                     #Histonvolumen!
     y = "Septin ring diameter (µm)"                                #Septin Ring Durchmesser!
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),                #Verhindert die Doppelung der X-Achse
    panel.background = element_blank(),          #Verhindert die Doppelung der X-Achse
    text = element_text(family = "jost"),        #Schriftart auf Jost ändern
    axis.title = element_text(family = "jost"),  #Schriftart auf Jost ändern
    axis.text = element_text(family = "jost"),   #Schriftart auf Jost ändern
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 20,
      margin = margin (t = 12)
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

print(p4)


##################################
#Kurze Zusammenfassung für Console
##################################

#Model Raw_Lin
print(summary(model_raw_lin)$coefficients)
print(r.squaredGLMM(model_raw_lin))

#Model Raw_Fac
print(summary(model_raw_fac)$coefficients)
print(r.squaredGLMM(model_raw_fac))

#Model Log_Lin
print(summary(model_log_lin)$coefficients)
slope_log <- fixef(model_log_lin)["Ploidie"]
cat(sprintf("  exp(%.4f) = %.3f -> +%.1f%% Volumen pro Ploidie-Einheit\n",slope_log, exp(slope_log), (exp(slope_log) - 1) * 100))
print(r.squaredGLMM(model_log_lin))

#Model Log_Fac
print(summary(model_log_fac)$coefficients)
print(r.squaredGLMM(model_log_fac))


#######################################################################################################################################

##########################################################
#Säulendiagramm Fold-Changes ohne beobachtete Fold-Changes
##########################################################

#Fold-Change-Kontraste gegenüber 1N  
fc_log_fac <- contrast(                #An Modellauswahl anpassen!
  emm_log_fac,                            #An Modellauswahl anpassen!
  method = "trt.vs.ctrl",                 
  ref = 1                                 
)


# Vorhandene Ploidiestufen außer Referenz 1N bestimmen
comparison_ploidies <- sort(
  setdiff(unique(df_long$Ploidie), 1)
)


#Fold-Changes, 95%-Konfidenzintervalle und p-Werte
fc_table <- summary(
  fc_log_fac,                          #An Modellauswahl anpassen!
  type = "response",
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

print(fc_table)
names(fc_table)


#Konfidenzintervall-Spalten vereinheitlichen
if ("asymp.LCL" %in% names(fc_table)) {
  names(fc_table)[names(fc_table) == "asymp.LCL"] <- "lower.CL"
}

if ("asymp.UCL" %in% names(fc_table)) {
  names(fc_table)[names(fc_table) == "asymp.UCL"] <- "upper.CL"
}


#Vorhandene Ploidiestufen außer Referenz 1N bestimmen
comparison_ploidies <- sort(
  setdiff(unique(df_long$Ploidie), 1)
)


#Tabelle für den Plot vorbereiten
fc_table <- fc_table %>%
  mutate(
    Ploidie = comparison_ploidies,
    Fold_Change = ratio,
    Modell = "Log-factor LMM"
  ) %>%
  select(
    Ploidie,
    Fold_Change,
    lower.CL,
    upper.CL,
    p_value = p.value,
    Modell
  )


#Referenzgruppe 1N ergänzen
fc_fac_complete <- data.frame(              #An Modellauswahl anpassen!
  Ploidie = 1,
  Fold_Change = 1,
  lower.CL = 1,
  upper.CL = 1,
  p_value = NA_real_,
  Modell = "Log-factor LMM"
) %>%
  bind_rows(fc_table)


# P-Werte hinzufügen
fc_final <- fc_fac_complete %>%            #An Modellauswahl anpassen!
  mutate(
    Ploidie_Label = paste0(Ploidie, "N")
  )


# Signifikanzsterne
fc_final <- fc_final %>%
  mutate(
    significance = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(fc_final)


# Anzahl der Vergleiche
n_comparisons <- length(comparison_ploidies)


# Höhe der Signifikanzklammern
# Höchste obere CI-Grenze bestimmen
max_ci <- max(
  fc_final$upper.CL,
  na.rm = TRUE
)

# Signifikanzklammern oberhalb der Konfidenzintervalle
bracket_y <- seq(
  from = max_ci + 0.3,
  by = 0.3,
  length.out = n_comparisons
)

y_max <- max(bracket_y) + 0.25

# Tabelle für Signifikanzklammern
brackets <- data.frame(
  x_start = rep(1, n_comparisons),
  x_end = comparison_ploidies,
  y = bracket_y,
  significance = fc_final %>%
    filter(Ploidie != 1) %>%
    arrange(Ploidie) %>%
    pull(significance)
)

print(brackets)

#Position der Fold-Change-Beschriftungen
fc_final <- fc_final %>%
  mutate(
    label_y = ifelse(
      Ploidie == 1,
      Fold_Change + 0.08,
      upper.CL + 0.08
    )
  )


#Plot
p_fc <- ggplot() +
  
# Säulen
  geom_col(
    data = fc_final,
    aes(
      x = Ploidie,
      y = Fold_Change,
      fill = Ploidie_Label
    ),
    
# Rahmenlinie um jede Säule
    color = "black",
    width = 0.6,
    linewidth = 0.6
  ) +
  
#95%-Konfidenzintervalle der modellgeschätzten Fold-Changes
  geom_errorbar(
    data = fc_final %>%
      filter(Ploidie != 1),
    aes(
      x = Ploidie,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.10,
    linewidth = 0.8,
    color = "black"
  ) +
  
  
  
#Fold-Change-Wert innerhalb der Säule
  geom_text(
    data = fc_final,
    aes(
      x = Ploidie,
      y = 0.30,                                     #Position der Fold Changes verändern
      label = sprintf("%.2f", Fold_Change)
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
      yend = y - 0.08
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
      yend = y - 0.08
    ),
    color = "black",
    linewidth = 0.6
  ) +
  
  
#Signifikanzangaben über den Klammern
  geom_text(
    data = brackets,
    aes(
      x = (x_start + x_end) / 2,
      y = y + 0.05,
      label = significance
    ),
    family = "jost",
    size = 5.0,
    color = "black",
    vjust = 0
  ) +
  
  
#X-Achse
  scale_x_continuous(
    breaks = fc_final$Ploidie,
    labels = fc_final$Ploidie_Label
  ) +
  
  
#Y-Achse
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(y_max),
      by = 0.5
    ),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  coord_cartesian(
    ylim = c(0, y_max)
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
#    y = "Fold change in nuclear volume"                #Zellkernvolumen!
#    y = "Fold change in histone volume"                 #Histonvolumen!
    y = "Fold change in septin ring diameter"          #Septin Ring Durchmesser!
  ) +
  
  
  theme_bw() +
  
  theme(
    # Schriftart
    text = element_text(
      family = "jost"
    ),
    
    axis.title = element_text(
      family = "jost"
    ),
    
    axis.text = element_text(
      family = "jost"
    ),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    
    # X-Achsenwerte
    axis.text.x = element_text(
      size = 19,
      colour = "black"
    ),
    
    # Y-Achsenwerte
    axis.text.y = element_text(
      size = 19,
      colour = "black"
    ),
    
    # X-Achsentitel
    axis.title.x = element_text(
      size = 22,
      margin = margin(t = 12),
      face = "plain"
    ),
    
    # Y-Achsentitel
    axis.title.y = element_text(
      size = 22,
      margin = margin(r = 10),
      face = "plain"
    ),
    
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    
    legend.position = "none"
  )

print(p_fc)
print(fc_final)

##############################################################################################################################


