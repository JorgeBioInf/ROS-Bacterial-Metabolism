"""
Computes distribution curves of the following metrics:
  - Mean minimum distance between residues susceptible to ROS and Active Sites
  - Mean minimum distance between residues susceptible to ROS and Binding Sites
for each protein within the experimentally evidenced data and the rest of the
proteome
"""


### DISTRIBUTION PLOT 
rm(list = ls()) 
library(dplyr)


setwd("C:/Users/JorgeXD/Documents/TFM/Curacion_dataset")
file <- readr::read_tsv('ROS_summary_full.tsv')
evi_file <- readr::read_tsv('Pputida_evidenced.txt')
no_evi <- anti_join(
  file,
  evi_file,
  by = "KEGG ID"
)


# 0) CALCULATING PARAMETERS
# Active sites
ASR = no_evi[["# SR in AS"]]
pos_ASR = ASR[ASR > 0]

ASR_perc = length(pos_ASR) / length (ASR)
ASR_perc

filt_ASR = no_evi[["# SR in AS"]]
pos_filt_ASR = filt_ASR[filt_ASR > 0]

filt_ASR_perc = length(pos_filt_ASR) / length(filt_ASR)
filt_ASR_perc


# Binding sites
BSR = no_evi[["# SR in BS"]]
pos_BSR = BSR[BSR > 0]

BSR_perc = length(pos_BSR) / length (BSR)
BSR_perc

filt_BSR = no_evi[["# SR in BS"]]
pos_filt_BSR = filt_BSR[filt_BSR > 0]

filt_BSR_perc = length(pos_filt_BSR) / length(filt_BSR)
filt_BSR_perc


# Both
ASR_BSR <- no_evi[, c("# SR in AS", "# SR in BS")]
pos_ASR_BSR <- ASR_BSR[ASR_BSR[["# SR in AS"]] > 0 & ASR_BSR[["# SR in BS"]] > 0,] 

ASR_BSR_perc = dim(pos_ASR_BSR) / dim(ASR_BSR) * 100
ASR_BSR_perc

filt_ASR_BSR <- evi_file[, c("# SR in AS", "# SR in BS")]
filt_pos_ASR_BSR <- filt_ASR_BSR[filt_ASR_BSR[["# SR in AS"]] > 0 & filt_ASR_BSR[["# SR in BS"]] > 0,]

filt_ASR_BSR_perc = dim(filt_pos_ASR_BSR) / dim(filt_ASR_BSR) * 100
filt_ASR_BSR_perc


## PLOTS

png("ROS_residues_plots.png", width = 1200, height = 800)
par(mfrow = c(2, 2))


# 1) ACTIVE SITES
# Residue count
AS_dens_all <- density(no_evi[["# SR in AS"]], na.rm = TRUE)
AS_dens_filtered <- density(evi_file[["# SR in AS"]], na.rm = TRUE)

plot(AS_dens_all,
     main = "Active sites",
     xlab = "Susceptible residues count",
     col = "darkblue",
     lwd = 2)

lines(AS_dens_filtered, col = "darkgreen", lwd = 2)

polygon(AS_dens_all, col = rgb(0.1, 0, 1, 0.3), border = NA)
polygon(AS_dens_filtered, col = rgb(0.1, 1, 0, 0.3), border = NA)

legend("topright",
       legend = c("No evidence", "With evidence"),
       col = c("darkblue", "darkgreen"),
       lwd = 2,
       bty = "n")


# Susceptible residue distances
no_null_all <- no_evi[no_evi["Mean SR - AS minimum distance"] != 0,]
no_null_filt <- evi_file[evi_file["Mean SR - AS minimum distance"] != 0,]

ASD_dens_all <- density(no_null_all[["Mean SR - AS minimum distance"]], na.rm = TRUE)
ASD_dens_filtered <- density(no_null_filt[["Mean SR - AS minimum distance"]], na.rm = TRUE)

plot(ASD_dens_all,
     main = "Active sites",
     xlab = "Average susceptible residue's distance",
     col = "darkblue",
     lwd = 2)

lines(ASD_dens_filtered, col = "darkgreen", lwd = 2)

polygon(ASD_dens_all, col = rgb(0.1, 0, 1, 0.3), border = NA)
polygon(ASD_dens_filtered, col = rgb(0.1, 1, 0, 0.3), border = NA)

legend("topright",
       legend = c("No evidence", "With evidence"),
       col = c("darkblue", "darkgreen"),
       lwd = 2,
       bty = "n")


# 2) BINDING SITES
# Susceptible residue counts
BS_dens_all <- density(no_evi[["# SR in BS"]], na.rm = TRUE)
BS_dens_filtered <- density(evi_file[["# SR in BS"]], na.rm = TRUE)

plot(BS_dens_all,
     main = "Binding sites",
     xlab = "Susceptible residues count",
     col = "darkblue",
     lwd = 2)

lines(BS_dens_filtered, col = "darkgreen", lwd = 2)

polygon(BS_dens_all, col = rgb(0.1, 0, 1, 0.3), border = NA)
polygon(BS_dens_filtered, col = rgb(0.1, 1, 0, 0.3), border = NA)

legend("topright",
       legend = c("No evidence", "With evidence"),
       col = c("darkblue", "darkgreen"),
       lwd = 2,
       bty = "n")


# Susceptible residue distances
no_null_all_BS <- no_evi[no_evi["Mean SR - BS minimum distance"] != 0,]
no_null_filt_BS <- evi_file[evi_file["Mean SR - BS minimum distance"] != 0,]

BSD_dens_all <- density(no_null_all_BS[["Mean SR - BS minimum distance"]], na.rm = TRUE)
BSD_dens_filtered <- density(no_null_filt_BS[["Mean SR - BS minimum distance"]], na.rm = TRUE)


plot(BSD_dens_all,
     main = "Binding sites",
     xlab = "Average susceptible residue's distance",
     col = "darkblue",
     lwd = 2)

lines(BSD_dens_filtered, col = "darkgreen", lwd = 2)

polygon(BSD_dens_all, col = rgb(0.1, 0, 1, 0.3), border = NA)
polygon(BSD_dens_filtered, col = rgb(0.1, 1, 0, 0.3), border = NA)

legend("topright",
       legend = c("No evidence", "With evidence"),
       col = c("darkblue", "darkgreen"),
       lwd = 2,
       bty = "n")

dev.off()
