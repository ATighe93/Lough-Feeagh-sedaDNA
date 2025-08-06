library(reshape2)
library(zoo)
library(ggplot2)
library(cowplot)

###### species #####
scale_3Gen <- c("#FF8C00" #darkorange
								, "#9400D3" #darkviolet
								, "#FF1493" #deeppink
								)

scale_3Gen <- c("#ff8c00"
                , "#d94600"
                , "#b20000")

species <- "Canis_lupus"
#   species <- "Rangifer_tarandus_platyrhynchus"
#   species <- "Salmo_salar"
depth1 <- "S109" # -55 cm, 390 YBP 
depth2 <- "S144" # -405 cm, 4000 YBP
depth3 <- "S95" # -625cm, 9420 YBP

###### LENGTH ######
ext_length <- "_lgdistribution.txt"

file.species_depth1 <- Sys.glob(paste(depth1, "*", species, ext_length, sep = ""))
file.species_depth2 <- Sys.glob(paste(depth2, "*", species, ext_length, sep = ""))
file.species_depth3 <- Sys.glob(paste(depth3, "*", species, ext_length, sep = ""))

sample_depth1 <- read.table(file.species_depth1[1], header=T)
sample_depth2 <- read.table(file.species_depth2[1], header=T)
sample_depth3 <- read.table(file.species_depth3[1], header=T)
sample_depth1_ok <- aggregate(sample_depth1$Occurences, by=list(sample_depth1$Length), FUN=sum)
sample_depth2_ok <- aggregate(sample_depth2$Occurences, by=list(sample_depth2$Length), FUN=sum)
sample_depth3_ok <- aggregate(sample_depth3$Occurences, by=list(sample_depth3$Length), FUN=sum)
colnames(sample_depth1_ok) <- c("Length", "-55 cm")
colnames(sample_depth2_ok) <- c("Length", "-405 cm")
colnames(sample_depth3_ok) <- c("Length", "-625 cm")

sample_Total_1 <- merge(sample_depth1_ok, sample_depth2_ok, by = "Length", all = T)
sample_Total_2 <- merge(sample_Total_1, sample_depth3_ok, by = "Length", all = T)
colnames(sample_Total_2) = c("Length", "390 YBP", "4,000 YBP", "9,420 YBP")

sample_Total_2[is.na(sample_Total_2)] <- 0
sample_Total_3 <- melt(sample_Total_2, id.vars = "Length")
sample_Total_4 <- unlist(sapply(1:nrow(sample_Total_3), function(x) rep(sample_Total_3$Length[x], e = sample_Total_3$value[x])))
sample_Total_5 <- unlist(sapply(1:nrow(sample_Total_3), function(x) rep(sample_Total_3$variable[x], e = sample_Total_3$value[x])))
sample_Total <- data.frame(length = sample_Total_4, Samples = sample_Total_5)

sample_length_depth1 <- sample_Total[sample_Total$Samples=="390 YBP",]
sample_length_depth1$length <- as.numeric(sample_length_depth1$length)
sample_length_depth2 <- sample_Total[sample_Total$Samples=="4,000 YBP",]
sample_length_depth2$length <- as.numeric(sample_length_depth2$length)
sample_length_depth3 <- sample_Total[sample_Total$Samples=="9,420 YBP",]
sample_length_depth3$length <- as.numeric(sample_length_depth3$length)

plot_length <- ggplot(sample_Total, aes(x=length, fill=Samples, after_stat(scaled))) + 
  geom_density(aes(group=Samples, colour=Samples), alpha = 0.05) +
  xlab("Length (bp)") +
  ylab("Frequency") +
  xlim(0,300) +
  scale_color_manual(values = scale_3Gen) +
  scale_fill_manual(values = scale_3Gen)

ggsave(
  "Canis_lupus_length_distribution.pdf",
  plot = plot_length,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)
ggsave(
  "Canis_lupus_length_distribution.tiff",
  plot = plot_length,
  device = "tiff",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)

##### SUBSTITUTION #####
ext_3 <- "_misincorporation.txt"
ext_5 <- "_misincorporation.txt"

L25 <- 1:25

#### 3' end ####
file.mis_depth1_3p <- Sys.glob(paste(depth1, "*", species, ext_3, sep = ""))
file.mis_depth2_3p <- Sys.glob(paste(depth2, "*", species, ext_3, sep = ""))
file.mis_depth3_3p <- Sys.glob(paste(depth3, "*", species, ext_3, sep = ""))

sample_depth1_3p <- read.table(file.mis_depth1_3p[1], skip = 0, header=T)
sample_depth1_3p$Depth <- "Depth -55cm"
sample_depth1_3p$Std <- as.character(sample_depth1_3p$Std)
sample_depth1_3p$Pos <- as.integer(sample_depth1_3p$Pos)
sample_depth1_3p_25 <- sample_depth1_3p[sample_depth1_3p$Pos %in% L25 & sample_depth1_3p$End=="3p",]
sample_depth1_3p_25_ok <- sample_depth1_3p_25[, 5:30]
sample_depth1_3p_25_ok2 <- aggregate.data.frame(sample_depth1_3p_25_ok, sample_depth1_3p_25_ok['Pos'] , FUN=sum)
sample_depth1_3p_25_fin <- as.data.frame(sample_depth1_3p[L25, 4])
sample_depth1_3p_25_fin$'3pGtoA' <- sample_depth1_3p_25_ok2$G.A/sample_depth1_3p_25_ok2$G
sample_depth1_3p_25_fin$`sample_depth1_3p[L25, 4]` <- NULL
sample_depth1_3p_25_fin$Position <- 1:25

sample_depth2_3p <- read.table(file.mis_depth2_3p[1], skip = 0, header=T)
sample_depth2_3p$Depth <- "Depth -405cm"
sample_depth2_3p$Std <- as.character(sample_depth2_3p$Std)
sample_depth2_3p$Pos <- as.integer(sample_depth2_3p$Pos)
sample_depth2_3p_25 <- sample_depth2_3p[sample_depth2_3p$Pos %in% L25 & sample_depth2_3p$End=="3p",]
sample_depth2_3p_25_ok <- sample_depth2_3p_25[, 5:30]
sample_depth2_3p_25_ok2 <- aggregate.data.frame(sample_depth2_3p_25_ok, sample_depth2_3p_25_ok['Pos'] , FUN=sum)
sample_depth2_3p_25_fin <- as.data.frame(sample_depth2_3p[L25, 4])
sample_depth2_3p_25_fin$'3pGtoA' <- sample_depth2_3p_25_ok2$G.A/sample_depth2_3p_25_ok2$G
sample_depth2_3p_25_fin$`sample_depth2_3p[L25, 4]` <- NULL
sample_depth2_3p_25_fin$Position <- 1:25

sample_depth3_3p <- read.table(file.mis_depth3_3p[1], skip = 0, header=T)
sample_depth3_3p$Depth <- "Depth -625cm"
sample_depth3_3p$Std <- as.character(sample_depth3_3p$Std)
sample_depth3_3p$Pos <- as.integer(sample_depth3_3p$Pos)
sample_depth3_3p_25 <- sample_depth3_3p[sample_depth3_3p$Pos %in% L25 & sample_depth3_3p$End=="3p",]
sample_depth3_3p_25_ok <- sample_depth3_3p_25[, 5:30]
sample_depth3_3p_25_ok2 <- aggregate.data.frame(sample_depth3_3p_25_ok, sample_depth3_3p_25_ok['Pos'] , FUN=sum)
sample_depth3_3p_25_fin <- as.data.frame(sample_depth3_3p[L25, 4])
sample_depth3_3p_25_fin$'3pGtoA' <- sample_depth3_3p_25_ok2$G.A/sample_depth3_3p_25_ok2$G
sample_depth3_3p_25_fin$`sample_depth3_3p[L25, 4]` <- NULL
sample_depth3_3p_25_fin$Position <- 1:25

sample_3p_inter1 <- merge(sample_depth1_3p_25_fin, sample_depth2_3p_25_fin, by = "Position", all=T)
sample_3p_inter2 <- merge(sample_3p_inter1, sample_depth3_3p_25_fin, by = "Position", all=T)
colnames(sample_3p_inter2) <- c("Position", "-55cm", "-405cm", "-625cm")
sample_3p_fin <- melt(sample_3p_inter2, id.vars = "Position")
sample_3p_fin$Type <- "3'end G/A"
sample_3p_fin$Main_mutation_type <- "Others"
sample_3p_fin$Main_mutation_type[sample_3p_fin$Depth=="-55cm"] <- "G/A"
sample_3p_fin$Main_mutation_type[sample_3p_fin$Depth=="-55cm"] <- "Others"
colnames(sample_3p_fin) <- c("Position", "Depth", "Frequency", "Type", "Main_mutation_type")
sample_3p_fin$Position_trick <- -sample_3p_fin$Position

plot_3p <- ggplot(sample_3p_fin, aes(x = Position_trick, y = Frequency*100, colour = Depth, linetype = Main_mutation_type)) +
  geom_line(aes(color = Depth)) +
  scale_color_manual(values = scale_3Gen) +
  scale_linetype_manual(values = c("Others"="dotted", "G/A"="solid"))+
  xlim(-15, 0) +
  ylim(0, 5) +
  labs(x = "Position in read from 3'end (nt)", y = "3' end G to A deamination percentage", color = "Depth", linetype = "Main mutation type on the 3 
       most external positions")

ggsave(
  "Canis_lupus_3end_deamination_pattern.pdf",
  plot = plot_3p,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)
ggsave(
  "Rangifer_tarandus_platyrhynchus_3end_deamination_pattern.pdf",
  plot = plot_3p,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)
ggsave(
  "Salmo_salar_3end_deamination_pattern.pdf",
  plot = plot_3p,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)

#### 5' End ####
file.mis_depth1_5p <- Sys.glob(paste(depth1, "*", species, ext_5, sep = ""))
file.mis_depth2_5p <- Sys.glob(paste(depth2, "*", species, ext_5, sep = ""))
file.mis_depth3_5p <- Sys.glob(paste(depth3, "*", species, ext_5, sep = ""))


sample_depth1_5p <- read.table(file.mis_depth1_5p[1], skip = 0, header=T)
sample_depth1_5p$Depth <- "Depth -55 cm"
sample_depth1_5p$Std <- as.character(sample_depth1_5p$Std)
sample_depth1_5p$Pos <- as.integer(sample_depth1_5p$Pos)
sample_depth1_5p_25 <- sample_depth1_5p[sample_depth1_5p$Pos %in% L25 & sample_depth1_5p$End=="5p",]
sample_depth1_5p_25_ok <- sample_depth1_5p_25[, 5:30]
sample_depth1_5p_25_ok2 <- aggregate.data.frame(sample_depth1_5p_25_ok, sample_depth1_5p_25_ok['Pos'] , FUN=sum)
sample_depth1_5p_25_fin <- as.data.frame(sample_depth1_5p[L25, 4])
sample_depth1_5p_25_fin$'5pCtoT' <- sample_depth1_5p_25_ok2$C.T/sample_depth1_5p_25_ok2$C
sample_depth1_5p_25_fin$`sample_depth1_5p[L25, 4]` <- NULL
sample_depth1_5p_25_fin$Position <- 1:25

sample_depth2_5p <- read.table(file.mis_depth2_5p[1], skip = 0, header=T)
sample_depth2_5p$Depth <- "Depth -405 cm"
sample_depth2_5p$Std <- as.character(sample_depth2_5p$Std)
sample_depth2_5p$Pos <- as.integer(sample_depth2_5p$Pos)
sample_depth2_5p_25 <- sample_depth2_5p[sample_depth2_5p$Pos %in% L25 & sample_depth2_5p$End=="5p",]
sample_depth2_5p_25_ok <- sample_depth2_5p_25[, 5:30]
sample_depth2_5p_25_ok2 <- aggregate.data.frame(sample_depth2_5p_25_ok, sample_depth2_5p_25_ok['Pos'] , FUN=sum)
sample_depth2_5p_25_fin <- as.data.frame(sample_depth2_5p[L25, 4])
sample_depth2_5p_25_fin$'5pCtoT' <- sample_depth2_5p_25_ok2$C.T/sample_depth2_5p_25_ok2$C
sample_depth2_5p_25_fin$`sample_depth2_5p[L25, 4]` <- NULL
sample_depth2_5p_25_fin$Position <- 1:25

sample_depth3_5p <- read.table(file.mis_depth3_5p[1], skip = 0, header=T)
sample_depth3_5p$Depth <- "Depth -625 cm"
sample_depth3_5p$Std <- as.character(sample_depth3_5p$Std)
sample_depth3_5p$Pos <- as.integer(sample_depth3_5p$Pos)
sample_depth3_5p_25 <- sample_depth3_5p[sample_depth3_5p$Pos %in% L25 & sample_depth3_5p$End=="5p",]
sample_depth3_5p_25_ok <- sample_depth3_5p_25[, 5:30]
sample_depth3_5p_25_ok2 <- aggregate.data.frame(sample_depth3_5p_25_ok, sample_depth3_5p_25_ok['Pos'] , FUN=sum)
sample_depth3_5p_25_fin <- as.data.frame(sample_depth3_5p[L25, 4])
sample_depth3_5p_25_fin$'5pCtoT' <- sample_depth3_5p_25_ok2$C.T/sample_depth3_5p_25_ok2$C
sample_depth3_5p_25_fin$`sample_depth3_5p[L25, 4]` <- NULL
sample_depth3_5p_25_fin$Position <- 1:25

sample_5p_inter1 <- merge(sample_depth1_5p_25_fin, sample_depth2_5p_25_fin, by = "Position", all=T)
sample_5p_inter2 <- merge(sample_5p_inter1, sample_depth3_5p_25_fin, by = "Position", all=T)
colnames(sample_5p_inter2) <- c("Position", "390 YBP", "4,000 YBP", "9,420 YBP")
sample_5p_fin <- melt(sample_5p_inter2, id.vars = "Position")
sample_5p_fin$Type <- "5'end C/T"
sample_5p_fin$Main_mutation_type <- "Others"
colnames(sample_5p_fin) <- c("Position", "Samples", "Frequency", "Type", "Main_mutation_type")
sample_5p_fin$Main_mutation_type[sample_5p_fin$Samples=="390 YBP"] <- "C>T"

plot_5p <- ggplot(sample_5p_fin, aes(x = Position, y = Frequency*100, colour = Samples, linetype = Main_mutation_type)) +
  geom_line(aes(color = Samples)) +
  scale_color_manual(values = scale_3Gen) +
  scale_linetype_manual(values = c("Others"="dotted", "C>T"="solid"))+
  xlim(0, 16) +
  ylim(0, 5) +
  labs(x = "Position in read from 5' end (bp)", y = "5' end C to T deamination percentage", color = "Samples", linetype = "Main mutation type on the \n3 most external positions") +
  theme(legend.justification=1, legend.text.align = 0)

ggsave(
  "Canis_lupus_5end_deamination_pattern_v2.pdf",
  plot = plot_5p,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)
ggsave(
  "Canis_lupus_5end_deamination_pattern_v2.tiff",
  plot = plot_5p,
  device = "tiff",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)

ggsave(
  "Rangifer_tarandus_platyrhynchus_5end_deamination_pattern.pdf",
  plot = plot_5p,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)
ggsave(
  "Salmo_salar_5end_deamination_pattern.pdf",
  plot = plot_5p,
  device = "pdf",
  scale = 1,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
)


#### ONE PLOT WITH 3 GRAPHS ####
plots <- ggdraw() +
  draw_plot(plot = plot_length, x = -0.007, y = 0.2, width = 0.5, height = 0.65) +
  draw_plot(plot = plot_5p, x = 0.305, y = 0.2, width = 0.45, height = 0.65) +
  draw_plot(plot = plot_3p, x = 0.565, y = 0.2, width = 0.45, height = 0.65)

ggsave(filename <- paste(species
                         ,"_damage_patterns.tiff"
                         , sep = "")
       , plot = plots
       , device = "tiff"
       , scale = 1
       , width = 21
       , height = 11.5
       , units = "cm"
       , dpi = 300
       , limitsize = TRUE
)

ggdraw() +
  draw_plot(plot = plot_5p, x = 0.000, y = 0.2, width = 0.45, height = 0.65) +
  draw_plot(plot = plot_3p, x = 0.300, y = 0.2, width = 0.45, height = 0.65)

#### ONE PLOT WITH 2 GRAPHS ####
plots_min <- ggdraw() +
  draw_plot(plot = plot_length, x = -0.007, y = 0.2, width = 0.5, height = 0.65) +
  draw_plot(plot = plot_5p, x = 0.35, y = 0.2, width = 0.62, height = 0.65) 
ggsave(filename <- paste(species
                            ,"_damage_patterns_v2.pdf"
                            , sep = "")
       , plot = plots_min
       , device = "pdf"
       , scale = 1
       , width = 21
       , height = 11.5
       , units = "cm"
       , dpi = 300
       , limitsize = TRUE
)
ggsave(filename <- paste(species
                         ,"_damage_patterns_v2.tiff"
                         , sep = "")
       , plot = plots_min
       , device = "tiff"
       , scale = 1
       , width = 21
       , height = 11.5
       , units = "cm"
       , dpi = 300
       , limitsize = TRUE
)
