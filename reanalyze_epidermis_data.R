## clear environment
rm(list=ls())



########################################################################################
# GSE139301 EPIDERMIS DATASET

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139301
# https://pubmed.ncbi.nlm.nih.gov/30377266/

## OUTLIER TO REMOVE: 115 (only 3 samples)

########################################################################################


## probe_to_gene mapping
platform_info <- read.table('epidermis_data/GPL13667-15572_processed.csv', header = TRUE, sep =',')
platform_info$Gene.Symbol[platform_info$Gene.Symbol==""] <- platform_info$ID[platform_info$Gene.Symbol==""]
platform_info$Gene.Symbol[platform_info$Gene.Symbol=="---"] <- platform_info$ID[platform_info$Gene.Symbol=="---"]
probe_gene_labels <- unlist(as.vector(platform_info["Gene.Symbol"]))

data <- read.table('epidermis_data/GSE139301_processed.txt', header = TRUE, sep ='\t',row.names=1)

colnames <- colnames(data)

ids <- c()
times <- c()

for(sample in colnames){
  tmp <- unlist(strsplit(sample,".",fixed=T))

  id <- unlist(strsplit(tmp[1],"_",fixed=T))[2]

  time <- as.numeric(unlist(strsplit(tmp[2],"_",fixed=T))[2])

  ids <- append(ids,id)
  times <- append(times,time)

}


times <- times/100



times <- times[ids != "115"]
data <- data[,ids != "115"]
ids <- ids[ids != "115"]
id_factors <- factor(ids)




######################################################
# filtering low expression (mean < 4) probes

data_means <- unname(rowMeans(data))
data <- data[data_means >= 4, ]
probe_gene_labels <- probe_gene_labels[data_means >= 4]

######################################################

data <- 2^data

unique_genes <- unlist(levels(factor(probe_gene_labels)))

gene_data <- matrix(NA, nrow = length(unique_genes), ncol = dim(data)[2])

for (gene_ind in 1:length(unique_genes)){

  gene_to_find <- unique_genes[gene_ind]

  sample_averages <- unname(colMeans(data[(probe_gene_labels==gene_to_find),]))

  gene_data[gene_ind, ] <- sample_averages

}

rownames(gene_data) <- unique_genes
colnames(gene_data) <- colnames(data)





########################################################################################
# COSINOR ANALYSIS




p_values <- c()
acrophases <- c()


for(gene_ind in 1:dim(gene_data)[1]){
  
  
  expression <- unlist(as.vector(data[gene_ind,]))
  
  trig_model <- lm(expression ~ cos(times*pi/12) + sin(times*pi/12) + id_factors) 
  base_model <- lm(expression ~ id_factors) 
  
  
  anova_results <- anova(trig_model,base_model)
  p_value <- anova_results[2,6]
  
  p_values <- append(p_values, p_value)
  
  
  cos_coef <- trig_model$coefficients[2]
  sin_coef <- trig_model$coefficients[3]
  


  
  # phase fit
  phase_rad <- atan2(sin_coef,cos_coef)
  
  ## convert radians to hours on 24 hour clock
  phase_hour <- unname(phase_rad * (24 / (2 * pi)) )
  
  if(phase_hour < 0){
    phase_hour <- 24 + phase_hour
  }
  
  acrophases <- append(acrophases, phase_hour)
  
}




# select cycling genes
bh_q <- p.adjust(p_values,method="BH")

sig_genes <- rownames(gene_data)[bh_q <= 0.05]
sig_acrophases <- acrophases[bh_q<=0.05]




num_bins <- 48

df <- data.frame(
  phase = sig_acrophases
)

rose_plot <- ggplot(df, aes(x = phase)) +
  geom_histogram(
    binwidth = 24 / num_bins,
    boundary = 0,
    closed   = "left",
    position = "identity",
    alpha    = 0.5,
    fill     = "#D55E00",
    color    = "#D55E00" # outline
    
  ) +
  coord_polar(theta = "x", start = 4*pi/2) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = 0:23,
    labels = paste0(0:23, "h"),
    expand = c(0, 0)
  ) +
  labs(title = "Cosinor Acrophases", x = "Acrophase (hours)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5)
  )




if (!dir.exists("plots")) {
  dir.create("plots")
}


plot_path <- paste0("plots/epidermis_cosinor.png")
ggsave(plot_path, rose_plot, width = 8, height = 6, dpi = 300)




















stop()




########################################################################################
# METACYCLE-ARSER ANALYSIS





if (!dir.exists("input")) {
  dir.create("input")
}



if (!dir.exists("input/epidermis")) {
  dir.create("input/epidermis")
}



write.table(gene_data,
            file = "input/epidermis/datafile.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA)



design <- data.frame(
  libraryID = colnames(gene_data),
  subjectID = ids,
  hour = times,
  stringsAsFactors = FALSE
)


write.table(design,
            file = "input/epidermis/designfile.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)





if (!dir.exists("output")) {
  dir.create("output")
}

if (!dir.exists("output/epidermis")) {
  dir.create("output/epidermis")
}



meta3d(datafile           = "input/epidermis/datafile.txt",
       designfile         = "input/epidermis/designfile.txt",
       outdir             = "output/epidermis",
       filestyle          = "txt",
       design_libColm     = 1,
       design_subjectColm = 2,
       design_hrColm      = 3,
       minper             = 24,
       maxper             = 24,
       cycMethodOne       = "ARS",
       timeUnit           = "hour")




metacycle_results <- read.table("output/epidermis/meta3dGroupID_AllSubjects_designfile.txt", header = TRUE, sep = "\t")



metacycle_results_to_plot <- metacycle_results[metacycle_results$CycID %in% sig_genes, ]



num_bins <- 48

df <- data.frame(
  phase = metacycle_results_to_plot$meta3d_Phase
)

rose_plot <- ggplot(df, aes(x = phase)) +
  geom_histogram(
    binwidth = 24 / num_bins,
    boundary = 0,
    closed   = "left",
    position = "identity",
    alpha    = 0.5,
    fill     = "#D55E00",
    color    = "#D55E00" # outline
    
  ) +
  coord_polar(theta = "x", start = 4*pi/2) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = 0:23,
    labels = paste0(0:23, "h"),
    expand = c(0, 0)
  ) +
  labs(title = "MetaCycle-ARS Acrophases", x = "Acrophase (hours)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5)
  )




if (!dir.exists("plots")) {
  dir.create("plots")
}


plot_path <- paste0("plots/epidermis_metacycle.png")
ggsave(plot_path, rose_plot, width = 8, height = 6, dpi = 300)


