
## clear environment
rm(list=ls()) 


# install.packages("MetaCycle")
library(MetaCycle)
# install.packages("ggplot2")
library(ggplot2)

# seed random number generator for reproducibility
set.seed(20251012)


num_genes <- 3000
num_subjects <- 20
amp <- 2
mesor <- 5
samples_times <-c(0, 6, 12, 18)
times <- rep(samples_times, num_subjects)
sample_count <- length(times)

ids <- c()
for (num in 1:num_subjects){
  ids <- append(ids, rep(num,length(samples_times)))
}


true_acrophases <- c()
expression_matrix <- matrix(NA, nrow = num_genes, ncol = sample_count)

for (gene_ind in 1:num_genes){
  
  acrophase <- runif(1, min = 0, max = 24)
  expression <- amp*cos((times-acrophase)*(2*pi/24)) + mesor
  expression <- expression + rnorm(length(expression), mean = 0, sd = 0.15)
  expression_matrix[gene_ind,] <- expression
  true_acrophases <- append(true_acrophases, acrophase)
  
}

rownames(expression_matrix) <- paste0("Gene", seq_len(num_genes))
colnames(expression_matrix) <- paste0("Sample", seq_len(sample_count))





if (!dir.exists("input")) {
  dir.create("input")
}


if (!dir.exists("input/simple_example")) {
  dir.create("input/simple_example")
}




write.table(expression_matrix,
            file = "input/simple_example/datafile.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA)


design <- data.frame(
  libraryID = colnames(expression_matrix),
  subjectID = ids,
  hour = times,
  stringsAsFactors = FALSE
)


write.table(design,
            file = "input/simple_example/designfile.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



if (!dir.exists("output")) {
  dir.create("output")
}


if (!dir.exists("output/simple_example")) {
  dir.create("output/simple_example")
}



meta3d(datafile           = "input/simple_example/datafile.txt",
       designfile         = "input/simple_example/designfile.txt",
       outdir             = "output/simple_example",
       filestyle          = "txt",
       design_libColm     = 1,
       design_subjectColm = 2,
       design_hrColm      = 3,
       minper             = 24,
       maxper             = 24,
       cycMethodOne       = "ARS",
       timeUnit           = "hour")


# metacycle predicted phases
results <- read.table("output/simple_example/meta3dGroupID_AllSubjects_designfile.txt", header = TRUE, sep = "\t")
meta_phases <- results$meta3d_Phase







data_label1 <- "True Acrohases"
data_label2 <- "Metacycle-ARS Estimates"
num_bins    <- 48

df <- data.frame(
  phase   = c(true_acrophases, meta_phases),
  dataset = c(rep(data_label1, length(true_acrophases)),
              rep(data_label2,  length(meta_phases)))
)

rose_plot <- ggplot(df, aes(x = phase, fill = dataset)) +
  geom_histogram(
    binwidth = 24 / num_bins,
    boundary = 0,
    closed   = "left",
    position = "identity",
    alpha    = 0.5
  ) +
  scale_fill_manual(values = c(
    "True Acrohases"          = "#0072B2",
    "Metacycle-ARS Estimates" = "#D55E00"
  )) +
  coord_polar(theta = "x", start = 4*pi/2) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = 0:23,
    labels = paste0(0:23, "h"),
    expand = c(0, 0)
  ) +
  labs(title = "Metacycle-ARS Phase Estimates", x = NULL, y = NULL, fill = "") +
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
    plot.title = element_text(size = 14, hjust = 0.5)
  )





if (!dir.exists("plots")) {
  dir.create("plots")
}



ggsave("plots/simple_example.png", rose_plot, width = 8, height = 6, dpi = 300)



############################################################################################################

# predicted amp percent error plot

amp_percent_errors <- (results$meta3d_AMP - amp) / amp * 100


df_plot <- data.frame(
  true_acrophases = true_acrophases,
  amp_percent_errors = amp_percent_errors
)

scatter_plot <- ggplot(df_plot, aes(x = true_acrophases, y = amp_percent_errors)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(x = "True Acrophase (hours)", y = "Predicted Amplitude % Error ") + 
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24))



filename <- paste("plots/simple_example_predicted_amps.png",sep="")
ggsave(filename, plot = scatter_plot, width = 4, height = 4, dpi = 300)














