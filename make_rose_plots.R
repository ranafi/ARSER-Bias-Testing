
## clear environment
rm(list=ls()) 





# install.packages("MetaCycle")
library(MetaCycle)
# install.packages("ggplot2")
library(ggplot2)




make_rose_plot <- function(amp = 2, mesor = 5, num_genes = 3000, num_subjects = 20, sample_times = c(0, 6, 12, 18), example_num = 0) {
  
  times <- rep(sample_times, num_subjects)
  sample_count <- length(times)
  
  ids <- c()
  for (num in 1:num_subjects){
    ids <- append(ids, rep(num,length(sample_times)))
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
  
  input_path <- paste0("input/example",example_num,"/")
  
  

  if (!dir.exists(input_path)) {
    dir.create(input_path)
  }
  
  
  write.table(expression_matrix,
              file = paste0(input_path,"datafile.txt"),
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
              file = paste0(input_path,"designfile.txt"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  
  
  
  
  
  
  if (!dir.exists("output")) {
    dir.create("output")
  }
  
  
  
  
  if (!dir.exists(paste0("output/example",example_num))) {
    dir.create(paste0("output/example",example_num))
  }
  
  
  meta3d(datafile           = paste0(input_path,"datafile.txt"),
         designfile         = paste0(input_path,"designfile.txt"),
         outdir             = paste0("output/example",example_num),
         filestyle          = "txt",
         design_libColm     = 1,
         design_subjectColm = 2,
         design_hrColm      = 3,
         minper             = 24,
         maxper             = 24,
         cycMethodOne       = "ARS",
         timeUnit           = "hour")
  
  
  # metacycle predicted phases
  results <- read.table(paste0("output/example",example_num,"/meta3dGroupID_AllSubjects_designfile.txt"), header = TRUE, sep = "\t")
  meta_phases <- results$meta3d_Phase
  
  
  
  
  
  
  
  data_label1 <- "True Acrophases"
  data_label2 <- "MetaCycle-ARS Estimates"
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
      "True Acrophases"          = "#0072B2",
      "MetaCycle-ARS Estimates" = "#D55E00"
    )) +
    coord_polar(theta = "x", start = 4*pi/2) +
    scale_x_continuous(
      limits = c(0, 24),
      breaks = 0:23,
      labels = paste0(0:23, "h"),
      expand = c(0, 0)
    ) +
    labs(title = NULL, x = "Acrophase (hours)", y = NULL, fill = NULL) +
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
  
  
  plot_path <- paste0("plots/rose_plot",example_num,".png")
  ggsave(plot_path, rose_plot, width = 8, height = 6, dpi = 300)
  
}







# seed random number generator for reproducibility
set.seed(20251013)


# FIG 3A - adding 2 to each sample time, modes also shift by 2
make_rose_plot(example_num = 1, sample_times = c(2, 8, 14, 20))

# FIG 3B - increased sampling frequency doesn't really help
make_rose_plot(example_num = 2, sample_times = seq(0,22,2))

# FIG 3C - more than one full period DOES help
make_rose_plot(example_num = 3, sample_times = c(0, 6, 12, 18, 24, 30, 36, 42))

