## clear environment
rm(list=ls()) 


# install.packages("MetaCycle")
library(MetaCycle)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("patchwork")
library(patchwork)


# seed random number generator for reproducibility
set.seed(20251012)

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

# make simulated expression data


num_subjects <- 20
amp <- 2
# acrophase <- 9
acrophase <- 13
mesor <- 5

samples_times <-c(0, 6, 12, 18)
times <- rep(samples_times, num_subjects)
sample_count <- length(times)

ids <- c()
for (num in 1:num_subjects){
  ids <- append(ids, rep(num,length(samples_times)))
}



expression <- amp*cos((times-acrophase)*(2*pi/24)) + mesor
expression <- expression + rnorm(length(expression), mean = 0, sd = 0.15)




###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
# cosinor 


trig_model <- lm(expression ~ cos(times*pi/12) + sin(times*pi/12) ) 

cos_fit <- unname(trig_model$coefficients[2])
sin_fit <- unname(trig_model$coefficients[3])
# mesor_fit <- unname(model$coefficients[1])

amp_fit <- unname(sqrt(cos_fit**2 + sin_fit**2))
phase_rad_fit <- atan2(sin_fit,cos_fit)
phase_hour_fit <- unname(phase_rad_fit * (24 / (2 * pi)) )

if(phase_hour_fit < 0){
  phase_hour_fit <- 24 + phase_hour_fit
}



###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
# metacycle


expression_matrix <- matrix(expression, nrow = 1)
rownames(expression_matrix) <- "Gene1"
colnames(expression_matrix) <- paste0("Sample", seq_along(expression))


if (!dir.exists("input")) {
  dir.create("input")
}

input_path <- paste0("input/amp_test_acrophase",acrophase,"/")


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



output_path <- paste0("output/amp_test_acrophase",acrophase,"/")


if (!dir.exists(output_path)) {
  dir.create(output_path)
}



meta3d(datafile           = paste0(input_path,"datafile.txt"),
       designfile         = paste0(input_path,"designfile.txt"),
       outdir             = output_path,
       filestyle          = "txt",
       design_libColm     = 1,
       design_subjectColm = 2,
       design_hrColm      = 3,
       minper             = 24,
       maxper             = 24,
       cycMethodOne       = "ARS",
       timeUnit           = "hour")


# metacycle predicted phases
results <- read.table(paste0(output_path,"/meta3dGroupID_AllSubjects_designfile.txt"), header = TRUE, sep = "\t")
# meta_phases <- results$meta3d_Phase


meta_phase <- results$meta3d_Phase[1]
meta_amp   <- results$meta3d_AMP[1]

cat("MetaCycle Phase:", meta_phase, "hours\n")
cat("MetaCycle Amplitude:", meta_amp, "\n")

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
# timeseries plot





duplicate_timepoints <- function(data, original_times, shift = 24) {
  dup_data <- subset(data, times %in% original_times)
  
  dup_data$times <- dup_data$times + shift
  
  return(dup_data)
}



adjusted_expression <- expression - mean(expression)

y_min <- min(adjusted_expression)
y_max <- max(adjusted_expression)



fit_x <- seq(0, 48, by = 1)
fit_y <- cos_fit * cos(fit_x * (2 * pi / 24)) + sin_fit * sin(fit_x * (2 * pi / 24))





data <- data.frame(times = times, adjusted_expression = adjusted_expression)

original_times <- c(0, 6, 12, 18)
dup_data <- duplicate_timepoints(data, original_times, shift = 24)

fit_data <- data.frame(
  fit_x = seq(0, 48, by = 1),
  fit_y = cos_fit * cos(fit_x * (2 * pi / 24)) + sin_fit * sin(fit_x * (2 * pi / 24)) 
  
)

fit_data_solid <- subset(fit_data, fit_x <= 18)
fit_data_dotted <- subset(fit_data, fit_x > 18)

p1 <- ggplot() +
  geom_point(data = data, aes(x = times, y = adjusted_expression), color = "black", size = 2.2) +
  
  geom_point(data = dup_data, aes(x = times, y = adjusted_expression), color = "gray", alpha = 0.9, size = 2.2) +
  
  geom_line(data = fit_data_solid, aes(x = fit_x, y = fit_y), color = "blue", linetype = "solid",linewidth = 1.1) +
  
  geom_line(data = fit_data_dotted, aes(x = fit_x, y = fit_y), color = "blue", linetype = "dotted",linewidth = 1.1) +
  
  
  scale_x_continuous(limits = c(0, 48), 
                     breaks = seq(0, 48, by = 6)) +
  
  labs(x = "Time (hours)", y = "MESOR-Corrected Expression") +
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
    
  )










fit_meta <- data.frame(
  fit_x = fit_x,
  fit_y = meta_amp * cos( 2*pi * (fit_x - meta_phase) / 24 )
)
fit_meta_solid <- subset(fit_meta, fit_x <= 18)
fit_meta_dotted<- subset(fit_meta, fit_x >  18)





p2 <- ggplot() +
  geom_point(data = data, aes(x = times, y = adjusted_expression), color = "black", size = 2.2) +
  
  geom_point(data = dup_data, aes(x = times, y = adjusted_expression), color = "gray", alpha = 0.9, size = 2.2) +
  
  # geom_line(data = fit_data_solid, aes(x = fit_x, y = fit_y), color = "blue", linetype = "solid") +
  # 
  # geom_line(data = fit_data_dotted, aes(x = fit_x, y = fit_y), color = "blue", linetype = "dotted") +
  
  # MetaCycle fits (red)
  geom_line(data = fit_meta_solid, aes(x=fit_x, y=fit_y),
            inherit.aes = FALSE, color="blue",   linetype="solid",linewidth = 1.1) +
  geom_line(data = fit_meta_dotted, aes(x=fit_x, y=fit_y),
            inherit.aes = FALSE, color="blue",   linetype="dotted",linewidth = 1.1) +
  
  
  scale_x_continuous(limits = c(0, 48), 
                     breaks = seq(0, 48, by = 6)) +
  
  labs(x = "Time (hours)", y = "MESOR-Corrected Expression") +
  
  theme_minimal() 



cosinor_title <- paste0("Cosinor\nTrue Acrophase: ",acrophase,", Fit Acrophase: ",round(phase_hour_fit,2),
                        "\nTrue Amplitude: ",amp,", Fit Amplitude: ",round(amp_fit,2))

metacycle_title <- paste0("MetaCycle-ARS\nTrue Acrophase: ",acrophase,", Fit Acrophase: ",round(meta_phase,2),
                        "\nTrue Amplitude: ",amp,", Fit Amplitude: ",round(meta_amp,2))


p1 <- p1 + 
  ggtitle(cosinor_title) +
  scale_y_continuous(limits = c(y_min - 0.5, y_max + 0.5))

p2 <- p2 + 
  ggtitle(metacycle_title) +
  scale_y_continuous(limits = c(y_min - 0.5, y_max + 0.5))

# Combine the plots side by side
combined_plot <- p1 + p2 + 
  plot_layout(ncol = 2) & 
  theme(axis.title.y = element_text(vjust = 1)) & 
  theme(
    axis.title.x = element_text(size = 16,  margin = margin(t = 6)),
    axis.title.y = element_text(size = 16,  margin = margin(r = 6)),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 16)
  )

combined_plot <- combined_plot + 
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )




plot_path <- paste0("plots/amp_test_acrophase",acrophase,".png")


ggsave(plot_path,
       combined_plot,
       width = 10, height = 6, units = "in", dpi = 300,
       bg = "white")


cat("Cosinor Phase:", phase_hour_fit, "hours\n")
cat("Cosinor Amplitude:", amp_fit, "\n")

cat("MetaCycle Phase:", meta_phase, "hours\n")
cat("MetaCycle Amplitude:", meta_amp, "\n")

