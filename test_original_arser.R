
## clear environment
rm(list=ls()) 


# install.packages("ggplot2")
library(ggplot2)

# seed random number generator for reproducibility
set.seed(20251021)


num_genes <- 3000

num_subjects <- 20
amp <- 2
mesor <- 5

samples_times <-seq(0,22,2)


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



#####################################################################################
# make inputs for python arser

if (!dir.exists("input")) {
  dir.create("input")
}

if (!dir.exists("input/python_arser")) {
  dir.create("input/python_arser")
}


if (!dir.exists("output")) {
  dir.create("output")
}

if (!dir.exists("output/python_arser")) {
  dir.create("output/python_arser")
}




design <- data.frame(
  libraryID = colnames(expression_matrix),
  subjectID = ids,
  hour = times,
  stringsAsFactors = FALSE
)



time_grid <- sort(unique(design$hour))
collapsed_mat <- sapply(time_grid, function(h) {
  cols <- which(design$hour == h)
  if (length(cols) == 1L) expression_matrix[, cols] else rowMeans(expression_matrix[, cols, drop = FALSE])
})

colnames(collapsed_mat) <- as.character(time_grid)  # e.g., "0" "6" "12" "18"

arser_df <- data.frame(Probes = rownames(expression_matrix),
                       collapsed_mat,
                       check.names = FALSE)

out_path <- "input/python_arser/python_arser_input.txt"
write.table(arser_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)


###################################################################################################################################

# at this point you need to run the original ARSER code on python_arser_input.txt
# this is a bit difficult because it's written in python2 with several old dependencies.
# especially tough if you have a new mac with an apple chip.
# only way i could get it running was in a docker container.

# first, clone/download original ARSER from here
# https://github.com/cauyrd/ARSER


# then, make two folders called "input" and "output" in the ARSER directory, and put python_arser_input.txt in the "input" folder

# then, create a DOCKERFILE with this text inside:
######################################################
# # file: Dockerfile
# # (Keep using --platform=linux/amd64 at build/run time instead of pinning it here.)
# 
# FROM ubuntu:18.04
# 
# # system deps + R + headers for building rpy2
# RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
# build-essential curl git ca-certificates \
# python2.7 python2.7-dev python-pip \
# r-base=3.4.4-1ubuntu1 r-base-dev=3.4.4-1ubuntu1 \
# libreadline-dev libbz2-dev liblzma-dev libpcre3-dev \
# libcurl4-openssl-dev libxml2-dev zlib1g-dev libssl-dev \
# && rm -rf /var/lib/apt/lists/*
#   
#   # make sure R is discoverable
#   ENV R_HOME=/usr/lib/R
# ENV PATH="/usr/lib/R/bin:${PATH}"
# 
# # py2-compatible toolchain + libs
# RUN python2.7 -m pip install --upgrade 'pip<21' 'setuptools<45' 'wheel<1' \
# && python2.7 -m pip install \
# 'numpy<1.17' 'scipy<1.3' 'matplotlib<3.1' \
# 'rpy2==2.8.6'
# 
# WORKDIR /opt/ARSER
# CMD ["python2.7", "arser.py"]
######################################################


# then build the docker image this with line in the terminal:
# docker build --no-cache --platform=linux/amd64 -t arser-legacy .


# then run ARSER on our input data with this line in the terminal:
# docker run --rm -it --platform=linux/amd64 \
# -e MPLBACKEND=Agg \
# -v "$PWD":/opt/ARSER -w /opt/ARSER arser-legacy \
# python2.7 arser.py input/python_arser_input.txt output/python_arser_output.txt 24 24 24


# and after it runs, the results will be in output/python_arser_output.txt
# just grab that file and put it back in output/python_arser/ and then continue with this R code


###################################################################################################################################

arser_results <- read.table("output/python_arser/python_arser_output.txt", header = TRUE, sep = "\t")

arser_phases <- as.numeric(arser_results$phase)


data_label1 <- "True Acrohases"
data_label2 <- "ARSER Estimates"
num_bins    <- 48

df <- data.frame(
  phase   = c(true_acrophases, arser_phases),
  dataset = c(rep(data_label1, length(true_acrophases)),
              rep(data_label2,  length(arser_phases)))
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
    "ARSER Estimates" = "#D55E00"
  )) +
  coord_polar(theta = "x", start = 4*pi/2) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = 0:23,
    labels = paste0(0:23, "h"),
    expand = c(0, 0)
  ) +
  labs(title = "ARSER Phase Estimates", x = NULL, y = NULL, fill = "") +
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

ggsave("plots/python_arser.png", rose_plot, width = 8, height = 6, dpi = 300)

