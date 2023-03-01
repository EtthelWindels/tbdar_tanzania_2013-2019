## ------------------------------------------------------------------------
## Create posterior density plots for birth-death analyses
## 2022-05-19 Etthel Windels
## ------------------------------------------------------------------------


# Load libraries ----------------------------------------------------------

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(MCMCvis)))
suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(dplyr)))


# Parser ------------------------------------------------------------------

parser <- argparse::ArgumentParser()
parser$add_argument("--logfile_1", type="character", help="Log file L1")
parser$add_argument("--logfile_2", type="character", help="Log file L2")
parser$add_argument("--logfile_3", type="character", help="Log file L3")
parser$add_argument("--logfile_4", type="character", help="Log file L4")
parser$add_argument("--output_path", type="character", help="Path to output figures")
parser$add_argument("--Re_out", type="character", help="Outplot figure Re")
parser$add_argument("--lambda_out", type="character", help="Outplot figure lambda")
parser$add_argument("--infperiod_out", type="character", help="Outplot figure infectious period")
parser$add_argument("--Re_prior_out", type="character", help="Outplot figure Re + prior")
parser$add_argument("--bUR_prior_out", type="character", help="Outplot figure bUR + prior")

args <- parser$parse_args()


# Read arguments ----------------------------------------------------------

LOGFILE_1 <- args$logfile_1
LOGFILE_2 <- args$logfile_2
LOGFILE_3 <- args$logfile_3
LOGFILE_4 <- args$logfile_4
OUTPUT_PATH <- args$output_path
RE_OUT <- args$Re_out
LAMBDA_OUT <- args$lambda_out
INFPERIOD_OUT <- args$infperiod_out
RE_PRIOR_OUT <- args$Re_prior_out
bUR_PRIOR_OUT <- args$bUR_prior_out

print(paste("logfile_1: ", LOGFILE_1))
print(paste("logfile_2: ", LOGFILE_2))
print(paste("logfile_3: ", LOGFILE_3))
print(paste("logfile_4: ", LOGFILE_4))
print(paste("output_path: ", OUTPUT_PATH))
print(paste("Re_out: ", RE_OUT))
print(paste("infperiod_out: ", INFPERIOD_OUT))
print(paste("lambda_out: ", LAMBDA_OUT))
print(paste("Re_prior_out: ", RE_PRIOR_OUT))
print(paste("bUR_prior_out: ", bUR_PRIOR_OUT))


# Function to load log files ----------------------------------------------

loadLog <- function(filename, burninFrac=0.1, subsample=NA) {
  
  df_in <- as.matrix(read.table(filename, header=T))
  
  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }
  
  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }
  
  return(df_in)
}


# Load log files ----------------------------------------------------------

logfile_L1 <- loadTrajectories(LOGFILE_1, burninFrac=0) %>%
              as.data.frame()
logfile_L1$lambda <- logfile_L1$reproductiveNumber*logfile_L1$becomeUninfectiousRate

logfile_L2 <- loadTrajectories(LOGFILE_2, burninFrac=0) %>%
  as.data.frame()
logfile_L2$lambda <- logfile_L2$reproductiveNumber*logfile_L2$becomeUninfectiousRate

logfile_L3 <- loadTrajectories(LOGFILE_3, burninFrac=0) %>%
  as.data.frame()
logfile_L3$lambda <- logfile_L3$reproductiveNumber*logfile_L3$becomeUninfectiousRate

logfile_L4 <- loadTrajectories(LOGFILE_4, burninFrac=0) %>%
  as.data.frame()
logfile_L4$lambda <- logfile_L4$reproductiveNumber*logfile_L4$becomeUninfectiousRate


log_data <- data.frame(lineage=c(rep("L1",dim(logfile_L1)[1]), rep("L2",dim(logfile_L2)[1]), rep("L3",dim(logfile_L3)[1]), rep("L4",dim(logfile_L4)[1])),
                      Re=c(logfile_L1$reproductiveNumber, logfile_L2$reproductiveNumber, logfile_L3$reproductiveNumber, logfile_L4$reproductiveNumber),
                      lambda=c(logfile_L1$lambda, logfile_L2$lambda, logfile_L3$lambda, logfile_L4$lambda),
                      bUR=c(logfile_L1$becomeUninfectiousRate, logfile_L2$becomeUninfectiousRate, logfile_L3$becomeUninfectiousRate, logfile_L4$becomeUninfectiousRate),
                      inf_period=c(1/logfile_L1$becomeUninfectiousRate, 1/logfile_L2$becomeUninfectiousRate, 1/logfile_L3$becomeUninfectiousRate, 1/logfile_L4$becomeUninfectiousRate))


# Define prior distributions ----------------------------------------------

prior_R0 <- rlnorm(100000,0,1)    
prior_bUR <- rlnorm(100000,0,0.5)


prior_data <- data.frame(lineage=c(rep("L1",dim(logfile_L1)[1]), rep("L2",dim(logfile_L2)[1]), rep("L3",dim(logfile_L3)[1]), rep("L4",dim(logfile_L4)[1]), rep('prior', length(prior_R0))),
                    R0=c(logfile_L1$reproductiveNumber, logfile_L2$reproductiveNumber, logfile_L3$reproductiveNumber, logfile_L4$reproductiveNumber, prior_R0),
                    bUR=c(logfile_L1$becomeUninfectiousRate, logfile_L2$becomeUninfectiousRate, logfile_L3$becomeUninfectiousRate, logfile_L4$becomeUninfectiousRate, prior_bUR))
prior_data$lineage <- as.factor(prior_data$lineage)
prior_data$lineage <- factor(prior_data$lineage, levels=c('prior','L1','L2','L3','L4'))


# Plot posterior densities ------------------------------------------------

Re <- ggplot(log_data, aes(Re,fill=lineage, colour=lineage)) +
  geom_density(stat='density',alpha=0.5)+
  labs(title=expression(atop(bold("Reproductive number"))), x=expression(R[e]), y='Posterior density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0.5,1.5)+
  scale_fill_manual(labels=c('L1','L2','L3','L4'),values=c('#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))+
  scale_colour_manual(labels=c('L1','L2','L3','L4'),values=c('#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))

lambda <- ggplot(log_data, aes(lambda,fill=lineage, colour=lineage)) +
  geom_density(stat='density',alpha=0.5)+
  labs(title=expression(atop(bold("Transmission rate"))), x='Transmission rate', y='Posterior density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,1.5)+
  scale_fill_manual(labels=c('L1','L2','L3','L4'),values=c('#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))+
  scale_colour_manual(labels=c('L1','L2','L3','L4'),values=c('#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))

inf_period <- ggplot(log_data, aes(inf_period,fill=lineage, colour=lineage)) +
  geom_density(stat='density',alpha=0.5)+
  labs(title=expression(atop(bold("Infectious period"))), x='Infectious period (years)', y='Posterior density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,7.5)+
  scale_fill_manual(labels=c('L1','L2','L3','L4'),values=c('#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))+
  scale_colour_manual(labels=c('L1','L2','L3','L4'),values=c('#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))

ggsave(RE_OUT, plot=Re, path=OUTPUT_PATH, width=159, height=178, units="mm")
ggsave(LAMBDA_OUT, plot=lambda, path=OUTPUT_PATH, width=159, height=178, units="mm")
ggsave(INFPERIOD_OUT, plot=inf_period, path=OUTPUT_PATH, width=159, height=178, units="mm")

 
# Plot prior + posterior densities ----------------------------------------

R0prior <- ggplot(prior, aes(R0,fill=lineage, colour=lineage)) +
  geom_density(stat='density',alpha=0.5)+
  labs( x=expression(R[e]~susceptible), y='Density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,2.0)+
  scale_fill_manual(labels=c('prior','posterior_L1','posterior_L2','posterior_L3','posterior_L4'),values=c('lightgrey','#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))+
  scale_colour_manual(labels=c('prior','posterior_L1','posterior_L2','posterior_L3','posterior_L4'),values=c('lightgrey','#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))

deltaprior <- ggplot(prior, aes(bUR,fill=lineage, colour=lineage)) +
  geom_density(stat='density',alpha=0.5)+
  labs( x="Becoming uninfectious rate", y='Density')+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust=0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=25, hjust=0.5)) +
  xlim(0,3)+
  scale_fill_manual(labels=c('prior','posterior_L1','posterior_L2','posterior_L3','posterior_L4'),values=c('lightgrey','#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))+
  scale_colour_manual(labels=c('prior','posterior_L1','posterior_L2','posterior_L3','posterior_L4'),values=c('lightgrey','#f4b79fff','#8a96a3ff','#bfcbdbff','#c08168ff'))


ggsave(RE_PRIOR_OUT, plot=R0prior, path=OUTPUT_PATH, width=159, height=178, units="mm")
ggsave(bUR_PRIOR_OUT, plot=deltaprior, path=OUTPUT_PATH, width=159, height=178, units="mm")
