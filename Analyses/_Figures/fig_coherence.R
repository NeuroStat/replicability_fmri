####################
#### TITLE:     Plot index of coherence under increasing sample size.
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/Analyses/_Figures
#### First Modified: 31/05/2018
#### Notes: 
#################



##
###############
### Notes
###############
##

# Data analysis:
# Make independent batches of subjects.
# Each run, we estimate the mixture of two binomials on the summed thresholded maps within a sample size.

# Here we read in the processed data, calculate the coherence using Cohen's Kappa and plot the results.
# We do this for each sample size and each run.



##
###############
### Preparation
###############
##

# Location of intermediate results
LocIntRes <- '../_IntData/'

# Possible contrasts: default = MATH > LANGUAGE
contrast <- c('ML', 'Faces','Incentive', 'StopGo')

# Save locations depending on contrast
contrSave <- c('1_cognitive', '2_faces', '3_incentive', '4_stopgo')

# Load in libraries
library(tidyverse)
library(magrittr)
library(cowplot)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(reshape2)
library(RColorBrewer)
library(segmented)
library(Hmisc)
library(NeuRRoStat)

# Function described by Thirion et al. (2007) to calculate Cohen's Kappa
ThirionKappa <- function(lambda, piA1, piI1){
  # Check if the values are within 0 and 1
  if(any(c(lambda, piA1, piI1) < 0 | c(lambda, piA1, piI1) > 1)) error('Values should be between 0 and 1!')
  
  # Now calculate the piA0 and piI0
  piA0 <- 1 - piA1
  piI0 <- 1 - piI1
  
  # Calculate the proportion of correct classifications (p0)
  p0 <- lambda*piA1 + (1 - lambda) * piI0
  
  ## Agreement by chance
  # pi1 or also denoted as tau
  pi0 <- lambda * piA0 + (1 - lambda) * piI0
  # pi1 or also given by 1 - tau
  pi1 <- 1 - pi0
  # pC
  pC <- lambda * pi0 + (1 - lambda) * pi1
  
  # Cohen's Kappa, defined as Thirion
  CohKap_Thi <- (p0 - pC)/(1 - pC)
  
  return(CohKap_Thi)
}

# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs
NRUNS <- 50

# Starting amount of subjects in group analysis
startingTotal <- 10

# Database total
DATTOTAL <- 1400

# Total amount of possible subjects if we have a maximum of 3 disjoint subsets
NTOT <- 460

# Steps in the sequence (sample sizes)
steps <- seq(startingTotal, NTOT, by = 10)

# Create the sequence
sequence <- data.frame()
for(k in 1:NRUNS){
  for(i in 1:length(steps)){
    # Number of disjoint sets in this sample size
    numDS <- floor(DATTOTAL/steps[i])
    for(s in 1:numDS){
      sequence <- rbind(sequence,
                        data.frame('run' = k, 'step' = i, 'set' = s))
    }
  }
}

# Dimension of this data frame
dimSeq <- dim(sequence)[1]

# Variables for plotting
subjBreak <- c(seq(10,110,by=20), seq(150,NTOT, by=50))

##
###############
### Read in data
###############
##

# Data frame with estimated parameters using EM: voxelwise FDR and uncorrected
EMParam <- data.frame() %>% as_tibble()
for(s in 1:length(contrast)){
  # Select the contrast
  contr <- contrast[s]
  
  # Select the possible thresholding scenario's
  if(contr == 'ML'){
    scenario_pos <- c('unc', 'fdr')
  }
  if(contr != 'ML'){
    scenario_pos <- c('fdr')
  }

  # Loop over the thresholding scenario's and bind to data frame
  for(r in 1:length(scenario_pos)){
    EMParam <-
      readRDS(paste0(LocIntRes, contr, '/', 'EM_params_', scenario_pos[r],'.rda')) %>% 
        mutate(contrast = contr,
          threshold = scenario_pos[r]) %>%
      bind_rows(EMParam, .)
  }
}

# Add sequence information: don't know what this is anymore...
EMParam$sequence <- factor(rep(rep(1:(dim(EMParam)[1]/4), each = 2), 2))

##
###############
### Check results of EM: FDR and MATH contrast only
###############
##

# Let us plot the amount of iterations needed to estimate the parameters
EMParam %>% filter(final == TRUE & contrast == 'ML') %>%
  filter(threshold == 'fdr') %>%
  ggplot(., aes(sequence, num.iter)) +
  geom_col()
# Not very informative...

# Better have some kind of barplot
EMParam %>% filter(threshold == 'fdr' & contrast == 'ML') %>%
  ggplot(., aes(num.iter, group = final)) +
    geom_bar(aes(fill = final), position = 'dodge2') +
    scale_y_continuous("Count over all steps * runs") +
    scale_x_continuous('Number of iterations') +
    ggtitle('Number of iterations needed in the two-step EM algorithm',
            subtitle = 'Split according to first run or final estimation.') +
    theme_classic()
ggsave(filename = paste(getwd(), '/1_cognitive/EM_Checks/iterations_EM.png', sep = ''), plot = last_plot())

# Difference in coefficients
# Values of first run
valF <- EMParam %>% filter(final == FALSE & contrast == 'ML') %>% 
  filter(threshold == 'fdr') %>%
  select(lambda, PI1, PI2)
# Values of second run
valS <- EMParam %>% filter(final == TRUE & contrast == 'ML') %>% 
  filter(threshold == 'fdr') %>%
  select(lambda, PI1, PI2)

# Plot difference
plotRuns <- seq(1,NRUNS,by = 6)
for(r in 1:length(plotRuns)){
  #print(plotRuns[r] : (plotRuns[r] + 5))
  ToPrint <- plotRuns[r] : (plotRuns[r] + 5)
  if(max(ToPrint) > 50){
    ToPrint <- min(ToPrint):50
  }
  ToPlot <- data.frame(valF - valS, step = EMParam$step, run = EMParam$run) %>% 
    filter(run %in% ToPrint) %>%
    gather(., key = 'parameter', value = 'value', 1:3) %>%
    ggplot(., aes(x = step, y = value, group = parameter)) +
    geom_line(aes(colour = parameter), size = 0.7) +
    facet_wrap( ~ run) +
    ggtitle(paste0('Differences in parameter estimates between first and second attempt
                   in the two step EM algorithm'))
  ggsave(filename = paste(getwd(), '/1_cognitive/EM_Checks/DiffEM_part', r, '.png', sep = ''), plot = ToPlot)
}

# Seems ok! Go with final estimates.

##
###############
### Calculate Cohen's Kappa: FDR and MATH contrast
###############
##

# Filter the final estimates
Kappa <- EMParam %>% filter(final == TRUE & contrast == 'ML') %>%
  filter(threshold == 'fdr') %>%
  # Remove columns that don't provide information
  select(-num.iter, -final, -sequence) %>%
  # Transform step to sample size
  mutate(SampleSize = step * 10) %>%
  # Calculate Kappa
  mutate(kappa = NeuRRoStat::CohenKappa(lambda = lambda, piA1 = PI1, piI1 = PI2)) %>%
  # Put NaN to 0
  mutate(kappa = ifelse(is.nan(kappa), 0, kappa)) %>%
  # Remove columns
  select(run, SampleSize, kappa)

# Average, min and max per sample size
Kappa %>% group_by(SampleSize) %>%
  summarise(AvgKap = mean(kappa),
            minKap = min(kappa),
            maxKap = max(kappa)) 

# Min and max median
Kappa %>% group_by(SampleSize) %>%
  summarise(MedKap = median(kappa)) %>%
  filter(MedKap == min(MedKap) | 
           MedKap == max(MedKap))

# Fitted two regressions with a knot at N = 60
Kappa %>% 
  filter(SampleSize  < 60) %>%
  lm(kappa ~ SampleSize, data = .) %>%
  coef(.) %>%
  data.frame("values" = .) %>%
  mutate(valPer10 = values * 10)
Kappa %>% 
  filter(SampleSize  >= 60) %>%
  lm(kappa ~ SampleSize, data = .) %>%
  coef(.) %>%
  data.frame("values" = .) %>%
  mutate(valPer10 = values * 10)

# Let us use an iterative technique from the semgented package
fit <- lm(kappa ~ SampleSize, data = Kappa)
segm <- segmented(fit, seg.Z = ~ SampleSize)
slope(segm)
plot(segm)

ggplot(Kappa, aes(x = SampleSize, y = kappa)) + 
  geom_point()

# Plot results
ggplot(Kappa, aes(x = SampleSize, y = kappa)) +
geom_point() +
geom_smooth()

# Using boxplots
ggplot(Kappa, aes(x=factor(SampleSize), y = kappa)) + 
  geom_boxplot(outlier.size = .7, outlier.color = 'orange', size = 0.3) +
  scale_x_discrete(breaks = subjBreak, name = "Sample size") +
  scale_y_continuous(name='Kappa', breaks = seq(0,1,0.2)) +
  labs(caption = 'FDR = 0.05') +
  theme_bw()

# Set window 
quartz.options(width=18,height=12)

# Version OHBM 2018
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
KappaPlot <- ggplot(Kappa, aes(x=factor(SampleSize), y = kappa)) + 
  geom_boxplot(outlier.size = .7, outlier.color = 'orange', size = 0.3) +
  scale_x_discrete(breaks = subjBreak, name = "Sample size") +
  scale_y_continuous(name=expression(Coherence~~(kappa)), breaks = seq(0,1,0.2)) +
  labs(title = 'Concordance/coherence',
    subtitle = 'FDR = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
KappaPlot

# When is median K >= 0.80?
Kappa %>%
  group_by(SampleSize) %>%
  summarise(MedKap = median(kappa)) %>%
  filter(MedKap >= 0.80)

# Median K when N = 30?
Kappa %>%
  group_by(SampleSize) %>%
  summarise(MedKap = median(kappa)) %>%
  filter(SampleSize == 30)


##
###############
### Calculate Cohen's Kappa: FDR + all contrasts
###############
##

# Filter the final estimates
KappaC <- EMParam %>% filter(final == TRUE) %>%
  filter(threshold == 'fdr') %>%
  # Remove columns that don't provide information
  select(-num.iter, -final, -sequence) %>%
  # Transform step to sample size
  mutate(SampleSize = step * 10) %>%
  # Calculate Kappa
  mutate(kappa = NeuRRoStat::CohenKappa(lambda = lambda, piA1 = PI1, piI1 = PI2)) %>%
  # Put NaN to 0
  mutate(kappa = ifelse(is.nan(kappa), 0, kappa)) %>%
  # Remove columns
  select(run, SampleSize, kappa, contrast)

# Labels for contrasts
KappaC$contrastL <- factor(KappaC$contrast, levels = contrast,
                          labels = c('cogntive', 
                                     'faces',
                                     'MID',
                                     'Stop and Signal'))

# Plot with boxplots for all contrasts
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
KappaPlotC <- ggplot(KappaC, aes(x=factor(SampleSize), y = kappa)) + 
  geom_boxplot(outlier.size = .2, outlier.color = 'orange', size = 0.25,
               aes(fill = contrastL), width = 1) +
  scale_x_discrete(breaks = subjBreak, name = "Sample size") +
  scale_y_continuous(name=expression(Coherence~~(kappa)), breaks = seq(0,1,0.2)) +
  scale_fill_brewer('task ', type = 'qual', palette = 2) +
  labs(title = 'Concordance/coherence',
       subtitle = 'FDR = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        legend.text = element_text(size = 9, face = 'plain'),
        legend.title = element_text(size = 10, face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
KappaPlotC

# Plot with median +/- 1e and 3e quantile
cohQPlot4P <- 
  KappaC %>%
  group_by(SampleSize, contrastL) %>%
  summarise(avgCoh = mean(kappa, na.rm = TRUE),
            medCoh = median(kappa, na.rm = TRUE),
            Q1 = quantile(kappa, probs = 0.25, na.rm = TRUE),
            Q3 = quantile(kappa, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>% 
  ggplot(., aes(x = SampleSize, y = medCoh)) +
  geom_line(aes(colour = contrastL), size = 0.9) +
  geom_crossbar(aes(x = SampleSize, ymin = Q1, ymax = Q3,
                    fill = contrastL),
                colour = 'black', 
                size = 0.2,
                alpha = 0.5) +
  scale_x_continuous(name="Sample size") +
  scale_y_continuous(name=expression(Coherence~~(kappa))) +
  scale_fill_brewer('contrast ', type = 'qual', palette = 6) +
  scale_colour_brewer('contrast ', type = 'qual', palette = 6) +
  labs(title = 'Concordance/coherence',
       subtitle = 'FDR = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 13, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 11, face = 'bold'),
        legend.position = 'bottom')
cohQPlot4P

##### Some numbers #####
# When is median kappa >= 80?
KappaC %>%
  group_by(SampleSize, contrast) %>%
  summarise(medCoh = median(kappa, na.rm = TRUE)) %>% 
  filter(medCoh >= 0.80) 

# Maximum median kappa for MID and StopGo?
KappaC %>%
  group_by(SampleSize, contrast) %>%
  summarise(medCoh = median(kappa, na.rm = TRUE)) %>% 
  ungroup() %>% group_by(contrast) %>%
  filter(medCoh == max(medCoh))

# Values at N = 30
KappaC %>%
  group_by(SampleSize, contrast) %>%
  summarise(medCoh = median(kappa, na.rm = TRUE)) %>% 
  ungroup() %>% group_by(contrast) %>%
  filter(SampleSize == 30)


##
###############
### Save plots
###############
##

# Only the cognitive, FDR plot
ggsave(filename = paste0(getwd(), '/1_cognitive/Kappa_FDR_05.png'),
       plot = KappaPlot,
       width = 20, height = 14, units = 'cm', scale = 0.9)

# All contrasts in boxplots
ggsave(filename = paste0(getwd(), '/Kappa_FDR_05_contrasts.png'),
       plot = KappaPlotC,
       width = 20, height = 14, units = 'cm', scale = 0.9)


# Using the quantiles
ggsave(filename = paste0(getwd(), '/Kappa_FDR_05_quantiles.png'),
       plot = cohQPlot4P,
       width = 24, height = 22, units = 'cm', scale = 1)




##
###############
### Calculate Cohen's Kappa: comparing uncorrected, FDR and the one from Thirion
###############
##

# First calculate Kappa as written in text of Thirion et al. (2007)
# Filter the final estimates
Kappa_both <- EMParam %>% filter(final == TRUE & contrast == 'ML') %>% 
  # Remove columns that don't provide information
  select(-num.iter, -final, -sequence) %>%
  # Transform step to sample size
  mutate(SampleSize = step * 10) %>%
  # Calculate Kappa
  mutate(kappa = NeuRRoStat::CohenKappa(lambda = lambda, piA1 = PI1, piI1 = PI2)) %>%
  # Calculate Kappa according to Thirion
  mutate(kappa_Thirion = ThirionKappa(lambda = lambda, piA1 = PI1, piI1 = PI2)) %>%
  # Put NaN to 0
  mutate(kappa = ifelse(is.nan(kappa), 0, kappa)) %>%
  # Remove columns
  select(run, SampleSize, kappa, kappa_Thirion, threshold) %>%
  # gather results
  gather(key = 'typeOfkappa', value = 'value', 3:4)


# Plot the results
ggplot(Kappa_both, aes(x=factor(SampleSize), y = value, fill = threshold)) + 
  geom_boxplot(aes(x = factor(SampleSize), fill = threshold), 
               outlier.size = .4, outlier.color = 'orange', size = 0.3) +
  scale_x_discrete(breaks = subjBreak, name = "Sample size") +
  scale_y_continuous(name='Kappa', breaks = seq(0,1,0.2)) +
  scale_fill_brewer('Threshold', labels = c('FDR 0.05', 'Uncorrected 0.001'), 
                    type = 'qual') +
  facet_grid(typeOfkappa ~ .) +
  theme_bw() +
  theme(legend.position = 'bottom')


# Figure seems not good. Probably the text contained an error, but not the caluclation.

##
###############
### Calculate Cohen's Kappa: comparing uncorrected with FDR
###############
##

Kappa_both %>% 
  filter(typeOfkappa == 'kappa' & contrast == 'ML') %>%
  ggplot(., aes(x=factor(SampleSize), y = value, fill = threshold)) + 
    geom_boxplot(aes(x = factor(SampleSize), fill = threshold), 
                 outlier.size = .4, outlier.color = 'orange', size = 0.3) +
    scale_x_discrete(breaks = subjBreak, name = "Sample size") +
    scale_y_continuous(name='Kappa', breaks = seq(0,1,0.2)) +
    scale_fill_brewer('Threshold', labels = c('FDR 0.05', 'Uncorrected 0.001'), 
                      type = 'qual', palette = 2) +
    theme_bw() +
    theme(legend.position = 'bottom')


# Median kappa >= 0.80 in both thresholding levels?
Kappa_both %>% 
  filter(typeOfkappa == 'kappa' & contrast == 'ML') %>%
  group_by(SampleSize, threshold) %>%
  summarise(MedKap = median(value)) %>%
  filter(MedKap >= 0.80)



##
###############
### Plot lamda: estimated proportion of truly active voxels
###############
##

subjBreak <- c(10, seq(50, NTOT, by = 50))
EMParam %>% filter(final == TRUE & contrast == 'ML') %>% 
  filter(threshold == 'fdr') %>%
  # Remove columns that don't provide information
  select(-num.iter, -final, -sequence) %>%
  # Transform step to sample size
  mutate(SampleSize = step * 10) %>%
  # now plot
  ggplot(., aes(x = factor(SampleSize), y = lambda)) + 
  geom_boxplot(outlier.size = .7, outlier.color = 'orange', size = 0.3) +
  scale_x_discrete(name = "Sample size", breaks = subjBreak) +
  scale_y_continuous(name=
              expression(Proportion~of~truly~active~voxels~~(lambda))) +
  labs(title = 'Proportion of activated voxels',
       subtitle = 'FDR = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')


  