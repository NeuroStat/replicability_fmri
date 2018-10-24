####################
#### TITLE:     Plot figure with perfect versus non-perfect coherence
#### Contents: 	
#### 
#### Source Files: NIMPP_Freddie
#### First Modified: 13/09/2018
#### Notes: 
#################

##
###############
### Notes
###############
##

# Create observations for two plots, one with perfect and other with
# mixed coherence.

##
###############
### Preparation
###############
##

# Seed
set.seed(1112)

# Load in libraries
library(tidyverse)
library(RColorBrewer)
library(NeuRRoStat)

# Number of observations
N <- 10000

# Number of classes
Nc <- 10

##
###############
### Create data
###############
##

# First create perfect coherence
PCoh <- sample(x = c(0,10), size = N, replace = TRUE, prob = c(.60,.40))
perfCoh <- data.frame(value = PCoh) %>% as_tibble() %>%
  mutate(class = ifelse(value == 0, 'True Inactive', 'True Active'))
perfCoh$class <- as.factor(perfCoh$class)
# True kappa with these values
TrueK_PC <- NeuRRoStat::CohenKappa(lambda = .60, piA1 = 1, piI1 = 0)

# Now create less perfect coherence
#   to achieve, we need a mixture distribution
# Therefore, we first create the label, then each replication with a probability
# of success, then sum all replications
# lambda equals prob of belonging to null
lambda <- 0.75
# Prob of success in null
prob_null <- 0.8
# prob of success in alternative
prob_act <- 0.5
# True kappa with these values
TrueK_MC <- NeuRRoStat::CohenKappa(lambda = lambda, piA1 = prob_act, piI1 = 1 - prob_null)



# First create vector with labels
labels <- data.frame(class = c(rep('True Inactive', length = lambda * N), 
                     rep('True Active', length = (1 - lambda) * N))) %>% as_tibble()

# Empty data frame
mixCoh <- labels %>% mutate(observation = 0)
# For loop to create the data
for(i in 1:Nc){
  # Create data for this replication
  repl <- labels %>% 
    mutate(realization = c(rbinom(n = lambda * N, size = 1, 
                                  prob = (1 - prob_null)),
                           rbinom(n = (1- lambda) * N, size = 1, 
                                  prob = prob_act)))
  
  # Add to empty data frame
  mixCoh$observation <- mixCoh$observation + repl$realization
}


##
###############
### Plot results
###############
##

# Save with dimensions: 570 x 550

# Plot perfect coherence
ggplot(data = perfCoh, aes(x = value)) +
  geom_histogram(aes(fill = class), binwidth = 0.5) + 
  scale_x_continuous('Sum of binary - replicated images') +
  scale_y_continuous('Number of voxels') +
  labs(title = 'Perfect coherence',
       subtitle = 'Number of independent replications = 10') +
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

# Plot mixed coherence
ggplot(data = mixCoh, aes(x = observation)) +
  geom_histogram(aes(fill = class), position = 'dodge2', binwidth = .5) + 
  scale_x_continuous('Sum of binary - replicated images') +
  scale_y_continuous('Number of voxels') +
  labs(title = 'Without perfect coherence',
       subtitle = 'Number of independent replications = 10') +
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

# In one panel
# First create labels of the expressions
labels_text <- data.frame(
  label = c('lambda == 0.60', 'lambda == 0.75',
            'pi[A]^1 == 1', 'pi[A]^1 == 0.5',
            'pi[I]^1 == 0', 'pi[I]^1 == 0.2'),
  sourceF = c('kappa == 1', 'kappa == 0.21',
              'kappa == 1', 'kappa == 0.21',
              'kappa == 1', 'kappa == 0.21'),
  x = c(9,9,9,9,9,9),
  y = c(5900,5900,5500,5500,5100,5100)
)
# Now plot
mixCoh %>% 
	rename(value = observation) %>%
	mutate(source = 'kappa == 0.21') %>%
	bind_rows(cbind(perfCoh, 'source' = 'kappa == 1')) %>%
	# Add order variable
	mutate(order = ifelse(source == 'kappa == 1', 1, 2)) %>%
	# now add order to source
	mutate(sourceF = reorder(source, order)) %>%
  ggplot(data = ., aes(x = value)) +
  geom_histogram(aes(fill = class), position = 'dodge2', binwidth = 1) + 
  scale_x_continuous('Sum of binarized - replicated images',
  breaks = c(0,2,4,6,8,10)) +
  scale_y_continuous('Number of voxels') +
  scale_fill_manual('class', labels = c('True Active', 'True Inactive'),
                    values = c('#00BFC4', '#F8766D')) +
  facet_wrap(~ sourceF, labeller = label_parsed) + 
  labs(title = 'Illustration of coherence') +
  # Add parameters to the panels
  geom_text(
    data = labels_text,
    mapping = aes(x = x, y = y, label = label),
    size = 3,
    fontface = 'bold',
    hjust = 0,
    parse = TRUE) +
  # Add box
  annotate("rect", xmin = 8.5, xmax = 10.5, 
           ymin = 4900, ymax = 6100,
           alpha = .2) + 
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text.y = element_text(size = 0),
        axis.text = element_text(size = 13, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 13),
        plot.subtitle = element_text(hjust = 0, vjust = -2),
        legend.key.size = unit(.80,'cm'),
        legend.text = element_text(size = 11),
        legend.position = 'bottom',
        legend.title = element_text(size = 11))  

# Save the figure
ggsave(filename = paste0(getwd(),'/IllustrationCoh.png'), plot = last_plot(),
scale = 1, width = 20, height = 14, units = 'cm', dpi = 600)





































