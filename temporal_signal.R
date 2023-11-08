rm(list = ls())

setwd('/Users/ltj8/Documents/burk_projects/ms_burk/better_subpanel/')

library(reshape2)
library(scales)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(stringr)
library(plyr)
library(stats)

temp <- read.table(file = "ms_burk_subpanel_tempest_data.txt", 
                   sep = '\t', header = TRUE)
head(temp)
attach(temp)


slope <- lm(formula = distance ~ date, data = temp)$coefficients[2]

intercept <- lm(formula = distance ~ date, data = temp)$coefficients[1]

ggplot(temp, aes(x=date, y=distance)) + theme_light() + 
  geom_point() + theme_classic() +
  geom_abline(aes(slope = slope, intercept=intercept, colour = "red")) + 
  xlab("Year") + ylab("Root-to-tip Divergence") + 
  # scale_x_continuous(breaks = c(1990,1995,2000,2005,2010,2015,2020)) +
  scale_y_continuous(position = "left") + 
  theme(axis.title.y = element_text(size=25)) +
  theme(axis.title.x = element_text(size=25)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.text.x = element_text(size=10)) 

dev.print(pdf, 'ms_burk_pruned_snps_temporal_signal.pdf')

# png(file="2023_oct_burk_subpanel_temporal_signal.png",
#     width=600, height=350)

