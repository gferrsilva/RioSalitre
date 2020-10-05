#####
# Code of Quality Control processing of portable X-Ray Fluorescence analysis
# Script follows
# 
# version: 1.2 (2020/10/05)
#
# Last modifications: Code cleaning, import('ggplot2')
#
# -----
# Process flow: 
#               
#               
# -----
# Guilherme Ferreira, (guilherme.ferreira@cprm.gov.br)
# June, 2020
#####

#####
# Setting up the enviroment
#####

# Working Directory
setwd("~/GitHub/RioSalitre")


# Objects
recr <- data.frame(xmin = 22.5, xmax = 36.5, ymin = -Inf, ymax = Inf)

#####
# Import Packages
#####

library(tidyverse)
library(Cairo)
library(janitor)
library(purrr)
library(ggplot2)

#####
# Built-in Functions
#####

qc_trans <- function(x, breaks, keep = NULL){
  # Function to corret translational bias
  # 3 arguments, x: Dataframe, matrix or similar
  # breaks: integer, index imediatly before the translational bias
  # keep: string vector with one or more numeric variable to keep out from analysis
  require(tidyverse)
  corfactor <- NULL
  mean1 <- NULL
  mean2 <- NULL
  df <- select(x, where(is.numeric)) %>%
    select(-all_of(keep))
  nome <- names(df)
    for (i in seq_along(df)) {
      a <- sum(df[1:breaks,i],na.rm = T)/length(which(!is.na(df[1:breaks,i])))
      b <- sum(df[((breaks+1):nrow(df)),i], na.rm = T)/length(which(!is.na(df[((breaks+1):nrow(df)),i])))
      corfactor[i] <- a/b
      mean1[i] <- a
      mean2[i] <- b
    }
  return(tibble(Element = nome,
                Factor = round(corfactor,digits = 2),
                Mean_Before_Break = round(mean1,digits = 2),
                Mean_After_Break = round(mean2,digits = 2)))
}
#####
# Import and Process Data
#####

# File with standard analysis
CRM <- read_tsv("Data/CRM_standard.txt") %>%
  mutate_at(18:length(.), as.numeric) %>%
  select(-c(7:17)) %>%
  mutate(Time = excel_numeric_to_date(Time,include_time = T),
         Index = seq.int(nrow(.)),
         `Reading_N` = as.factor(`Reading No`), `Reading No` = NULL) %>%
  select(-ends_with(c('O','O2','O3','O5'),ignore.case = F)) %>%
  select(1,53,2:52)

# File
df <- read_tsv("Data/QAQC.txt") %>%
  mutate_at(c(20,22,44,46,50,58), as.numeric)

# File with measurements
raw <- read_tsv("Data/Salitre_raw.txt")

#####
# Correcting the bias on CRM analysis
#####

# Creating the table of correction factors
cortable <- qc_trans(CRM %>%
                       select(-ends_with('rror')), breaks = 22, keep = 'Duration')

# Sr
CRM$Sr[23:36] <- (CRM$Sr[23:36])* cortable[[4,2]]

Sr_upper <- CRM$Sr + CRM$Sr_Error
Sr_lower <- CRM$Sr - CRM$Sr_Error

# Fe
CRM$Fe[23:36] <- (CRM$Fe[23:36])*cortable[[18,2]]

Fe_upper <- CRM$Fe + CRM$Fe_Error
Fe_lower <- CRM$Fe - CRM$Fe_Error

# Zr
CRM$Zr[23:36] <- (CRM$Zr[23:36])*cortable[[3,2]]

Zr_upper <- CRM$Zr + CRM$Zr_Error
Zr_lower <- CRM$Zr - CRM$Zr_Error

# Cu
CRM$Cu[23:36] <- (CRM$Cu[23:36])*cortable[[15,2]]

Cu_upper <- CRM$Cu + CRM$Cu_Error
Cu_lower <- CRM$Cu - CRM$Cu_Error

# Mo
CRM$Mo[23:36] <- (CRM$Mo[23:36])*cortable[[2,2]]

Mo_upper <- CRM$Mo + CRM$Mo_Error
Mo_lower <- CRM$Mo - CRM$Mo_Error

#####
# Data Vis
#####
CairoPDF(file = 'Figures/QAQC.PDF', width = 6, height = 4)
(Sr <- ggplot(CRM, aes(x = Index, y = Sr)) +
   geom_rect(inherit.aes = F, data = recr,
             mapping = aes(xmin = xmin,
                           xmax = xmax,
                           ymin = ymin,
                           ymax = ymax),
             fill = 'gray',
             alpha = .7) +
   geom_hline(yintercept = 109, size = .5, col = 'red') +
   geom_hline(yintercept = 87.20, linetype = 'dashed', col = 'red') +
   geom_hline(yintercept = 130.8, linetype = 'dashed', col = 'red') +
   geom_errorbar(data = CRM, aes(x = Index, ymin = Sr_lower, ymax = Sr_upper)) +
   geom_point(col = 'purple') +
   ylim(c(50,150)) +
   ylab('Sr (ppm)') +
   xlab('Standard Analytical Position') +
   labs(subtitle = 'RM 180-646 (Till-4)') +
   annotate('text', x = 29.5, y = 75,
            label = 'Corr. x1.021',
            col = 'black',
            size = 3)
)

(Fe <- ggplot(CRM, aes(x = Index, y = Fe)) +
    geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'grey',
              alpha = .7) +
    geom_hline(yintercept = 39700, size = .5, col = 'red') +
    geom_hline(yintercept = 35730, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 43670, linetype = 'dashed', col = 'red') +
    geom_errorbar(data = CRM, aes(x = Index, ymin = Fe_lower, ymax = Fe_upper)) +
    geom_point(col = 'purple') + 
    ylim(c(30000,45000)) +
    ylab('Fe (ppm)') +
    xlab('Standard Analytical Position') +
    labs(subtitle = 'RM 180-646 (Till-4)') +
    annotate('text', x = 29.5, y = 34000,
             label = 'Corr. x1.025',
             col = 'black',
             size = 3)
)

(Zr <- ggplot(CRM, aes(x = Index, y = Zr)) +
    geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'gray',
              alpha = .7) +
    geom_hline(yintercept = 385, size = .5, col = 'red') +
    geom_hline(yintercept = 346.5, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 423.5, linetype = 'dashed', col = 'red') +
    geom_errorbar(data = CRM, aes(x = Index, ymin = Zr_lower, ymax = Zr_upper)) +
    geom_point(col = 'purple') + 
    ylim(c(300,450)) +
    ylab('Zr (ppm)') +
    xlab('Standard Analytical Position') +
    labs(subtitle = 'RM 180-646 (Till-4)') +
    annotate('text', x = 29.5, y = 320,
             label = 'Corr. x1.099',
             col = 'black',
             size = 3)
)

(Cu <- ggplot(CRM, aes(x = Index, y = Cu)) +
    geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'gray',
              alpha = .7) +
    geom_hline(yintercept = 237, size = .5, col = 'red') +
    geom_hline(yintercept = 189.6, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 284.4, linetype = 'dashed', col = 'red') +
    geom_errorbar(data = CRM, aes(x = Index, ymin = Cu_lower, ymax = Cu_upper)) +
    geom_point(col = 'purple') + 
    ylim(c(150,350)) +
    ylab('Cu (ppm)') +
    xlab('Standard Analytical Position') +
    labs(subtitle = 'RM 180-646 (Till-4)') +
    annotate('text', x = 29.5, y = 335,
             label = 'Corr. x1.099',
             col = 'black',
             size = 3)
)

(Mo <- ggplot(CRM, aes(x = Index, y = Mo)) +
    geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'gray',
              alpha = .7) +
    geom_hline(yintercept = 16, size = .5, col = 'red') +
    geom_hline(yintercept = 19.2, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 12.8, linetype = 'dashed', col = 'red') +
    geom_errorbar(data = CRM, aes(x = Index, ymin = Mo_lower, ymax = Mo_upper)) +
    geom_point(col = 'purple') + 
    ylim(c(0,25)) +
    ylab('Mo (ppm)') +
    xlab('Standard Analytical Position') +
    labs(subtitle = 'RM 180-646 (Till-4)') +
    annotate('text', x = 29.5, y = 5,
             label = 'Corr. x1.084',
             col = 'black',
             size = 3)
)
dev.off()
#####
# Applying the corrections on other analysis
#####

## Correcting Fe Values
df[725:1021,50] <-  as.double(df$Fe[725:1021])* cortable[[18,2]]

## Correcting Zr Values
df[725:1021,20] <-  as.double(df$Zr[725:1021]) * cortable[[3,2]]

## Correcting Sr Values
df[725:1021,22] <-  as.double(df$Sr[725:1021]) * cortable[[4,2]]

## Correcting Cu Values
df[725:1021,44] <-  as.double(df$Cu[725:1021]) * cortable[[15,2]]

## Correcting Ni Values
df[725:1021,46] <-  as.double(df$Ni[725:1021]) * cortable[[16,2]]

## Correcting Ti Values
df[725:1021,58] <-  as.double(df$Ti[725:1021]) * cortable[[22,2]]


# Preparing data to export ----
df <- df %>%
  filter(Index %in% raw$Index)

df <- bind_cols(raw[1:7], df[18:133]) %>%
  select(-ends_with('Error')) %>%
  select(-c(55:63))

# Writing file
write_delim(df,'Data/salitre_corr.txt',delim = '/t')
