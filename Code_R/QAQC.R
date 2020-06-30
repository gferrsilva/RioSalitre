setwd("~/Salitre")
#####################################
# Carregando e preparando dados
#####################################

library(readr)
library(tidyverse)
library(Cairo)


CRM <- read_tsv("CRM_salitre.txt")

#Corrigindo o tipo de variável no banco de dados

for (i in 20:length(CRM)) {
  CRM[[i]] <- as.numeric(CRM[[i]])
}

CRM$Time <- as.POSIXct.Date(CRM$Time, TZ = 'GMT -3')

CRM$Index <- seq.int(nrow(CRM))
## Getting the values for the correction

TILL <- CRM

TILL_label <- TILL[1:5]
TILL_num <- TILL[18:133] %>% select(-ends_with('Error'))

TILL <- bind_cols(TILL_label, TILL_num[-c(48:56)])

corfactor <- NULL
for (i in 6:(length(TILL))) {
  a <- sum(TILL[1:22,i],na.rm = T)/length(which(!is.na(TILL[1:22,i])))
  b <- sum(TILL[23:36,i], na.rm = T)/length(which(!is.na(TILL[23:36,i])))
  corfactor[i-5] <- a/b
}

mean1 <- NULL
for (i in 6:(length(TILL))) {
  a <- sum(TILL[1:22,i],na.rm = T)/length(which(!is.na(TILL[1:22,i])))
  mean1[i-5] <- a
}

mean2 <- NULL
for (i in 6:(length(TILL))) {
  b <- sum(TILL[23:36,i], na.rm = T)/length(which(!is.na(TILL[23:36,i])))
  mean2[i-5] <- b
}
elem_list <- names(TILL_num[-c(48:56)])

cortable <- data.frame(Element = elem_list, Correction = corfactor, Mean1 = mean1, Mean2 = mean2)
elems <- c("Zr","Sr","Cu","Ni","Fe","Cr","V","Ti","S","Nd","Pr","Ba","Sb","Sn","Cd","Cl","Cs","Te")

table <- subset(cortable, Element %in% elems)
write.csv(table, 'tables/qaqc.csv')

library(ggplot2)

recr <- data.frame(xmin = 22.5, xmax = 36.5, ymin = -Inf, ymax = Inf)
### QAQC para Sr
Sr_corfac <- mean(CRM$Sr[1:22])/mean(CRM$Sr[23:36])
CRM$Sr[23:36] <- Sr_corfac*(CRM$Sr[23:36])

Sr_upper <- CRM$Sr + CRM$Sr_Error
Sr_lower <- CRM$Sr - CRM$Sr_Error


(Sr <- ggplot(CRM, aes(x = Index, y = Sr)) +
  geom_rect(inherit.aes = F, data = recr,
            mapping = aes(xmin = xmin,
                          xmax = xmax,
                          ymin = ymin,
                          ymax = ymax),
            fill = 'white', col = 'gray',
            alpha = .8) +
  geom_hline(yintercept = 109, size = .5, col = 'red') +
  geom_hline(yintercept = 87.20, linetype = 'dashed', col = 'red') +
  geom_hline(yintercept = 130.8, linetype = 'dashed', col = 'red') +
  # geom_hline(yintercept = 103.55, linetype = 'dashed', col = 'red') +
  # geom_hline(yintercept = 114.45, linetype = 'dashed', col = 'red') +
  geom_errorbar(data = CRM, aes(x = Index, ymin = Sr_lower, ymax = Sr_upper)) +
  geom_point(col = 'purple') +
  ylim(c(50,150)) +
  # geom_hline(yintercept = mean(CRM$Sr), size = 1, col = 'purple')
  ylab('Sr (ppm)') +
  xlab('Standard Analytical Position') +
  labs(subtitle = 'RM 180-646 (Till-4)') +
  # annotate('text', x = 18,
  #          y = c(112,123,96),
  #          label = c('CRM Value', '+10%','-10%'),
  #          col = 'red', size = 3) +
  annotate('text', x = 29.5, y = 75,
           label = 'Corr. x1.021',
           col = 'black',
           size = 3)
)

CairoPDF(file = 'Figures/Vetores/QAQC/Sr.PDF', width = 4, height = 3)
print(Sr)
dev.off()

print(mean(CRM$Sr[1:22])/mean(CRM$Sr[23:36]))

### QAQC para Fe

Fe_corfac <- mean(CRM$Fe[1:22])/mean(CRM$Fe[23:36])
CRM$Fe[23:36] <- Fe_corfac*(CRM$Fe[23:36])
 
Fe_upper <- CRM$Fe + CRM$Fe_Error
Fe_lower <- CRM$Fe - CRM$Fe_Error

(Fe <- ggplot(CRM, aes(x = Index, y = Fe)) +
    geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'grey',
              alpha = .5) +
    geom_hline(yintercept = 39700, size = .5, col = 'red') +
    geom_hline(yintercept = 35730, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 43670, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = 37715, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = 41685, linetype = 'dashed', col = 'red') +
    geom_errorbar(data = CRM, aes(x = Index, ymin = Fe_lower, ymax = Fe_upper)) +
  geom_point(col = 'purple') + 
  ylim(c(30000,45000)) +
  ylab('Fe (ppm)') +
  xlab('Standard Analytical Position') +
  labs(subtitle = 'RM 180-646 (Till-4)') +
  # annotate('text', x = 18,
  #          y = c(38800,44300,35000),
  #          label = c('CRM Value',
  #                    '+10%','-10%'),
  #          col = 'red', size = 3) +
  annotate('text', x = 29.5, y = 34000,
           label = 'Corr. x1.025',
           col = 'black',
           size = 3)
)

CairoPDF(file = 'Figures/Vetores/QAQC/Fe.PDF', width = 4, height = 3)
print(Fe)
dev.off()

print(mean(CRM$Fe[1:22])/mean(CRM$Fe[23:36]))

### QAQC para Zr
Zr_corfac <- mean(CRM$Zr[1:22])/mean(CRM$Zr[23:36])
CRM$Zr[23:36] <- Zr_corfac*(CRM$Zr[23:36])


Zr_upper <- CRM$Zr + CRM$Zr_Error
Zr_lower <- CRM$Zr - CRM$Zr_Error

(Zr <- ggplot(CRM, aes(x = Index, y = Zr)) +
    geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'gray',
              alpha = .5) +
    geom_hline(yintercept = 385, size = .5, col = 'red') +
    geom_hline(yintercept = 346.5, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 423.5, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = 365.75, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = 404.25, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = mean(CRM$Fe), size = 1, col = 'purple')
    geom_errorbar(data = CRM, aes(x = Index, ymin = Zr_lower, ymax = Zr_upper)) +
  geom_point(col = 'purple') + 
  ylim(c(300,450)) +
  ylab('Zr (ppm)') +
  xlab('Standard Analytical Position') +
  labs(subtitle = 'RM 180-646 (Till-4)') +
  # annotate('text', x = c(25,18,18),
  #          y = c(391,430,340),
  #          label = c('CRM Value',
  #                    '+10%','-10%'),
  #          col = 'red', size = 3) +
  annotate('text', x = 29.5, y = 320,
           label = 'Corr. x1.099',
           col = 'black',
           size = 3)
)

CairoPDF(file = 'Figures/Vetores/QAQC/Zr.PDF', width = 4, height = 3)
print(Zr)
dev.off()

print(mean(CRM$Zr[1:22])/mean(CRM$Zr[23:36]))

### QAQC para Cu
Cu_corfac <- mean(CRM$Cu[1:22])/mean(CRM$Cu[23:36])
CRM$Cu[23:36] <- Cu_corfac*(CRM$Cu[23:36])

Cu_upper <- CRM$Cu + CRM$Cu_Error
Cu_lower <- CRM$Cu - CRM$Cu_Error

(Cu <- ggplot(CRM, aes(x = Index, y = Cu)) +
  geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'white', col = 'gray',
              alpha = .8) +
  geom_hline(yintercept = 237, size = .5, col = 'red') +
  geom_hline(yintercept = 189.6, linetype = 'dashed', col = 'red') +
  geom_hline(yintercept = 284.4, linetype = 'dashed', col = 'red') +
  geom_errorbar(data = CRM, aes(x = Index, ymin = Cu_lower, ymax = Cu_upper)) +
  geom_point(col = 'purple') + 
  ylim(c(200,350)) +
  ylab('Cu (ppm)') +
  xlab('Standard Analytical Position') +
  labs(subtitle = 'RM 180-646 (Till-4)') +
  # annotate('text', x = c(30,18,18),
  #          y = c(232,205,255),
  #          label = c('CRM Value',
  #                    '+10%','-10%'),
  #          col = 'red', size = 3) +
  annotate('text', x = 29.5, y = 335,
             label = 'Corr. x1.099',
             col = 'black',
             size = 3)
)

CairoPDF(file = 'Figures/Vetores/QAQC/Cu.PDF', width = 4, height = 3)
print(Cu)
dev.off()


### QAQC para Mo
Mo_corfac <- mean(CRM$Mo[1:22])/mean(CRM$Mo[23:36])
CRM$Mo[23:36] <- Mo_corfac*(CRM$Mo[23:36])

Mo_upper <- CRM$Mo + CRM$Mo_Error
Mo_lower <- CRM$Mo - CRM$Mo_Error

(Mo <- ggplot(CRM, aes(x = Index, y = Mo)) +
  geom_rect(inherit.aes = F, data = recr,
              mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              fill = 'white', col = 'gray',
              alpha = .8) +
    geom_hline(yintercept = 16, size = .5, col = 'red') +
    geom_hline(yintercept = 19.2, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 12.8, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = 22.4, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = 9.6, linetype = 'dashed', col = 'red') +
    # geom_hline(yintercept = mean(CRM$Fe), size = 1, col = 'purple')
    geom_errorbar(data = CRM, aes(x = Index, ymin = Mo_lower, ymax = Mo_upper)) +
    geom_point(col = 'purple') + 
    ylim(c(0,25)) +
    ylab('Mo (ppm)') +
    xlab('Standard Analytical Position') +
    labs(subtitle = 'RM 180-646 (Till-4)') +
    # annotate('text', x = c(25,18,18),
    #          y = c(17.2,8.5,24),
    #          label = c('CRM Value',
    #                    '-2sigma','+2sigma'),
    #          col = 'red', size = 3)+
    annotate('text', x = 29.5, y = 5,
             label = 'Corr. x1.084',
             col = 'black',
             size = 3)
)

CairoPDF(file = 'Figures/Vetores/QAQC/Cu.PDF', width = 4, height = 3)
print(Mo)
dev.off()

df <- read_tsv("QAQC.txt")

## Correcting Fe Values
df[725:1021,50] <-  as.double(df$Fe[725:1021]) * 1.0259685

## Correcting Zr Values
df[725:1021,20] <-  as.double(df$Zr[725:1021]) * 1.0993535

## Correcting Sr Values
df[725:1021,22] <-  as.double(df$Sr[725:1021]) * 1.0217436

## Correcting Cu Values
df[725:1021,44] <-  as.double(df$Cu[725:1021]) * 1.0242120

## Correcting Ni Values
df[725:1021,46] <-  as.double(df$Ni[725:1021]) * 1.1082672

## Correcting Ti Values
df[725:1021,58] <-  as.double(df$Ti[725:1021]) * 1.0668877

raw <- read_tsv("Salitre_full.txt")

ind <- raw$Index

df1 <- df %>% filter(Index %in% ind)

labs <- raw[1:7]

df <- bind_cols(labs, df1[18:133])

df <- df %>% 
  select(-ends_with('Error'))

df <- df[-c(55:63)]

write_delim(df,'salitre_full_corr.txt')
