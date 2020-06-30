#####
# Code of data processing of the article: Comparison of K-Means and Model-Based
# Clustering methods for drill core pseudo-log generation based on X-Ray Fluorescence Data
# 
# version: 2.0 (2020/06/23)
#
# Last modifications: Code cleaning
#
# -----
# Process flow: Import pXRF data (after QAQC), Data-cleaning, Missing value imputation,
#               Moving Average Filtering, PCA, Elbow Analysis, Silhouette,
#               K-Means, Model-Based-Cluster
# -----
# Guilherme Ferreira, (guilherme.ferreira@cprm.gov.br)
# June, 2020
#####

#####
# Setting up the enviroment
#####

setwd("~/Salitre")
set.seed(123)

#####
# Import Packages
#####

library(tidyverse) # ggplot, dplyr, readr, tibble, readr
library(hrbrthemes) # Color palletes and themes
library(Cairo) # Export Plots
library(TTR) # Exponential Movel Average
library(corrplot) # Correlation analysis
library(factoextra) # K-Means Clustering and PCA
library(mclust) # Model-Based Clustering

#####
# Built-in Functions
#####

minmax.norm <- function(df, Keep = NULL) {
  # Função para normalizar dados pelo range de mínimo e máximo
  # Argumento df = dataframe com dados a serem normalizados, por coluna
  # Argumento Keep = vetor de strings ou números com nome ou posição de colunas que não serão normalizadas
  require(tidyverse)
  normalize <- function(x) {
    return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = T) - min(x, na.rm = T)))
  }
  if(is.null(Keep)){
    index <- names(df)
    dflabs <- df %>% discard(is.numeric)
    dfnum <- df %>% keep(is.numeric)
    dfnum <- as_tibble(sapply(dfnum, normalize))
    df <- as_tibble(cbind(dflabs, dfnum))
    df %>%
      select(all_of(index))
  } else {
    index <- names(df)
    k <- df %>%
      select(all_of(Keep))
    df <- df %>%
      select(-all_of(Keep))
    dflabs <- df %>% discard(is.numeric)
    dfnum <- df %>% keep(is.numeric)
    dfnum <- as_tibble(sapply(dfnum, normalize))
    df <- as_tibble(cbind(dflabs, k, dfnum))
    df %>%
      select(all_of(index))
  }
}

filter.var <- function(df, cut = .95) { 
  # Filter the columns by the number of non-NA rows.
  # 2 args: df and cut (default .95)
    require(tidyverse)
  df %>%
    select_if(~sum(!is.na(.x)) >= (cut * nrow(df)))  
}

imput.var <- function(df, coef = 0.5) {
  # Imputation in missing values
  require(tidyverse)
  index <- names(df)
  dflabs <- df %>% discard(is.numeric)
  dfnum <- df %>% keep(is.numeric)
  
  min.var <- {}
  
  for (i in 1:length(dfnum)) {
    min.var[i] <- min(dfnum[i], na.rm = TRUE)  
  }
  
  for (i in 1:length(dfnum)) {
    dfnum[[i]] <- replace_na(dfnum[[i]], coef*min.var[[i]])
  }
  
  df <- as_tibble(cbind(dflabs, dfnum))
  df %>%
    select(all_of(index)) %>%
    mutate_if(is.numeric, function(x) ifelse(is.infinite(x), NA, x))
}

# Extending ggplot functions. Source [1]:
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     required_aes = c("x", "y"))

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#####
# Import and Process Data
#####

# Importing ----

# Raw data from pXRF txt file
raw <- read_tsv("Salitre_full_corr.txt") %>% # Raw data from txt file
  separate(col = "SAMPLE", into = c("CORE", "FROM", "ID"), sep = "-") %>%
  select(-c("Date", "Duration", "Time", "ID", "Bal", "K")) %>%
  mutate_at(c(3,5:51), as.numeric, length = 2) %>% # Converting columns to numeric
  mutate_at(c(1:2,4), as.factor)
# Key values of lithotypes
litho <- read_tsv("litho_apa3001.txt",col_types = list('f','d','f'))

# Summary ----

# Statistical Summary
df_sum <- as.data.frame(sapply(sapply(raw, summary)[5:51], cbind))
row.names(df_sum) <- c('Min', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.', 'NA')
write.csv(df_sum, "statisticalsummary.csv")

# Data wrangling ----

## df_long: untransformed dataset long format  
df_long <-  raw %>%
  filter(rowMeans(is.na(.)) < .7) %>%
  filter.var() %>% # See line 75
  imput.var() %>% # See line 83
  mutate(FROM = as.numeric(FROM, length = 2)) %>%
  left_join(litho, by = c('CORE','FROM')) %>%
  gather(key = "element", value = "measure", 5:24, na.rm = F)
## df_norm: transformed and normalized dataset
df_norm <-  raw %>%
  filter(rowMeans(is.na(.)) < .7) %>%
  filter.var() %>% # See line 75
  mutate(FROM = as.numeric(FROM, length = 2)) %>%
  left_join(litho, by = c('CORE','FROM')) %>%
  imput.var() %>% # See line 83
  minmax.norm(Keep = c('FROM','Index')) %>% # See line 44
  select(1:4, 25, 5:24)
## long: transformed dataset
long <- df_norm %>%
  gather(key = "element", value = "measure", 6:25, na.rm = F)

# Filtering Data ----

df_fil <- df_norm %>%
  filter(TYPE == 'S') %>%
  select(6:25) %>%
  lapply(EMA, n = 5) %>%
  bind_cols(df_norm[1:5] %>%
              filter(TYPE == 'S')) %>%
  filter(!is.na(Fe)) %>%
  select(21:25, 1:24)


#####
# Data Analysis
#####

# Hypothesis Test ----

## Normality Test, Shapiro-Wilk (Shapiro & Wilk, 1965; Razali & Wah, 2011)
t(do.call(rbind,
          lapply(X = df_norm %>% select(-CORE, -FROM, -TYPE, -Index, -LITHO),
                 function(x) shapiro.test(x)[c("statistic", "p.value")])))

## Variance comparison of groups, Kruskal-Wallis (Kruskal & Wallis, 1952)
t(do.call(rbind,
          lapply(X = df_norm %>% select(-CORE, -FROM, -Index, -LITHO, - TYPE),
                 g = factor(df_norm$TYPE),
                 function (x, g) kruskal.test(x, g)[c("statistic",
                                                      "p.value")])))

# Correlation Matrix ----

## Raw data
M_r <- cor(df_norm %>% 
             filter(TYPE == 'S') %>%
             select(6:25, -Si, -P), method = "spearman")
corrplot(M_r, method="color",# col=col(200),  
         type="lower", order= 'original',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=FALSE, tl.pos = 'n',cl.pos = 'n',
         tl.cex = 1, # tamanho da fonte dos labels
         tl.offset = 1.1, # offset do nome das colunas em relação a matriz
         number.font = 1.5,
         number.digits = 1)

## Filtered data
M_f <- cor(df_fil %>% 
             filter(TYPE == 'S') %>%
             select(6:25, -Si, -P),
           method = "spearman")
corrplot(M_f, method="color",# col=col(200),  
         type="upper", order= 'original', 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=FALSE, 
         tl.cex = 1, # tamanho da fonte dos labels
         tl.offset = 1, # offset do nome das colunas em relação a matriz
         number.font = 1.5,
         number.digits = 1)

# Principal Component Analysis ----

## Raw data
pca_raw <- prcomp(df_norm %>% 
                    filter(TYPE == 'S') %>%
                    select(6:25, -Si, -P), center = TRUE, scale. = F)
summary(pca_raw)
### Plot
fviz_pca_ind(pca_raw,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F,      # Avoid text overlapping
             axes = c(1,2), #controla quais eixos devem ser mostrados na figura
             geom = c("point")
)

fviz_pca_var(pca_raw,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

## Filtered data
pca_fil <- prcomp(df_fil[c(6:20,23:25)], center = TRUE, scale. = F)
summary(pca_fil)
### Plot
fviz_pca_ind(pca_fil,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F,      # Avoid text overlapping
             axes = c(1,2), #controla quais eixos devem ser mostrados na figura
             geom = c("point")
)
fviz_pca_var(pca_fil,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)



#####
# Unsupervised Learning
#####

# Hyperparameter tunning ----

## Elbow method
### Raw data
fviz_nbclust(df_norm %>% 
               filter(TYPE == 'S') %>%
               select(6:25, -Si, -P), kmeans, method = "wss", k.max = 25)
### Filtered data
fviz_nbclust(df_fil %>%
               select(6:25, -Si, -P), kmeans, method = "wss", k.max = 25)

## Silhoutte method
### Raw data
fviz_nbclust(df_norm %>% 
               filter(TYPE == 'S') %>%
               select(6:25, -Si, -P), kmeans, method = "silhouette", k.max = 25)
### Filtered data
fviz_nbclust(df_fil %>%
               select(6:25, -Si, -P), kmeans, method = "silhouette", k.max = 25)

# K-Means ----

k2_raw <- kmeans(df_norm %>% 
                   filter(TYPE == 'S') %>%
                   select(6:25, -Si, -P), 
             centers = 2, 
             nstart = 50,iter.max = 100)
k2_fil <- kmeans(df_fil %>%
                   select(6:25, -Si, -P), 
                 centers = 2, 
                 nstart = 50,iter.max = 100)
## K-Means plot for raw data
fviz_cluster(k2_raw,
             data = pca_raw$x[,1:2],
             axes = c(1,2), stand = F, ellipse = T,
             geom = c('point'), main = 'K-Means Cluster (Raw Data)')
## K-Means plot for filtered data
fviz_cluster(k2_fil,
             data = pca_fil$x[,1:2],
             axes = c(1,2), stand = F, ellipse = T,
             geom = c('point'), main = 'K-Means Cluster (Filtered Data)')

# Model-Based Cluster ----

mclus_raw <- Mclust(pca_raw$x[,1:2])
summary(mclus_raw)
mclus_fil <- Mclust(pca_fil$x[,1:2])
summary(mclus_fil)
## MClus plot for raw data
fviz_cluster(mclus_raw,
             data = pca_raw$x,
             axes = c(1,2), stand = F, ellipse = T, 
             ellipse.type = 'convex', 
             choose.vars = c(1,2), geom = c('point'),
             main = 'MB Cluster (Raw Data)')
## MClus plot for filterd data
fviz_cluster(mclus_fil,
             data = pca_fil$x,
             axes = c(1,2), stand = F, ellipse = T, 
             ellipse.type = 'convex', 
             choose.vars = c(1,2), geom = c('point'),
             main = 'MB Cluster (Raw Data)')

#####
# Data Visualization
#####


# ---- Exploratory graphs ----

## QQplot of 1st and 2nd Spot
p1 <- ggplot(long, aes(sample = measure, col = TYPE)) + 
  geom_qq() + geom_qq_line() + theme_grey() + theme(
    legend.position= 'bottom',
    panel.spacing = unit(.1, "lines"),
    strip.text.x = element_text(size = 12)
  ) + xlab("") + ylab("") +
  labs(title = 'Q-Q plot by element') +
  scale_fill_continuous(name = "Group",
                        labels = c("Duplicate", "Measure")) +
  facet_wrap(~element, scale = "free") +
  scale_color_discrete(name = "Group",
                       labels = c("2nd Spot", "Measure"))

## Density plot of 1st and 2nd Spot
p2 <- ggplot(long, aes(x = measure, y = ..density.., fill = TYPE)) +
  geom_density(alpha = .4) +
  labs(title = 'Probability density by element') + theme_grey() +  theme(
    legend.position= 'bottom',
    panel.spacing = unit(.1, "lines"),
    strip.text.x = element_text(size = 10)
  ) + facet_wrap(~element, scale = "free") + xlab("") + ylab("") + 
  scale_fill_discrete(name = 'Group', labels = c('2nd Spot', '1st Spot')
  )
## Boxplot for selected elements classified by 1st and 2nd Spot

p3 <- ggplot(long, aes(y = measure, x = element, col = TYPE)) + 
  geom_boxplot() + theme_grey() + theme(
    legend.position= 'top',
    panel.spacing = unit(.1, "lines"),
    strip.text.x = element_text(size = 12)
  ) + #facet_wrap(. ~ element, scale = 'free') + 
  xlab("") + ylab("") +
  labs(title = 'Boxplot by element') +
  scale_color_discrete(name = '', labels = c("2nd Spot", "1st Spot"))

# Plot of K-Means solution for raw data

df_fil %>%
  filter(TYPE == 'S') %>%
  bind_cols(mclus_fil['classification'], pca_fil$x) %>%
  ggplot(aes(x = PC1,
             y = PC2)) +
  stat_chull(aes(fill = factor(classification)),
             alpha = .4,
             col = 'white') +
  geom_point(size = 3.55,
             col = 'pink',
             alpha = .5) +
  geom_point(aes(shape = LITHO,
                 col = LITHO), 
             size = 3.5) +
  scale_shape_manual(values = c(20,20,18,18)) +
  scale_fill_manual(values = c('#453781FF',
                               '#A65C85FF',
                               '#EB8055FF',
                               '#B8DE29FF')) +
  scale_color_manual(values = c('#0C2A50FF',
                                '#6B4596FF',
                                '#CC6A70FF',
                                '#F7CB44FF')) +
  theme_bw()

# Plot of MBClus solution for filtered data

df_norm %>%
  filter(TYPE == 'S') %>%
  bind_cols(mclus_raw['classification'], pca_raw$x) %>%
  ggplot(aes(x = PC1,
             y = PC2)) +
  stat_chull(aes(fill = factor(classification)),
             alpha = .3, col = 'black') +
  geom_point(aes(fill = LITHO), 
             size = 5,
             shape = 21,
             col = 'black',
             alpha = .6)


CairoPDF(file = 'Figures/Vetores/exploratory_graphs.PDF', width = 8, height = 6)
print(p1)
print(p2)
print(p3)
dev.off()


#####
# REFERENCES
#####

# Shapiro, S. S.; Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). Biometrika. 52 (3-4): 591-611. doi:10.1093/biomet/52.3-4.591. JSTOR 2333709. MR 0205384. p. 593
# Kruskal, William H.; Wallis, W. Allen (1 de dezembro de 1952). Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association. 47 (260): 583-621. ISSN 0162-1459. doi:10.1080/01621459.1952.10483441
# Razali, Nornadiah; Wah, Yap Bee (2011). "Power comparisons of Shapiro-Wilk, Kolmogorov-Smirnov, Lilliefors and Anderson-Darling tests". Journal of Statistical Modeling and Analytics. 2 (1): 21-33
# [1] https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
# [2] https://medium.com/apprentice-journal/pca-application-in-machine-learning-4827c07a61db
# [3] https://uc-r.github.io/kmeans_clustering
# [4] http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
