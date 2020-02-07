# K-Means applied to fine-grained rocks clustering

This codes is part of the "Principal Component Analysis and K-Means clustering based on X-Ray Fluorescence Data in fine-grained meta volcano-sedimentary rocks from the Rio Salitre Greenstone Belt, Brazil" written by the team* of the Geological Survey of Brazil.

*Guilherme Ferreira da Silva
correspondent author

Directory of Geology and Mineral Resources, Geological Survey of Brazil – CPRM, 
ORCiD https://orcid.org/0000-0002-3675-7289
guilherme.ferreira@cprm.gov.br

  # 1. SETTINGS
  ## 1.1 Reading Data
``` R
library(readr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(corrplot)
library(factoextra)
set.seed(123)
```
  ## 1.2 Data Preparation
Discarding irrelevant Factor Variables

``` R
df_raw <- read_tsv("Salitre.txt")

trash <- c("Index", "Date", "Duration", "CORE", "Time")

df_raw <- df_raw %>%
  select(-trash)

# Defining variables as numbers
for (i in 3:length(df_raw)) {
  df_raw[[i]] <- as.numeric(df_raw[[i]])
 }
```

Statistical summary

``` R
df_summary <- t(do.call(rbind,lapply(df_raw, summary)))
write.csv(df_summary, "df_summary.csv")
```

Filtering for variables that have at least 75% os values higher than lower detection limit
``` R
cut <- .75
df <- df_raw %>%
  select_if(~sum(!is.na(.x)) >= (cut * nrow(df_raw)))
```
Replacing <LDL values to 1/2 the minimun of each variable
``` R
minimo <- {}
for (i in 3:length(df)) {
  minimo[i] <- min(df[[i]], na.rm = TRUE)  
}

for (i in 3:length(df)) {
  df[[i]] <- replace_na(df[[i]], .5*minimo[[i]])
}
```

# 2. DATA SCALING 

## 2.1 Amplitude normalization
``` R
normalize <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = T) - min(x, na.rm = T)))
}
labels <- df %>% select(SAMPLE, TYPE)

df_norm <- df %>%
  select(-SAMPLE) %>%
  select(-TYPE) %>%
  sapply(FUN = normalize)

df_norm <- bind_cols(labels, as_data_frame(df_norm))

head(df_norm)
```
Output:
```
# A tibble: 6 x 27
  SAMPLE TYPE     Zr    Sr     W    Cu    Ni     Fe     Cr      V     Ti     Sc      K     S    Nd
  <chr>  <chr> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
1 apa20~ S     0.205 0.266 0.599 0.451 0.549 0.0793 0.0843 0.0862 0.0618 0.0943 0.0201 0.159 0.356
2 apa20~ D     0.198 0.255 0.590 0.452 0.628 0.0828 0.0857 0.114  0.0524 0.0650 0.0182 0.135 0.427
3 apa20~ S     0.137 0.173 0.529 0.368 0.475 0.129  0.0705 0.152  0.0726 0.0631 0.0210 0.140 0.524
4 apa20~ D     0.152 0.180 0.531 0.387 0.519 0.113  0.0752 0.144  0.0763 0.0854 0.0217 0.136 0.444
5 apa20~ S     0.203 0.300 0.538 0.416 0.473 0.118  0.0811 0.202  0.0990 0.136  0.0219 0.127 0.585
6 apa20~ D     0.161 0.446 0.503 0.372 0.483 0.0839 0.0719 0.0982 0.0639 0      0.0197 0.129 0.386
# ... with 12 more variables: Pr <dbl>, La <dbl>, Ba <dbl>, Sb <dbl>, Sn <dbl>, Cd <dbl>, Bal <dbl>,
#   P <dbl>, Si <dbl>, Cl <dbl>, Cs <dbl>, Te <dbl>
```
## 2.2 Log Transformation
``` R
labels <- df %>% select(SAMPLE, TYPE)

df_log <- df %>%
  select(-SAMPLE) %>%
  select(-TYPE) %>%
  sapply(FUN = log10)

df_log <- bind_cols(labels, as_data_frame(df_log))
```
# 3. DATA MANIPULATION
## 3.1 Splitting samples and duplicates
``` R
sam <- df_norm %>%
  filter(TYPE == "S")

dup <- df_norm %>%
  filter(TYPE == "D")
```

## 3.2 Converting in Long Data Frame
``` R
descarte <- c("Cl", "P", "Bal", "Sc", "Sn", "Cd", "Sb", "Te", "Ba", "Cs") # ver código abaixo

sam <- sam %>% select(-descarte)

dup <- dup %>% select(-descarte)

s_long <- sam %>%
  gather(key = "element", value = "measure", 3:length(sam), na.rm = TRUE)
d_long <- dup %>%
  gather(key = "element", value = "measure", 3:length(sam), na.rm = TRUE)

long <- full_join(s_long, d_long) %>%
  rename(Group = TYPE)
```

## 3.3 Selecting only the numerical variables

``` R
df_num <- select(sam,-SAMPLE) %>%
  select(-TYPE)

write.csv(df_num, "salitre_norm.csv")
```
# 4. STATISTICAL SIGNIFICANCE VERIFICATION

## 4.1 Check of distribution
Shapiro-Wilk test (Shapiro & Wilk, 1965; Razali & Wah, 2011)

For alpha defined as 5%, this test checks if the data sample is normally distributed
Null hypothesis (H0): The data sample is normally distributed
Alternative hypothesis (H1): The data sample has another distribution

``` R
t(do.call(rbind, 
          lapply(X = df %>% select(-SAMPLE, -TYPE), 
                 function(x) shapiro.test(x)[c("statistic", "p.value")]
          )))
```
Output:
```
          Zr           Sr           W            Cu           Ni           Fe           Cr          
statistic 0.8922746    0.6953068    0.8343817    0.8757685    0.9050847    0.8931482    0.6540958   
p.value   8.357584e-16 8.777475e-26 1.293297e-19 5.112535e-17 9.089255e-15 9.770249e-16 3.519157e-27
          V            Ti           Sc           K            S           Nd           Pr          
statistic 0.8077025    0.768408     0.6508929    0.7614976    0.6405334   0.9818664    0.9724514   
p.value   4.828141e-21 6.770853e-23 2.776015e-27 3.387864e-23 1.30379e-27 9.504289e-05 1.150519e-06
          La           Ba           Sb           Sn           Cd           Bal         P           
statistic 0.9364157    0.9769057    0.9512631    0.9617951    0.9707346    0.3565725   0.931491    
p.value   9.635029e-12 8.306087e-06 6.013767e-10 1.881325e-08 5.620156e-07 1.40641e-34 2.816869e-12
          Si           Cl           Cs           Te        
statistic 0.4730813    0.7995515    0.9845759    0.9905358 
p.value   4.281468e-32 1.893626e-21 0.0004077478 0.01435178

```

## 4.2 Equivalence of data and duplicates
Kruskal-Wallis (Kruskal & Wallis, 1952)

For alhpa defined as 5%, this test for non-parametric data checks if sample and duplicates are originated from the same distribution. 
Null hypothesis (H0): Sample and duplicate have the same dristibution
Alternative hypothesis (H1): Sample and duplicate do not come from the same dristibution
``` R
t(do.call(rbind,
          lapply(X = df %>% select(-SAMPLE),
                 g = factor(df$TYPE),
                 function (x, g) kruskal.test(x, g)[c("statistic",
                                                      "p.value")])))
```
Output:

```
          TYPE       Zr         Sr         W         Cu        Ni        Fe        Cr        
statistic 383        0.06917833 0.02232867 0.7173547 0.1936311 0.2933421 1.135893  0.07686608
p.value   2.7658e-85 0.7925369  0.8812161  0.397013  0.6599115 0.5880867 0.2865218 0.7815902 
          V         Ti        Sc        K         S         Nd        Pr         La         Ba       
statistic 0.1944266 0.2265028 2.487362  0.284446  2.002724  0.8635429 0.08772136 0.03776917 0.1948168
p.value   0.6592577 0.6341291 0.1147639 0.5938019 0.1570168 0.3527491 0.7670942  0.8459074  0.6589375
          Sb        Sn        Cd         Bal        P          Si        Cl        Cs        Te      
statistic 0.1516571 0.560968  0.07235625 3.591698   0.06819939 1.719424  1.859628  1.505758  1.815892
p.value   0.6969568 0.4538705 0.7879365  0.05806888 0.7939767  0.1897673 0.1726677 0.2197875 0.177803
```

# 5. DATA VISUALIZATION

# 5.1 QQ-plot for selected elements

``` R
ggplot(long, aes(sample = measure, col = Group)) + 
  geom_qq(alpha = .6) + geom_qq_line() + theme(
    legend.position= c(.85,.1),
    panel.spacing = unit(.1, "lines"),
    strip.text.x = element_text(size = 12)
  ) + xlab("") + ylab("") +
  labs(title = 'Q-Q plot by element') +
  scale_fill_continuous(name = "Group", labels = c("Duplicate", "Measure")) +
  facet_wrap(~element, scale = "free") + scale_color_discrete(name = "Group", labels = c("Duplicate", "Measure"))
```
[!png] 
# 5.2 Density plot of Sample and Duplicate

``` R
ggplot(long, aes(x = measure, y = ..density.., fill = Group)) +
  geom_density(alpha = .6) +
  labs(title = 'Probability density by element') +  theme(
    legend.position= c(.85,.05),
    panel.spacing = unit(.1, "lines"),
    strip.text.x = element_text(size = 10)
  ) + facet_wrap(~element, scale = "free") + xlab("") + ylab("") + 
  scale_fill_discrete(
    name = "Group",
    labels = c("Duplicate", "Measure")
  )
```




#####################################
# Correlograma
#####################################

# Rodando o correlograma: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
  ## Após uma primeira análise no correlograma, reparei que os elementos abaixo
  ## tinham uma alta correlação entre si, mas não tinha significado aparente
descarte <- c("Cl", "P", "Bal", "Sc", "Sn", "Cd", "Sb", "Te", "Ba", "Cs")

M <- cor(df_num, method = "spearman")

library(corrplot)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color",# col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         tl.cex = 1, # tamanho da fonte dos labels
         tl.offset = 1.1, # offset do nome das colunas em relação a matriz
         number.font = 2,
         number.digits = 2
)

#####################################
# Análise de Componentes Principais
#####################################

## Script gerado com auxílio do tutorial disponível em: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

df_pca <- prcomp(df_num, center = TRUE, scale. = TRUE)

summary(df_pca)

library(factoextra)

# Visualizando os Autovalores
fviz_eig(df_pca)

# Gráfico de indivíduos: aqueles de perfil similar serão agrupados conjuntamente

fviz_pca_ind(df_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F,      # Avoid text overlapping
             axes = c(1,2), #controla quais eixos devem ser mostrados na figura
             geom = c("point")
             )

# Gráfico de variáveis: variáveis correlacionáveis apontam para a mesma direção

fviz_pca_var(df_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Gráfico de indivíduos e variáveis

fviz_pca_biplot(df_pca, 
                #palette = "jco", 
                addEllipses = TRUE, label = "var")
               # col.var = "contrib", repel = TRUE,
                #legend.title = "Species") 

#####################################
# Hyperparameter tunning
#####################################

## Elbow method

set.seed(123)

fviz_nbclust(df_num, kmeans, method = "wss", k.max = 25)
## Silhoutte method

fviz_nbclust(df_num, kmeans, method = "silhouette", k.max = 25)

#####################################
# K-Means Clustering
#####################################

## Distance matrix

distance <- get_dist(df_num, method = "manhattan")


fviz_dist(distance, gradient = list(
  low = "#00AFBB", mid = "white", high = "#FC4E07"),
  show_labels = T, lab_size = 5, order = TRUE)

## Clustering

k3 <- kmeans(df_num, 
             centers = 3, 
             nstart = 50)

df_clus <- as.data.frame(k3["cluster"])

write.csv(k3["cluster"], "cluster.csv")

fviz_cluster(k3,
             data = df_num,
             axes = c(1,2))
# REFERENCES
## Papers
Kruskal, William H.; Wallis, W. Allen (1 de dezembro de 1952). Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association. 47 (260): 583–621. ISSN 0162-1459. doi:10.1080/01621459.1952.10483441

Razali, Nornadiah; Wah, Yap Bee (2011). "Power comparisons of Shapiro–Wilk, Kolmogorov–Smirnov, Lilliefors and Anderson–Darling tests". Journal of Statistical Modeling and Analytics. 2 (1): 21–33

Shapiro, S. S.; Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). Biometrika. 52 (3–4): 591–611. doi:10.1093/biomet/52.3-4.591. JSTOR 2333709. MR 0205384. p. 593

## Websites
Correlation Matrix
http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram

PCA
https://towardsdatascience.com/dimensionality-reduction-does-pca-really-improve-classification-outcome-6e9ba21f0a32

PCA for Machine Learning
https://medium.com/apprentice-journal/pca-application-in-machine-learning-4827c07a61db

K-Means Clustering
https://uc-r.github.io/kmeans_clustering
