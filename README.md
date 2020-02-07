# K-Means-clustering-based-in-pXRF-data-in-fine-grained-rocks

setwd("C:/Users/guilherme.ferreira/Desktop/Salitre")

## Carregando e preparando dados

library(readr)
library(tidyverse)

df_raw <- read_tsv("Salitre.txt")

#Discartando colunas que não interessam

trash <- c("Index", "Date", "Duration", "CORE", "Time")

df_raw <- df_raw %>%
  select(-trash)

# Gerando o Resumo estatístico

df_summary <- t(do.call(rbind,lapply(df_raw, summary)))
write.csv(df_summary, "df_summary.csv")

# Corrigindo o tipo de variável no banco de dados

for (i in 3:length(df_raw)) {
  df_raw[[i]] <- as.numeric(df_raw[[i]])
}

# Selecionando as variáveis que possuem mais de cut de valores acima do limite de detecção

cut <- .75
df <- df_raw %>%
  select_if(~sum(!is.na(.x)) >= (cut * nrow(df_raw)))

# Substituindo os valores ausentes por metade do mínimo de cada coluna

minimo <- {}
for (i in 3:length(df)) {
  minimo[i] <- min(df[[i]], na.rm = TRUE)  
}

for (i in 3:length(df)) {
  df[[i]] <- replace_na(df[[i]], .5*minimo[[i]])
}

# Normalizando os dados 

### Normalizando pela amplitude
normalize <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = T) - min(x, na.rm = T)))
}
labels <- df %>% select(SAMPLE, TYPE)

df_norm <- df %>%
  select(-SAMPLE) %>%
  select(-TYPE) %>%
  sapply(FUN = normalize)

df_norm <- bind_cols(labels, as_data_frame(df_norm))
# Normalizando os dados pela pela Transformação Log10
labels <- df %>% select(SAMPLE, TYPE)

df_log <- df %>%
  select(-SAMPLE) %>%
  select(-TYPE) %>%
  sapply(FUN = log10)

df_log <- bind_cols(labels, as_data_frame(df_log))

#####################################
# Separando banco de dados
#####################################

dup <- df_norm %>%
  filter(TYPE == "D")

sam <- df_norm %>%
  filter(TYPE == "S")

# Tranformando em Long Data Frame
descarte <- c("Cl", "P", "Bal", "Sc", "Sn", "Cd", "Sb", "Te", "Ba", "Cs") # ver código abaixo

sam <- sam %>% select(-descarte)

dup <- dup %>% select(-descarte)

s_long <- sam %>%
  gather(key = "element", value = "measure", 3:length(sam), na.rm = TRUE)
d_long <- dup %>%
  gather(key = "element", value = "measure", 3:length(sam), na.rm = TRUE)

long <- full_join(s_long, d_long) %>%
  rename(Group = TYPE)

#####################################
# Testes de Significância estatística
#####################################

#Teste de Normalidade de Shapiro-Wilk test (Shapiro & Wilk, 1965; Razali & Wah, 2011)
t(do.call(rbind, 
          lapply(X = df %>% select(-SAMPLE, -TYPE), 
                 function(x) shapiro.test(x)[c("statistic", "p.value")]
          )))

#Teste de variância em postos para dados não paramétricos de Kruskal-Wallis (Kruskal & Wallis, 1952)
t(do.call(rbind,
          lapply(X = df %>% select(-SAMPLE),
                 g = factor(df$TYPE),
                 function (x, g) kruskal.test(x, g)[c("statistic",
                                                      "p.value")])))

#####################################
# Visualizando os dados
#####################################

library(ggplot2)
library(hrbrthemes)

## QQplot do conjunto de dados, comparados a uma normal
ggplot(long, aes(sample = measure, col = Group)) + 
  geom_qq(alpha = .6) + geom_qq_line() + theme(
    legend.position= c(.85,.1),
    panel.spacing = unit(.1, "lines"),
    strip.text.x = element_text(size = 12)
  ) + xlab("") + ylab("") +
  labs(title = 'Q-Q plot by element') +
  scale_fill_continuous(name = "Group", labels = c("Duplicate", "Measure")) +
  facet_wrap(~element, scale = "free") + scale_color_discrete(name = "Group", labels = c("Duplicate", "Measure"))

## Density plot comparando as análises e as duplicatas
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



# Trabalhando com as variáveis numéricas
df_num <- select(sam,-SAMPLE) %>%
  select(-TYPE)

write.csv(df_num, "salitre_norm.csv")

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
#####################################
##REFERENCES
#####################################

# Shapiro, S. S.; Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). Biometrika. 52 (3–4): 591–611. doi:10.1093/biomet/52.3-4.591. JSTOR 2333709. MR 0205384. p. 593
# Kruskal, William H.; Wallis, W. Allen (1 de dezembro de 1952). Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association. 47 (260): 583–621. ISSN 0162-1459. doi:10.1080/01621459.1952.10483441
# Razali, Nornadiah; Wah, Yap Bee (2011). "Power comparisons of Shapiro–Wilk, Kolmogorov–Smirnov, Lilliefors and Anderson–Darling tests". Journal of Statistical Modeling and Analytics. 2 (1): 21–33

# https://towardsdatascience.com/dimensionality-reduction-does-pca-really-improve-classification-outcome-6e9ba21f0a32
# https://medium.com/apprentice-journal/pca-application-in-machine-learning-4827c07a61db
# https://uc-r.github.io/kmeans_clustering
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
