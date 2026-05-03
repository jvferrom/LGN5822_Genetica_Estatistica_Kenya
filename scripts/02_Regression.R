# ==============================================================================
# ATIVIDADE EM GRUPO 02 - Regressões (Kenya)
# Dataset: Pan-African Soybean Variety Trials (PAT)
# Referência: Araújo et al. (2025), Scientific Data
# ==============================================================================

# PACOTES ----------------------------------------------------------------------
library(tidyverse)   
library(ggpubr)      
library(corrplot)    
library(Hmisc)       
library(knitr)      
library(kableExtra)  
library(broom)       
library(car)         
library(lmtest)      
library(patchwork)  
library(scales)

# Criação de pastas
dir_figuras <- "../results/figures/02_regression/"
dir_tabelas <- "../results/tables/02_regression/"

# Criar as pastas se elas não existirem (usando os nomes corretos das variáveis)
if (!dir.exists(dir_figuras)) dir.create(dir_figuras, recursive = TRUE)
if (!dir.exists(dir_tabelas)) dir.create(dir_tabelas, recursive = TRUE)


# ==============================================================================
# 1. IMPORTAÇÃO E PRÉ-PROCESSAMENTO
# ==============================================================================

dat <- read.csv(
  "https://raw.githubusercontent.com/mauricioaraujj/Pan_African_Trials_Network/refs/heads/main/data/data.csv",
  sep = ";"
)
# Dimensões
cat("Linhas:", nrow(dat), "| Colunas:", ncol(dat), "\n")

# Selecionar e limpar variáveis de interesse
vars_interesse <- c("GY", "PH_R8", "W100G", "NDM", "FLW_DAYS",
                    "PROT", "OIL", "ELEV", "LAT", "LON", "COUNTRY")

dat_clean <- dat %>%
  select(all_of(vars_interesse)) %>%
  # Filtrar apenas observações do Quênia
  filter(COUNTRY == "Kenya") %>%
  # Substituir zeros em traits que não podem ser zero por NA
  # (PROT e OIL são NA nas linhas sem análise NIR; mantidos para uso pairwise)
  mutate(
    GY       = ifelse(GY       == 0, NA, GY),
    PH_R8    = ifelse(PH_R8    == 0, NA, PH_R8),
    W100G    = ifelse(W100G    == 0, NA, W100G),
    NDM      = ifelse(NDM      == 0, NA, NDM),
    FLW_DAYS = ifelse(FLW_DAYS == 0, NA, FLW_DAYS)
  )
cat("País: Kenya\n")
cat("Observações após filtragem e limpeza:", nrow(dat_clean), "\n")
cat("Locais únicos:", n_distinct(interaction(dat_clean$LAT, dat_clean$LON)), "\n")

# ==============================================================================
# 2. Exploração dos Dados
# ==============================================================================
# Estatísticas descritivas
tabela_descritiva <- dat_clean %>%
  select(-COUNTRY) %>%
  pivot_longer(everything(), names_to = "variavel", values_to = "valor") %>%
  group_by(variavel) %>%
  summarise(
    n       = sum(!is.na(valor)),
    media   = mean(valor, na.rm = TRUE),
    dp      = sd(valor, na.rm = TRUE),
    min     = min(valor, na.rm = TRUE),
    q25     = quantile(valor, 0.25, na.rm = TRUE),
    mediana = median(valor, na.rm = TRUE),
    q75     = quantile(valor, 0.75, na.rm = TRUE),
    max     = max(valor, na.rm = TRUE),
    cv_pct  = (sd(valor, na.rm = TRUE) / mean(valor, na.rm = TRUE)) * 100
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

# Salvando a tabela em CSV
write.csv(tabela_descritiva, file.path(dir_tabelas, "estatisticas_descritivas.csv"), row.names = FALSE)

# Imprimindo a tabela no console
tabela_descritiva %>%
  kable(
    caption   = "Estatísticas descritivas das variáveis de interesse",
    col.names = c("Variável", "n", "Média", "DP", "Mín", "Q1",
                  "Mediana", "Q3", "Máx", "CV (%)"),
    align     = "lrrrrrrrrr"
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                latex_options     = c("striped", "hold_position"),
                full_width = FALSE)

# Distribuições das variáveis
p_hist <- dat_clean %>%
  select(-COUNTRY) %>%
  pivot_longer(everything(), names_to = "variavel", values_to = "valor") %>%
  ggplot(aes(x = valor, fill = variavel)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.85) +
  facet_wrap(~variavel, scales = "free", ncol = 4) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#2C3E50"),
        strip.text = element_text(color = "white", face = "bold")) +
  labs(title = "Distribuição das variáveis agronômicas e ambientais",
       x = "Valor", y = "Frequência")

print(p_hist)

# Salvando o gráfico de distribuição
ggsave(file.path(dir_figuras, "distribuicao_variaveis.png"), plot = p_hist, width = 12, height = 8, dpi = 300)

# Matriz de correlação
# Variáveis numéricas com pelo menos 30 observações válidas
traits_numericos <- names(dat_clean)[sapply(dat_clean, is.numeric)]
disponibilidade  <- sapply(traits_numericos, function(t)
  sum(!is.na(dat_clean[[t]])))
traits_ok <- names(disponibilidade[disponibilidade >= 30])
label_uso <- traits_ok

# Correlações de Pearson com p-valores via rcorr (Hmisc)
cor_data  <- dat_clean %>% select(all_of(traits_ok)) %>% as.matrix()
res_rcorr <- Hmisc::rcorr(cor_data, type = "pearson")
cor_matrix <- res_rcorr$r
p_matrix   <- res_rcorr$P

# Montar data.frame longo para ggplot2
cor_long <- as.data.frame(cor_matrix) %>%
  mutate(Var1 = rownames(cor_matrix)) %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "r")

p_long <- as.data.frame(p_matrix) %>%
  mutate(Var1 = rownames(p_matrix)) %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "p_val")

cor_final <- cor_long %>%
  left_join(p_long, by = c("Var1", "Var2")) %>%
  mutate(
    i = match(Var1, label_uso),
    j = match(Var2, label_uso)
  ) %>%
  filter(i <= j) %>%
  select(-i, -j) %>%
  mutate(
    p_val      = ifelse(is.na(p_val), 1, p_val),
    stars      = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    r_label    = ifelse(Var1 == Var2, "1.00",
                        ifelse(is.na(r), "NA",
                               sprintf("%.2f\n%s", r, stars))),
    text_color = ifelse(abs(r) > 0.6 & !is.na(r), "white", "black"),
    Var1 = factor(Var1, levels = label_uso),
    Var2 = factor(Var2, levels = rev(label_uso))
  )

p_cor <- ggplot(cor_final, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = r_label, color = text_color),
            size = 3.2, fontface = "bold", lineheight = 0.8) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#c0392b", mid = "white", high = "#2471a3",
    midpoint = 0, limits = c(-1, 1),
    name = "Correlação\n(r)", na.value = "grey85"
  ) +
  scale_x_discrete(position = "bottom") +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5),
    axis.title    = element_blank(),
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Correlação de Pearson entre Traits (Kenya)",
    subtitle = "Significância: *** p<0.001  ** p<0.01  * p<0.05"
  )

print(p_cor)

# Salvando a matriz de correlação
ggsave(file.path(dir_figuras, "correlacao_pearson.png"), plot = p_cor, width = 10, height = 8, dpi = 300)


# ==============================================================================
# 3. Regressão Linear Simples
# ==============================================================================

