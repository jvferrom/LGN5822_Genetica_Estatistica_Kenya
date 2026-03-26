# ==============================================================================
# ATIVIDADE EM GRUPO 01 - Análise Exploratória de Dados (Kenya)
# Dataset: Pan-African Soybean Variety Trials (PAT)
# Referência: Araújo et al. (2025), Scientific Data
# ==============================================================================

# PACOTES ----------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(Hmisc)

# ==============================================================================
# 1. IMPORTAÇÃO E INSPEÇÃO INICIAL DOS DADOS
# ==============================================================================
dat <- read.csv(
  "https://raw.githubusercontent.com/mauricioaraujj/Pan_African_Trials_Network/refs/heads/main/data/data.csv",
  sep = ";"
)

cat("Dimensões do dataset completo:", dim(dat), "\n")
names(dat) <- toupper(names(dat))

#==============================================================================
# 2. FILTRO PARA KENYA
# ==============================================================================
dat_kenya <- subset(dat, COUNTRY == "Kenya")
cat("\nDados do Kenya:", nrow(dat_kenya), "registros em",
    length(unique(dat_kenya$ENV)), "ambientes\n")

# ==============================================================================
# 3. AJUSTE DE DADOS
# ==============================================================================
# Substituição de 0 por NA
zero_is_na <- c("GY", "W100G", "FLW_DAYS", "PH_R8", "NDM", "PROT", "OIL",
                "LOD", "PL_EMERG_COUNT", "PL_EMERG_PCT", "POD_SHATTERING")
for (t in zero_is_na) {
  if (t %in% names(dat_kenya)) {
    n_zeros <- sum(dat_kenya[[t]] == 0, na.rm = TRUE)
    dat_kenya[[t]][dat_kenya[[t]] == 0] <- NA
    if (n_zeros > 0) cat("  [", t, "] →", n_zeros, "zeros convertidos para NA\n")
  }
}

# ==============================================================================
# 4. DICIONÁRIO DE VARIÁVEIS
# ==============================================================================
trait_dict <- c(
  "GY"             = "Produtividade de Grãos (kg/ha)",
  "W100G"          = "Peso de 100 Grãos (g)",
  "FLW_DAYS"       = "Dias até o Florescimento",
  "PH_R8"          = "Altura da Planta na Maturidade (cm)",
  "NDM"            = "Dias até a Maturidade",
  "PROT"           = "Teor de Proteína (%)",
  "OIL"            = "Teor de Óleo (%)",
  "LOD"            = "Acamamento (escala 1–5)",
  "PL_EMERG_COUNT" = "Contagem de Emergência",
  "PL_EMERG_PCT"   = "Porcentagem de Emergência (%)"
)

traits <- names(trait_dict)



# ==============================================================================
# 5. CONTROLE DE QUALIDADE E FILTRAGEM
# ==============================================================================
# Estatísticas de Qualidade Experimental (CV%) por Ambiente baseados em GY
qc_stats <- dat_kenya %>%
  group_by(ENV) %>%
  summarise(
    N_Parcelas  = sum(!is.na(GY)),
    Media_GY    = mean(GY, na.rm = TRUE),
    DP_GY       = sd(GY, na.rm = TRUE),
    CV_pct      = (DP_GY / Media_GY) * 100,
    .groups     = "drop"
  ) %>%
  mutate(
    CV_pct      = round(CV_pct, 1),
    Media_GY    = round(Media_GY, 1),
    Qualidade   = case_when(
      CV_pct <= 15              ~ "Excelente (≤15%)",
      CV_pct <= 30              ~ "Aceitável (15–30%)",
      TRUE                      ~ "Alto ruído (>30%)"
    )
  ) %>%
  arrange(desc(CV_pct))

cat("\n===== QUALIDADE EXPERIMENTAL POR AMBIENTE (GY) =====\n")
print(as.data.frame(qc_stats))

#Plot de CV%
p_cv <- ggplot(qc_stats, aes(x = reorder(ENV, CV_pct), y = CV_pct, fill = Qualidade)) +
  geom_col() +
  geom_hline(yintercept = 30, linetype = "dashed", color = "red",    linewidth = 0.8) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "orange", linewidth = 0.8) +
  annotate("text", x = 1, y = 32, label = "Limite 30%", hjust = 0, color = "red", size = 3) +
  annotate("text", x = 1, y = 17, label = "Limite 15%", hjust = 0, color = "orange", size = 3) +
  scale_fill_manual(values = c("Excelente (≤15%)" = "#2ecc71", "Aceitável (15–30%)" = "#f39c12", "Alto ruído (>30%)" = "#e74c3c")) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = "Coeficiente de Variação (CV%) por Ambiente – Kenya", subtitle = "Baseado em Produtividade de Grãos (GY)", x = "Ambiente", y = "CV (%)", fill = "Qualidade") +
  theme(plot.title = element_text(face = "bold"))
print(p_cv)

# Filtragem com base no CV% de GY
limite_cv <- 30
ambientes_aprovados <- qc_stats %>% filter(CV_pct <= limite_cv) %>% pull(ENV)
dat_kenya_filtrado <- dat_kenya %>% filter(ENV %in% ambientes_aprovados)

cat("\n===== RESUMO DA FILTRAGEM =====\n")
cat("Ambientes originais       :", length(unique(dat_kenya$ENV)), "\n")
cat("Ambientes mantidos (CV ≤", limite_cv, "%) :", length(ambientes_aprovados), "\n")


# ==============================================================================
# 6. EDA FOCADA em 4 Trials e 5 Traits)
# ==============================================================================
# Subset
trials_selecionados <- c("E011", "E0150", "E0194", "E001")
traits_selecionados <- c("GY", "W100G", "FLW_DAYS", "PH_R8", "NDM")
dat_foco_long <- dat_kenya_filtrado %>%
  filter(ENV %in% trials_selecionados) %>%
  select(ENV, all_of(traits_selecionados)) %>%
  pivot_longer(cols = all_of(traits_selecionados), names_to = "Trait", values_to = "Valor") %>%
  drop_na(Valor)

#Detecção de Outliers (Baseado em IQR)
outliers_detectados <- dat_foco_long %>%
  group_by(ENV, Trait) %>%
  mutate(
    Q1 = quantile(Valor, 0.25), Q3 = quantile(Valor, 0.75), IQR = Q3 - Q1,
    Limite_Inferior = Q1 - 1.5 * IQR, Limite_Superior = Q3 + 1.5 * IQR,
    Is_Outlier = Valor < Limite_Inferior | Valor > Limite_Superior 
  ) %>%
  filter(Is_Outlier == TRUE) %>%
  select(ENV, Trait, Valor) %>% arrange(Trait, ENV)

cat("\n===== RELATÓRIO DE DETECÇÃO DE OUTLIERS =====\n")
cat("Total de outliers encontrados nos 4 Trials selecionados:", nrow(outliers_detectados), "\n")
print(as.data.frame(outliers_detectados))

# Estatísticas de Tendência Central e Dispersão
stats_foco <- dat_foco_long %>%
  group_by(ENV, Trait) %>%
  summarise(Media = round(mean(Valor), 2), Mediana = round(median(Valor), 2), DP = round(sd(Valor), 2), .groups = "drop")
cat("\n===== ESTATÍSTICAS DE TENDÊNCIA CENTRAL E DISPERSÃO (4 Trials) =====\n")
print(as.data.frame(stats_foco))


# Histogramas e Gráficos de Densidade
p_hist_dens <- ggplot(dat_foco_long, aes(x = Valor, fill = ENV)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, alpha = 0.5, color = "black", position = "identity") +
  geom_density(alpha = 0.2, linewidth = 0.8, aes(color = ENV)) +
  facet_wrap(~ Trait, scales = "free", ncol = 3) +
  theme_minimal(base_size = 11) +
  labs(title = "Histogramas e Curvas de Densidade", subtitle = "Avaliação de 5 Traits em 4 Trials selecionados no Kenya", x = "Valor", y = "Densidade", fill = "Trial", color = "Trial") +
  theme(strip.text = element_text(face = "bold"))
print(p_hist_dens)


# ==============================================================================
# 7. RELAÇÕES E "ACROSS TRIALS" EDA (Usando os dados filtrados por CV%)
# ==============================================================================
# MATRIZ DE CORRELAÇÃO ENTRE TRAITS
traits_numericos <- intersect(traits, names(dat_kenya_filtrado))
disponibilidade <- sapply(traits_numericos, function(t) sum(!is.na(dat_kenya_filtrado[[t]])))
traits_ok <- names(disponibilidade[disponibilidade >= 30])

cor_data <- dat_kenya_filtrado %>% select(all_of(traits_ok)) %>% as.matrix()
res_rcorr <- rcorr(cor_data, type = "pearson")
cor_matrix <- res_rcorr$r
p_matrix <- res_rcorr$P
label_uso <- traits_ok

cor_long <- as.data.frame(cor_matrix) %>% mutate(Var1 = rownames(cor_matrix)) %>% pivot_longer(-Var1, names_to = "Var2", values_to = "r")
p_long <- as.data.frame(p_matrix) %>% mutate(Var1 = rownames(p_matrix)) %>% pivot_longer(-Var1, names_to = "Var2", values_to = "p_val")

cor_final <- cor_long %>%
  left_join(p_long, by = c("Var1", "Var2")) %>%
  mutate(i = match(Var1, label_uso), j = match(Var2, label_uso)) %>%
  filter(i <= j) %>% select(-i, -j) %>%
  mutate(
    p_val = ifelse(is.na(p_val), 1, p_val),
    stars = case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ ""),
    r_label = ifelse(Var1 == Var2, "1.00", ifelse(is.na(r), "NA", sprintf("%.2f\n%s", r, stars))),
    Var1 = factor(Var1, levels = label_uso), Var2 = factor(Var2, levels = rev(label_uso))
  )

p_cor <- ggplot(cor_final, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = r_label), size = 3.2, fontface = "bold", lineheight = 0.8, color = ifelse(abs(cor_final$r) > 0.6 & !is.na(cor_final$r), "white", "black")) +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2471a3", midpoint = 0, limits = c(-1, 1), name = "Correlação\n(r)", na.value = "grey85") +
  scale_x_discrete(position = "bottom") + scale_y_discrete() + coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5), axis.title = element_blank(), panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "grey40")) +
  labs(title = "Correlação de Pearson entre Traits", subtitle = "Significância: *** p<0.001, ** p<0.01, * p<0.05")
print(p_cor)

# BOXPLOTS DE DISTRIBUIÇÃO ACROSS TRIALS
for (trait in traits_ok) {
  temp_data <- dat_kenya_filtrado[!is.na(dat_kenya_filtrado[[trait]]), ]
  if (nrow(temp_data) == 0) next
  
  env_order <- temp_data %>% count(ENV) %>% arrange(desc(n)) %>% pull(ENV)
  temp_data$ENV <- factor(temp_data$ENV, levels = env_order)
  nome_y <- trait_dict[[trait]]
  
  p <- ggplot(temp_data, aes(x = ENV, y = .data[[trait]], fill = ENV)) +
    geom_boxplot(na.rm = TRUE, outlier.size = 0.8, outlier.alpha = 0.5) +
    theme_minimal(base_size = 11) +
    labs(title = paste0("Distribuição de ", nome_y), x = "Ambiente", y = nome_y) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none", plot.title = element_text(face = "bold"))
  print(p)
}

# INTERAÇÃO GENÓTIPO X AMBIENTE
col_gen <- "GEN" 
top_genotipos <- dat_kenya_filtrado %>% filter(!is.na(GY)) %>% count(.data[[col_gen]]) %>% arrange(desc(n)) %>% slice(1:10) %>% pull(.data[[col_gen]])

gxe_data <- dat_kenya_filtrado %>%
  filter(.data[[col_gen]] %in% top_genotipos, !is.na(GY)) %>%
  group_by(ENV, .data[[col_gen]]) %>%
  summarise(Media_GY = mean(GY, na.rm = TRUE), .groups = "drop") %>%
  rename(Genotipo = .data[[col_gen]])

ordem_ambientes <- gxe_data %>% group_by(ENV) %>% summarise(Media_Env = mean(Media_GY, na.rm = TRUE)) %>% arrange(Media_Env) %>% pull(ENV)
gxe_data$ENV <- factor(gxe_data$ENV, levels = ordem_ambientes)

p_gxe <- ggplot(gxe_data, aes(x = ENV, y = Media_GY, group = Genotipo, color = Genotipo)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  theme_minimal(base_size = 11) +
  labs(title = "Interação GxE: Normas de Reação no Kenya", subtitle = "Produtividade (GY) dos 10 genótipos mais testados", x = "Ambientes (Menor para Maior Potencial)", y = "Produtividade Média (kg/ha)", color = "Genótipo") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
print(p_gxe)

# ==============================================================================
# 8. EXPORTAÇÃO DE RESULTADOS 
# ==============================================================================
# Gráfico de Qualidade Experimental (CV%)
ggsave("results/figures/01_qualidade_cv_kenya.png", plot = p_cv, 
       width = 10, height = 7, dpi = 300)
# Matriz de Correlação de Pearson
ggsave("results/figures/02_matriz_correlacao_kenya.png", plot = p_cor, 
       width = 9, height = 8, dpi = 300)

# Normas de Reação (Interação GxE)
ggsave("results/figures/03_normas_reacao_gxe.png", plot = p_gxe, 
       width = 11, height = 7, dpi = 300)

# Histogramas e Densidades (4 Trials selecionados)
ggsave("results/figures/04_histogramas_densidades_foco.png", plot = p_hist_dens, 
       width = 12, height = 8, dpi = 300)

# Estatísticas Descritivas
write.csv(desc_stats, "results/tables/01_estatisticas_descritivas_geral.csv", 
          row.names = FALSE)

# Relatório de Qualidade por Ambiente
write.csv(qc_stats, "results/tables/02_qualidade_ambientes_qc.csv", 
          row.names = FALSE)

# Estatísticas dos 4 Trials Selecionados
write.csv(stats_foco, "results/tables/03_estatisticas_foco_4trials.csv", 
          row.names = FALSE)

# Lista de Outliers Detectados
write.csv(outliers_detectados, "results/tables/04_outliers_identificados.csv", 
          row.names = FALSE)

