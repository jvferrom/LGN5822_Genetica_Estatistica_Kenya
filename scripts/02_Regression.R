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
dir_figuras <- "results/figures/02_regression/"
dir_tabelas <- "results/tables/02_regression/"

# Criar as pastas se elas não existirem
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

# Produtividade × Altura de planta (GY ~ PH_R8) ----------------------------------------------------------------------

#MANUAL
# Remover NAs para o par específico
df1 <- dat_clean %>%
  select(GY, PH_R8) %>%
  drop_na()

# Vetores y e X
y  <- df1$GY
X  <- cbind(1, df1$PH_R8)          # matriz modelo [n x 2]

# Estimativa de beta via OLS: beta = (X'X)^-1 X'y
XtX     <- crossprod(X)             # X'X
XtX_inv <- solve(XtX)               # (X'X)^-1
beta    <- XtX_inv %*% t(X) %*% y   # estimativas

cat("=== Coeficientes estimados manualmente ===\n")
cat(sprintf("  Intercepto (b0): %.4f\n", beta[1]))
cat(sprintf("  Inclinação  (b1): %.4f\n", beta[2]))

n  <- nrow(df1)                   # tamanho amostral
p  <- ncol(X)                     # número de parâmetros (b0 + b1)

# Valores ajustados e resíduos
y_hat  <- X %*% beta
e      <- y - y_hat

# Estimativa da variância dos erros (sigma^2)
SQE    <- as.numeric(crossprod(e))       # soma de quadrados dos erros
sigma2 <- SQE / (n - p)

# Desvios-padrão dos coeficientes
var_beta <- sigma2 * XtX_inv
SE_beta  <- sqrt(diag(var_beta))

# Estatísticas t e p-valores
t_stat <- beta / SE_beta
p_val  <- 2 * pt(-abs(t_stat), df = n - p)

# R² e R² ajustado
SQT    <- sum((y - mean(y))^2)
R2     <- 1 - SQE / SQT
R2adj <- 1 - (SQE / (n - p)) / (SQT / (n - 1))

# Estatística F global
SQReg  <- SQT - SQE
F_stat <- (SQReg / (p - 1)) / (SQE / (n - p))
pF     <- pf(F_stat, df1 = p - 1, df2 = n - p, lower.tail = FALSE)

cat("=== Inferência manual ===\n\n")
tabela_manual <- data.frame(
  Parametro  = c("Intercepto", "PH_R8"),
  Estimativa = round(beta, 4),
  EP         = round(SE_beta, 4),
  t_val      = round(t_stat, 4),
  p_valor    = format.pval(p_val, digits = 3, eps = 0.001)
)
print(tabela_manual, row.names = FALSE)

# Salvamento da tabela de inferência manual
write.csv(tabela_manual, file.path(dir_tabelas, "regressao_inferencia_gy_ph_r8.csv"), row.names = FALSE)

cat(sprintf("\nR²         : %.4f\n", R2))
cat(sprintf("R² ajustado: %.4f\n", R2adj))
cat(sprintf("F(%d, %d)   : %.2f  (p = %s)\n",
            p-1, n-p, F_stat, format.pval(pF, digits = 3)))
cat(sprintf("σ² estimado: %.2f   (σ = %.2f kg/ha)\n", sigma2, sqrt(sigma2)))

# Estimação via lm() e verificação
mod1 <- lm(GY ~ PH_R8, data = df1)
summary(mod1)
# Variáveis pré-calculadas para expressões inline
coef_mod1_b1     <- round(coef(mod1)[2], 2)
coef_mod1_b1_int <- round(coef(mod1)[2], 0)
r2_mod1          <- round(summary(mod1)$r.squared, 3)
r2_mod1_pct      <- round(summary(mod1)$r.squared * 100, 1)
# Verificação numérica: os valores devem ser idênticos
cat("=== Verificação: manual vs lm() ===\n")
cat(sprintf("b0 manual: %.6f  |  b0 lm: %.6f\n",
            beta[1], coef(mod1)[1]))
cat(sprintf("b1 manual: %.6f  |  b1 lm: %.6f\n",
            beta[2], coef(mod1)[2]))
cat(sprintf("R² manual: %.6f  |  R² lm: %.6f\n",
            R2, summary(mod1)$r.squared))

# Intervalos de confiança
# Manual: IC 95% = beta +/- t(alpha/2, n-p) * SE
alpha  <- 0.05
t_crit <- qt(1 - alpha/2, df = n - p)

IC_manual <- data.frame(
  Parametro  = c("Intercepto", "PH_R8"),
  Estimativa = round(as.numeric(beta), 3),
  LI_95      = round(as.numeric(beta) - t_crit * SE_beta, 3),
  LS_95      = round(as.numeric(beta) + t_crit * SE_beta, 3)
)

cat("=== IC 95% (manual) ===\n")
print(IC_manual, row.names = FALSE)

# Salvamento da tabela de Intervalo de Confiança
write.csv(IC_manual, file.path(dir_tabelas, "regressao_ic__gy_ph_r8.csv"), row.names = FALSE)

cat("\n=== IC 95% via confint() ===\n")
print(round(confint(mod1), 3))
#| fig-height: 5
#| fig-width: 10

# Formata o p-valor da estatística F global dinamicamente
p_str <- format.pval(pF, digits = 3, eps = 0.001)
p_label <- ifelse(grepl("<", p_str), paste0("p ", p_str), paste0("p = ", p_str))

# Gráfico 1: dispersão com linha de regressão
p_disp <- ggplot(df1, aes(x = PH_R8, y = GY)) +
  geom_point(alpha = 0.25, color = "#2C3E50", size = 1.5) +
  geom_smooth(method = "lm", color = "#E74C3C", fill = "#E74C3C",
              alpha = 0.15, linewidth = 1.2) +
  annotate("text", x = Inf, y = Inf,
           # Uso do %s para inserir a string p_label dinamicamente
           label = sprintf("R² = %.3f\ny = %.0f + %.1f·x\n%s", 
                           R2, beta[1], beta[2], p_label),
           hjust = 1.1, vjust = 1.5, size = 4,
           color = "#E74C3C", fontface = "bold") +
  labs(title    = "Produtividade × Altura de planta",
       subtitle = "Regressão linear simples",
       x        = "Altura de planta em R8 (cm)",
       y        = "Produtividade de grãos (kg/ha)") +
  theme_bw(base_size = 12)

# Gráfico 2: resíduos vs valores ajustados
# Utiliza y_hat e e (calculados matricialmente) em vez de mod1
df_res <- data.frame(ajustados = as.numeric(y_hat),
                     residuos  = as.numeric(e))

p_res <- ggplot(df_res, aes(x = ajustados, y = residuos)) +
  geom_point(alpha = 0.25, color = "#2980B9", size = 1.5) +
  geom_hline(yintercept = 0, color = "#E74C3C",
             linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "#27AE60",
              linewidth = 0.8, linetype = "dotted") +
  labs(title    = "Resíduos vs Valores ajustados",
       subtitle = "Diagnóstico de homocedasticidade",
       x        = "Valores ajustados (kg/ha)",
       y        = "Resíduos (kg/ha)") +
  theme_bw(base_size = 12)

# Unificando os gráficos para salvar
p_combined <- p_disp + p_res
print(p_combined)

# Salvamento dos gráficos combinados
ggsave(file.path(dir_figuras, "regressao_diagnostico__gy_ph_r8.png"), plot = p_combined, width = 10, height = 5, dpi = 300)

#| echo: false
p_mod1_b1   <- summary(mod1)$coefficients[2, 4]
p_mod1_str  <- format.pval(p_mod1_b1, digits = 3, eps = 0.001)
sig_mod1    <- ifelse(p_mod1_b1 < 0.05, "estatisticamente significativo",
                      "não estatisticamente significativo")
sinal_mod1  <- ifelse(coef(mod1)[2] >= 0, "positivo", "negativo")
n_mod1      <- length(mod1$residuals)


# Proteína × Óleo (PROT ~ OIL) ----------------------------------------------------------------------
#MANUAL
df2 <- dat_clean %>%
  select(PROT, OIL) %>%
  drop_na()

y2  <- df2$PROT
X2  <- cbind(1, df2$OIL)

XtX2     <- crossprod(X2)
XtX2_inv <- solve(XtX2)
beta2    <- XtX2_inv %*% t(X2) %*% y2

n2       <- nrow(df2)
p2       <- ncol(X2)
y2hat    <- X2 %*% beta2
e2       <- y2 - y2hat
SQE2     <- as.numeric(crossprod(e2))
sigma2_2 <- SQE2 / (n2 - p2)
SE2      <- sqrt(diag(sigma2_2 * XtX2_inv))
t2       <- beta2 / SE2
pv2      <- 2 * pt(-abs(t2), df = n2 - p2)
SQT2     <- sum((y2 - mean(y2))^2)
R2_2     <- 1 - SQE2 / SQT2
F2       <- ((SQT2 - SQE2) / (p2-1)) / (SQE2 / (n2 - p2))
pF2      <- pf(F2, df1 = p2-1, df2 = n2-p2, lower.tail = FALSE)

cat("=== PROT ~ OIL: Resultados manuais ===\n")
cat(sprintf("b0 (Intercepto): %.4f\n", beta2[1]))
cat(sprintf("b1 (OIL):        %.4f\n", beta2[2]))
cat(sprintf("R²:              %.4f\n", R2_2))
cat(sprintf("F(%d,%d):          %.2f  (p = %s)\n",
            p2-1, n2-p2, F2, format.pval(pF2, digits=3)))
cat(sprintf("σ estimado:      %.4f %%\n", sqrt(sigma2_2)))

# Criação do data.frame
tabela_manual_2 <- data.frame(
  Parametro  = c("Intercepto", "OIL"),
  Estimativa = round(as.numeric(beta2), 4),
  EP         = round(as.numeric(SE2), 4),
  t_val      = round(as.numeric(t2), 4),
  p_valor    = format.pval(as.numeric(pv2), digits = 3, eps = 0.001)
)

# Salvamento da tabela de inferência manual
write.csv(tabela_manual_2, file.path(dir_tabelas, "regressao_inferencia_prot_oil.csv"), row.names = FALSE)

# Verificação com lm()
mod2 <- lm(PROT ~ OIL, data = df2)
summary(mod2)
# Variáveis pré-calculadas para expressões inline
coef_mod2_b1     <- round(coef(mod2)[2], 3)
abs_coef_mod2_b1 <- abs(round(coef(mod2)[2], 2))
r2_mod2          <- round(summary(mod2)$r.squared, 3)
r2_mod2_pct      <- round(summary(mod2)$r.squared * 100, 1)
r2_mod2_2dp      <- round(summary(mod2)$r.squared, 2)

# Gráficos
#| fig-height: 5
#| fig-width: 10
p1 <- ggplot(df2, aes(x = OIL, y = PROT)) +
  geom_point(alpha = 0.4, color = "#8E44AD", size = 1.8) +
  geom_smooth(method = "lm", color = "#F39C12", fill = "#F39C12",
              alpha = 0.15, linewidth = 1.2) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\nb1 = %.3f\np = %.3f", R2_2, beta2[2], pF2),
           hjust = 1.1, vjust = 1.5, size = 4,
           color = "#F39C12", fontface = "bold") +
  labs(title    = "Proteína × Óleo (Kenya)",
       subtitle = "Relação não significativa no subconjunto NIR (n~170)",
       x        = "Teor de óleo (%)",
       y        = "Teor de proteína (%)") +
  theme_bw(base_size = 12)

# QQ-plot dos resíduos
df_res2 <- data.frame(residuos = mod2$residuals)
p2_qq <- ggplot(df_res2, aes(sample = residuos)) +
  stat_qq(alpha = 0.4, color = "#8E44AD") +
  stat_qq_line(color = "#E74C3C", linewidth = 1) +
  labs(title    = "QQ-plot dos resíduos",
       subtitle = "Avaliação da normalidade",
       x        = "Quantis teóricos",
       y        = "Quantis amostrais") +
  theme_bw(base_size = 12)

# Unificando os gráficos para visualização e salvamento
p_combined_2 <- p1 + p2_qq
print(p_combined_2)

# Salvamento
ggsave(file.path(dir_figuras, "regressao_diagnostico_prot_oil.png"), plot = p_combined_2, width = 10, height = 5, dpi = 300)

#| echo: false
p_mod2_b1  <- summary(mod2)$coefficients[2, 4]
p_mod2_str <- format.pval(p_mod2_b1, digits = 3, eps = 0.001)
sig_mod2   <- ifelse(p_mod2_b1 < 0.05, "estatisticamente significativo",
                     "não significativo (estatisticamente nulo)")
n_mod2     <- length(mod2$residuals)

# Dias para maturidade × Dias para florescimento (NDM ~ FLW_DAYS) ----------------------------------------------------------------------
#MANUAL
df3 <- dat_clean %>%
  select(NDM, FLW_DAYS) %>%
  drop_na()

y3       <- df3$NDM
X3       <- cbind(1, df3$FLW_DAYS)
beta3    <- solve(crossprod(X3)) %*% t(X3) %*% y3
n3       <- nrow(df3)
p3       <- 2
e3       <- y3 - X3 %*% beta3
SQE3     <- as.numeric(crossprod(e3))
sigma2_3 <- SQE3 / (n3 - p3)
SE3      <- sqrt(diag(sigma2_3 * solve(crossprod(X3))))
SQT3     <- sum((y3 - mean(y3))^2)
R2_3     <- 1 - SQE3 / SQT3

# Cálculos adicionais para inferência
t3       <- beta3 / SE3
pv3      <- 2 * pt(-abs(t3), df = n3 - p3)
F3       <- ((SQT3 - SQE3) / (p3 - 1)) / (SQE3 / (n3 - p3))
pF3      <- pf(F3, df1 = p3 - 1, df2 = n3 - p3, lower.tail = FALSE)

cat("=== NDM ~ FLW_DAYS: Resultados manuais ===\n")
cat(sprintf("b0 (Intercepto): %.4f\n", beta3[1]))
cat(sprintf("b1 (FLW_DAYS):   %.4f\n", beta3[2]))
cat(sprintf("R²:              %.4f\n", R2_3))
cat(sprintf("F(%d,%d):          %.2f  (p = %s)\n",
            p3-1, n3-p3, F3, format.pval(pF3, digits=3)))
cat(sprintf("σ estimado:      %.4f\n\n", sqrt(sigma2_3)))

# Tabela de inferência manual e salvamento
tabela_manual_3 <- data.frame(
  Parametro  = c("Intercepto", "FLW_DAYS"),
  Estimativa = round(as.numeric(beta3), 4),
  EP         = round(as.numeric(SE3), 4),
  t_val      = round(as.numeric(t3), 4),
  p_valor    = format.pval(as.numeric(pv3), digits = 3, eps = 0.001)
)

print(tabela_manual_3, row.names = FALSE)
write.csv(tabela_manual_3, file.path(dir_tabelas, "regressao_inferencia_ndm_flwdays.csv"), row.names = FALSE)

# Verificação lm()
mod3 <- lm(NDM ~ FLW_DAYS, data = df3)
summary(mod3)

cat("\n=== Verificação: manual vs lm() ===\n")
cat(sprintf("b0 manual: %.6f  |  b0 lm: %.6f\n", beta3[1], coef(mod3)[1]))
cat(sprintf("b1 manual: %.6f  |  b1 lm: %.6f\n", beta3[2], coef(mod3)[2]))
cat(sprintf("R² manual: %.6f  |  R² lm: %.6f\n", R2_3, summary(mod3)$r.squared))

# Gráfico
#| fig-height: 5
#| fig-width: 6

p3_disp <- ggplot(df3, aes(x = FLW_DAYS, y = NDM)) +
  geom_point(alpha = 0.2, color = "#16A085", size = 1.5) +
  geom_smooth(method = "lm", color = "#C0392B", fill = "#C0392B",
              alpha = 0.15, linewidth = 1.2) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\nb1 = %.3f\np = %s",
                           summary(mod3)$r.squared, coef(mod3)[2], format.pval(pF3, digits=3, eps=0.001)),
           hjust = 1.1, vjust = 1.5, size = 4.5,
           color = "#C0392B", fontface = "bold") +
  labs(title = "Dias para maturidade × Dias para florescimento",
       x     = "Dias para florescimento (FLW_DAYS)",
       y     = "Dias para maturidade (NDM)") +
  theme_bw(base_size = 13)

print(p3_disp)

ggsave(file.path(dir_figuras, "regressao_dispersao_ndm_flwdays.png"), plot = p3_disp, width = 6, height = 5, dpi = 300)

# Variáveis pré-calculadas para expressões inline
r2_mod3      <- round(summary(mod3)$r.squared, 3)
coef_mod3_b1 <- round(coef(mod3)[2], 2)



#Gráfico
#| fig-height: 5
#| fig-width: 6
mod3 <- lm(NDM ~ FLW_DAYS, data = df3)
cat("=== Resumo lm() ===\n")
cat(sprintf("b0: %.4f  |  b1: %.4f  |  R²: %.4f\n",
            coef(mod3)[1], coef(mod3)[2], summary(mod3)$r.squared))

ggplot(df3, aes(x = FLW_DAYS, y = NDM)) +
  geom_point(alpha = 0.2, color = "#16A085", size = 1.5) +
  geom_smooth(method = "lm", color = "#C0392B", fill = "#C0392B",
              alpha = 0.15, linewidth = 1.2) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\nb1 = %.3f",
                           summary(mod3)$r.squared, coef(mod3)[2]),
           hjust = 1.1, vjust = 1.5, size = 4.5,
           color = "#C0392B", fontface = "bold") +
  labs(title = "Dias para maturidade × Dias para florescimento",
       x     = "Dias para florescimento (FLW_DAYS)",
       y     = "Dias para maturidade (NDM)") +
  theme_bw(base_size = 13)
# Variáveis pré-calculadas para expressões inline
r2_mod3      <- round(summary(mod3)$r.squared, 3)
coef_mod3_b1 <- round(coef(mod3)[2], 2)

# ==============================================================================
# 4. Modelos Polinomiais
# ==============================================================================
# ------------------------------------------------------------------------------
# Modelo polinomial: NDM ~ ELEV + ELEV²
# ------------------------------------------------------------------------------

# Implementação manual
df_ne_poly <- dat_clean %>%
  select(NDM, ELEV) %>%
  drop_na() %>%
  mutate(
    ELEV_c  = ELEV - mean(ELEV),
    ELEV_c2 = ELEV_c^2
  )
mean_elev_ne <- mean(df_ne_poly$ELEV)

y_np   <- df_ne_poly$NDM
X_np   <- as.matrix(cbind(1, df_ne_poly$ELEV_c, df_ne_poly$ELEV_c2))
beta_np <- solve(crossprod(X_np)) %*% t(X_np) %*% y_np
n_np   <- nrow(df_ne_poly); p_np <- ncol(X_np)
e_np   <- y_np - X_np %*% beta_np
SQE_np <- as.numeric(crossprod(e_np))
SQT_np <- sum((y_np - mean(y_np))^2)
R2_np  <- 1 - SQE_np / SQT_np
R2adj_np <- 1 - (SQE_np/(n_np-p_np)) / (SQT_np/(n_np-1))
sig2_np <- SQE_np / (n_np - p_np)
SE_np   <- sqrt(diag(sig2_np * solve(crossprod(X_np))))
t_np    <- as.numeric(beta_np) / SE_np
pv_np   <- 2 * pt(-abs(t_np), df = n_np - p_np)

cat("=== NDM ~ ELEV_c + ELEV_c² (manual, ELEV centrada) ===\n\n")
tabela_manual_poly <- data.frame(
  Param      = c("b0 (Intercepto)", "b1 (ELEV_c)", "b2 (ELEV_c²)"),
  Estimativa = round(as.numeric(beta_np), 6),
  EP         = round(SE_np, 6),
  t          = round(t_np, 3),
  p_valor    = format.pval(pv_np, digits = 3, eps = 0.001)
)
print(tabela_manual_poly, row.names = FALSE)
cat(sprintf("\nR²: %.4f  |  R²adj: %.4f\n", R2_np, R2adj_np))

vertice_np_c <- -beta_np[2] / (2 * beta_np[3])
vertice_np   <- vertice_np_c + mean_elev_ne
cat(sprintf("Vértice (escala original): %.0f m\n", vertice_np))

# Salvamento
write.csv(tabela_manual_poly, file.path(dir_tabelas, "regressao_inferencia_ndm_elev_quad.csv"), row.names = FALSE)

# Verificação com lm()
mod_ne_lin  <- lm(NDM ~ ELEV_c,               data = df_ne_poly)
mod_ne_quad <- lm(NDM ~ ELEV_c + I(ELEV_c^2), data = df_ne_poly)

cat("\n=== Comparação modelos ===\n")
cat("R² linear:    ", round(summary(mod_ne_lin)$r.squared,  4), "\n")
cat("R² quadrático:", round(summary(mod_ne_quad)$r.squared, 4), "\n\n")
print(anova(mod_ne_lin, mod_ne_quad))

# Gráficos
new_ne <- data.frame(
  ELEV_c = seq(min(df_ne_poly$ELEV_c), max(df_ne_poly$ELEV_c), length = 300)
) %>% mutate(ELEV = ELEV_c + mean_elev_ne)
new_ne$pred_lin  <- predict(mod_ne_lin,  newdata = new_ne)
new_ne$pred_quad <- predict(mod_ne_quad, newdata = new_ne)

p_ndm_elev <- ggplot(df_ne_poly, aes(x = ELEV, y = NDM)) +
  geom_point(alpha = 0.18, color = "#7F8C8D", size = 1.5) +
  geom_line(data = new_ne, aes(y = pred_lin,  color = "Linear"),     linewidth = 1.2) +
  geom_line(data = new_ne, aes(y = pred_quad, color = "Quadrático"), linewidth = 1.2) +
  scale_color_manual(values = c("Linear" = "#3498DB", "Quadrático" = "#E74C3C")) +
  labs(title    = "Dias para maturidade × Altitude (Kenya)",
       subtitle = sprintf("Vértice: ~%.0f m  |  R²_lin=%.3f  vs  R²_quad=%.3f",
                          vertice_np,
                          summary(mod_ne_lin)$r.squared,
                          summary(mod_ne_quad)$r.squared),
       x = "Altitude (m)", y = "NDM (dias)", color = "Modelo") +
  theme_bw(base_size = 12) + theme(legend.position = "top")

df_ne_faixa <- df_ne_poly %>%
  mutate(faixa = cut(ELEV,
                     breaks = c(0, 500, 1000, 1500, 2000, Inf),
                     labels = c("0–500","500–1000","1000–1500","1500–2000",">2000")))

p_ndm_elev_box <- ggplot(df_ne_faixa, aes(x = faixa, y = NDM, fill = faixa)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.15) +
  scale_fill_brewer(palette = "Blues") +
  labs(title = "NDM por faixa de altitude",
       x = "Faixa (m)", y = "NDM (dias)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_ndm_elev_combined <- p_ndm_elev + p_ndm_elev_box
print(p_ndm_elev_combined)

ggsave(file.path(dir_figuras, "regressao_quad_ndm_elev.png"), 
       plot = p_ndm_elev_combined, width = 10, height = 5, dpi = 300)


# ==============================================================================
# 5. Modelos com Interações
# ==============================================================================

# ------------------------------------------------------------------------------
# Interação: NDM ~ FLW_DAYS * ELEV
# ------------------------------------------------------------------------------

df_int <- dat_clean %>%
  select(NDM, FLW_DAYS, ELEV) %>%
  drop_na() %>%
  mutate(
    ELEV_c = ELEV - mean(ELEV),
    inter  = FLW_DAYS * ELEV_c
  )

# Implementação manual
y_int <- df_int$NDM
X_int <- as.matrix(cbind(1, df_int$FLW_DAYS, df_int$ELEV_c, df_int$inter))

beta_int  <- solve(crossprod(X_int)) %*% t(X_int) %*% y_int
n_int     <- nrow(df_int); p_int <- ncol(X_int)
e_int     <- y_int - X_int %*% beta_int
SQE_int   <- as.numeric(crossprod(e_int))
SQT_int   <- sum((y_int - mean(y_int))^2)
R2_int    <- 1 - SQE_int / SQT_int
R2adj_int <- 1 - (SQE_int/(n_int-p_int)) / (SQT_int/(n_int-1))
sig2_int  <- SQE_int / (n_int - p_int)
SE_int    <- sqrt(diag(sig2_int * solve(crossprod(X_int))))
t_int     <- as.numeric(beta_int) / SE_int
pv_int    <- 2 * pt(-abs(t_int), df = n_int - p_int)

cat("=== NDM ~ FLW_DAYS * ELEV_c (manual) ===\n\n")
tabela_int <- data.frame(
  Param = c("b0", "b1 (FLW_DAYS)", "b2 (ELEV_c)", "b3 (FLW_DAYS:ELEV_c)"),
  Est   = round(as.numeric(beta_int), 6),
  EP    = round(SE_int, 6),
  t     = round(t_int, 3),
  p     = format.pval(pv_int, digits = 3, eps = 0.001)
)
print(tabela_int, row.names = FALSE)
cat(sprintf("\nR²: %.4f  |  R²adj: %.4f\n", R2_int, R2adj_int))

# Salvamento
write.csv(tabela_int, file.path(dir_tabelas, "regressao_interacao_ndm_flwdays_elev.csv"), row.names = FALSE)

# Verificação via lm()
mod_adit <- lm(NDM ~ FLW_DAYS + ELEV_c, data = df_int)
mod_int  <- lm(NDM ~ FLW_DAYS * ELEV_c, data = df_int)

cat("\n=== Modelo aditivo ===\n")
cat(sprintf("R² = %.4f  |  R²adj = %.4f\n",
            summary(mod_adit)$r.squared, summary(mod_adit)$adj.r.squared))
cat("\n=== Modelo com interação ===\n")
print(summary(mod_int))
cat("\n=== Teste F: aditivo vs interação ===\n")
print(anova(mod_adit, mod_int))

# Visualização da interação
mean_elev_int <- mean(dat_clean$ELEV, na.rm = TRUE)
q_elev <- quantile(df_int$ELEV, c(0.10, 0.50, 0.90), na.rm = TRUE)
q_elev_c <- q_elev - mean_elev_int

pred_df <- expand.grid(
  FLW_DAYS = seq(min(df_int$FLW_DAYS), max(df_int$FLW_DAYS), length = 100),
  ELEV_c   = q_elev_c
) %>%
  mutate(
    NDM_pred = predict(mod_int, newdata = .),
    ELEV_lbl = factor(
      paste0("ELEV = ", round(ELEV_c + mean_elev_int, 0), " m"),
      levels = paste0("ELEV = ", round(q_elev, 0), " m"))
  )

p_int1 <- ggplot(pred_df, aes(x = FLW_DAYS, y = NDM_pred, color = ELEV_lbl)) +
  geom_line(linewidth = 1.3) +
  scale_color_brewer(palette = "Dark2") +
  labs(title    = "NDM ~ FLW_DAYS: efeito moderado pela altitude",
       subtitle = "Linhas de regressão para P10, mediana e P90 de ELEV",
       x = "Dias para florescimento (FLW_DAYS)", y = "NDM previsto (dias)",
       color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")

p_int2 <- df_int %>%
  mutate(elev_grp = cut(ELEV,
                        breaks = quantile(ELEV, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                        labels = c("Baixa (<P33)", "Média (P33–P67)", "Alta (>P67)"),
                        include.lowest = TRUE)) %>%
  ggplot(aes(x = FLW_DAYS, y = NDM, color = elev_grp)) +
  geom_point(alpha = 0.15, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("#3498DB","#E67E22","#C0392B")) +
  labs(title    = "Regressões por grupo de altitude",
       subtitle = "Diagnóstico visual da interação FLW_DAYS × ELEV",
       x = "FLW_DAYS (dias)", y = "NDM (dias)", color = "Altitude") +
  theme_bw(base_size = 12) + theme(legend.position = "top")

p_int_combined <- p_int1 + p_int2
print(p_int_combined)

ggsave(file.path(dir_figuras, "regressao_interacao_ndm_flwdays_elev.png"), 
       plot = p_int_combined, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# Interação: GY ~ NDM * ELEV
# ------------------------------------------------------------------------------

df_gy_int <- dat_clean %>% select(GY, NDM, ELEV) %>% drop_na()

mod_gy_adit <- lm(GY ~ NDM + ELEV, data = df_gy_int)
mod_gy_int  <- lm(GY ~ NDM * ELEV, data = df_gy_int)

cat("\n=== GY ~ NDM * ELEV ===\n")
cat("=== Aditivo: R²adj =", round(summary(mod_gy_adit)$adj.r.squared, 4), "\n")
cat("=== Interação: R²adj =", round(summary(mod_gy_int)$adj.r.squared, 4), "\n\n")
print(anova(mod_gy_adit, mod_gy_int))
print(summary(mod_gy_int))

# Salvamento
tabela_gy_int <- data.frame(
  Parametro = names(coef(mod_gy_int)),
  Estimativa = round(coef(mod_gy_int), 6),
  EP = round(summary(mod_gy_int)$coefficients[,2], 6),
  t_valor = round(summary(mod_gy_int)$coefficients[,3], 4),
  p_valor = format.pval(summary(mod_gy_int)$coefficients[,4], digits = 3, eps = 0.001)
)
write.csv(tabela_gy_int, file.path(dir_tabelas, "regressao_interacao_gy_ndm_elev.csv"), row.names = FALSE)

# Ponto de indiferença
b_NDM <- coef(mod_gy_int)["NDM"]
b_int_gy <- coef(mod_gy_int)["NDM:ELEV"]
ponto_indiferenca <- round(-b_NDM / b_int_gy, 0)

cat(sprintf("\nPonto de indiferença (altitude onde efeito de NDM sobre GY se anula): %d m\n", ponto_indiferenca))


# ==============================================================================
# 6. Diagnóstico Completo dos Modelos
# ==============================================================================

# Selecionar um modelo múltiplo para diagnóstico (GY ~ PH_R8 + W100G + NDM)
df_mult <- dat_clean %>%
  select(GY, PH_R8, W100G, NDM) %>%
  drop_na()

mod_mult <- lm(GY ~ PH_R8 + W100G + NDM, data = df_mult)

# Gráficos de diagnóstico
col_escuro <- adjustcolor("#2C3E50", alpha.f = 0.3)
col_verm   <- adjustcolor("#C0392B", alpha.f = 0.5)

png(file.path(dir_figuras, "diagnostico_modelo_multiplo.png"), 
    width = 11, height = 10, units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

plot(mod_mult, which = 1, main = "1. Resíduos vs Ajustados",
     col = col_escuro, pch = 16)
plot(mod_mult, which = 2, main = "2. Normal QQ-plot",
     col = col_escuro, pch = 16)
plot(mod_mult, which = 3, main = "3. Scale-Location",
     col = col_escuro, pch = 16)
plot(mod_mult, which = 4, main = "4. Distância de Cook",
     col = col_verm,   pch = 16)
plot(mod_mult, which = 5, main = "5. Resíduos vs Leverage",
     col = col_escuro, pch = 16)

hist(mod_mult$residuals, breaks = 60, col = "#3498DB",
     border = "white", main = "6. Distribuição dos resíduos",
     xlab = "Resíduo (kg/ha)", ylab = "Frequência")
curve(dnorm(x, mean = 0, sd = sd(mod_mult$residuals)) *
        length(mod_mult$residuals) *
        diff(range(mod_mult$residuals)) / 60,
      add = TRUE, col = "#E74C3C", lwd = 2)
legend("topright", "Normal teórica", col = "#E74C3C", lwd = 2, bty = "n")

dev.off()

# Testes de homocedasticidade
bp_mult  <- lmtest::bptest(mod_mult)
ncv_mult <- car::ncvTest(mod_mult)
bp_ph    <- lmtest::bptest(mod1)

cat("\n=== Testes de Homocedasticidade ===\n")
cat(sprintf("Breusch-Pagan (GY ~ PH_R8+W100G+NDM): BP = %.3f, df = %d, p = %s\n",
            bp_mult$statistic, bp_mult$parameter,
            format.pval(bp_mult$p.value, digits = 3, eps = 0.001)))
cat(sprintf("Breusch-Pagan (GY ~ PH_R8): BP = %.3f, df = %d, p = %s\n",
            bp_ph$statistic, bp_ph$parameter,
            format.pval(bp_ph$p.value, digits = 3, eps = 0.001)))
cat(sprintf("NCV (car), modelo múltiplo: Chi² = %.3f, df = %d, p = %s\n",
            ncv_mult$ChiSquare, ncv_mult$Df,
            format.pval(ncv_mult$p, digits = 3, eps = 0.001)))

# Salvamento dos testes
testes_homoc <- data.frame(
  Teste = c("Breusch-Pagan (múltiplo)", "Breusch-Pagan (simples)", "NCV"),
  Estatistica = c(round(bp_mult$statistic, 3), round(bp_ph$statistic, 3), round(ncv_mult$ChiSquare, 3)),
  gl = c(bp_mult$parameter, bp_ph$parameter, ncv_mult$Df),
  p_valor = c(format.pval(bp_mult$p.value, 3, eps=0.001),
              format.pval(bp_ph$p.value, 3, eps=0.001),
              format.pval(ncv_mult$p, 3, eps=0.001))
)
write.csv(testes_homoc, file.path(dir_tabelas, "testes_homocedasticidade.csv"), row.names = FALSE)

# Teste de normalidade (Shapiro-Wilk)
set.seed(42)
res_sample <- sample(mod_mult$residuals, min(5000, length(mod_mult$residuals)))
sw <- shapiro.test(res_sample)
cat(sprintf("\n=== Shapiro-Wilk (normalidade dos resíduos) ===\n"))
cat(sprintf("W = %.4f, p-valor = %s\n", sw$statistic, format.pval(sw$p.value, digits=3)))
cat(sprintf("n da amostra: %d\n", length(res_sample)))

# Salvamento
shapiro_result <- data.frame(
  Teste = "Shapiro-Wilk",
  W = round(sw$statistic, 4),
  p_valor = format.pval(sw$p.value, digits = 3, eps = 0.001),
  n_amostra = length(res_sample)
)
write.csv(shapiro_result, file.path(dir_tabelas, "teste_normalidade_residuos.csv"), row.names = FALSE)


# ==============================================================================
# 7. Tabela Resumo dos Modelos
# ==============================================================================

# Modelos adicionais para a tabela resumo
mod_elev <- lm(GY ~ ELEV, data = dat_clean)
mod_fe <- lm(FLW_DAYS ~ ELEV, data = dat_clean)
mod_prot_elev <- lm(PROT ~ ELEV, data = dat_clean)
mod_oil_elev <- lm(OIL ~ ELEV, data = dat_clean)
mod_lat <- lm(NDM ~ LAT, data = dat_clean)
mod_quad <- lm(GY ~ ELEV + I(ELEV^2), data = dat_clean)

modelos_lista <- list(
  "GY ~ PH_R8"              = mod1,
  "PROT ~ OIL"              = mod2,
  "NDM ~ FLW_DAYS"          = mod3,
  "GY ~ PH_R8+W100G+NDM"   = mod_mult,
  "GY ~ ELEV"               = mod_elev,
  "FLW_DAYS ~ ELEV"         = mod_fe,
  "PROT ~ ELEV"             = mod_prot_elev,
  "OIL ~ ELEV"              = mod_oil_elev,
  "GY ~ ELEV + ELEV²"      = mod_quad,
  "NDM ~ ELEV + ELEV²"     = mod_ne_quad,
  "NDM ~ LAT"               = mod_lat,
  "NDM ~ FLW_DAYS*ELEV"    = mod_int,
  "GY ~ NDM*ELEV"           = mod_gy_int
)

diag_tab <- map_dfr(modelos_lista, function(m) {
  s <- summary(m)
  tibble(
    n     = length(m$residuals),
    R2    = round(s$r.squared,     3),
    R2adj = round(s$adj.r.squared, 3),
    RMSE  = round(sqrt(mean(m$residuals^2)), 1),
    F_p   = format.pval(pf(s$fstatistic[1], s$fstatistic[2],
                           s$fstatistic[3], lower.tail=FALSE),
                        digits=2, eps=0.001)
  )
}, .id = "Modelo")

write.csv(diag_tab, file.path(dir_tabelas, "tabela_resumo_modelos.csv"), row.names = FALSE)

print(diag_tab)