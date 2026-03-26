
# 🌱 Projeto de Análises de Dados - Kenya | LGN5822 - Genética-Estatística I: Fundamentos

## 🧬 Visão Geral
- Este projeto investiga a **Interação Genótipo x Ambiente (GxE)** em soja (*Glycine max*) no continente africano, com foco em dados da Quênia. O objetivo é realizar diversas análises com os dados disponíveis para formular hipóteses sobre diferentes caracteres de genótipos e ambientes testados.
- A análise utiliza dados históricos da rede **Pan-African Soybean Variety Trials (PATs)**, abrangendo avaliações fenotípicas, climáticas e edáficas coletadas ao longo de 10 anos.

## 📖 Referência Científica

Os dados e metodologias seguem as diretrizes estabelecidas no artigo:
> Araújo, M. S., et al. (2025). **High-resolution soybean trial data supporting the expansion of agriculture in Africa**. *Scientific Data*.

## 🗺️ Roadmap do Projeto

O desenvolvimento está estruturado em quatro fases evolutivas:

| Fase | Descrição | Status |
|------|-----------|:------:|
| **1. EDA** | Limpeza de dados, controle de qualidade (CV%) e estatística descritiva. | ✅ |
| **2. ??** | A definir | ⏳ |

## 🗂️ Estrutura do Repositório

Organização baseada em princípios de reprodutibilidade:

-   `data/`: Conjuntos de dados brutos e processados.
-   `scripts/`: Arquivos `.R` e `.qmd` do pipeline.
-   `results/`: Gráficos (boxplots, correlações) e tabelas resultantes.
-   `docs/`: Documentação técnica e referências em PDF.

## 🛠️ Tecnologias e Ferramentas
O ambiente de análise utiliza **Quarto v1.4+** e **R v4.5.1**.

## 🚀 Reprodutibilidade

Para reproduzir as análises:

1. Clone o repositório.
2. Certifique-se de que a estrutura de diretórios foi preservada.
3. No RStudio, abra o arquivo de script da fase desejada.
4. Renderize o documento via comando Quarto:
   ```bash
   quarto render scripts/01_eda.qmd
