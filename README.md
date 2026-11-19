# Time-varying fiscal multipliers of social benefits in Brazil (1997–2019)

Replication code and data for the paper by **Marina da Silva Sanches, Daniel César Stumm, João Enrique Perez Diaz and Clara Zanon Brenck**.

The paper estimates time-varying fiscal multipliers for social benefits and total primary spending in Brazil with a kernel-smoothed time-varying VAR (TV-VAR) on monthly data from 1997 to 2019. It then relates the social benefits multiplier to inequality and poverty through correlation exercises and a threshold VAR (TVAR).

## Contents

| File | Description |
|---|---|
| `DynamicMultipliers.R` | Full script: estimation, robustness checks, diagnostics, tables and figures |
| `dadosBRFDesagreg.xlsx` | Monthly fiscal data (primary spending, revenue, social benefits) and deflator |
| `pib_mensal_4380.csv.csv` | Monthly GDP (Central Bank of Brazil, series 4380) |
| `Exogenas.xlsx` | Controls: real exchange rate and 12-month IPCA |
| `dadosBR_trimestral.xlsx` | Quarterly fiscal data and GDP |
| `GiniBR_trimestral.xlsx` | Gini index, quarterly (interpolated) |
| `Pobreza_trimestral.xlsx` | Poverty and extreme poverty, quarterly (interpolated) |
| `base_padronizada.xlsx` | Annual social variables used in the correlation analysis |

## What the script does

1. **Monthly TV-VAR** (Feb/1997–Dec/2019) for primary spending and social benefits: cumulative IRFs, dynamic multipliers (h = 12, G/Y as a 12-month moving average) and bootstrap confidence bands (68% and 90%).
2. **Robustness**: exogenous controls (real exchange rate and inflation), alternative bandwidths, alternative horizons and quarterly estimation (1997Q2–2019Q4).
3. **Residual diagnostics** of the auxiliary constant-coefficient VARs and of the TVARs.
4. **Quarterly TVAR** (1998Q1–2019Q4) for social benefits, using the Gini index, poverty and extreme poverty as threshold variables, with regime-specific multipliers and bootstrap bands.
5. **Correlation analysis** between the annual social benefits multiplier and Gini, poverty, extreme poverty, Top 1% and Top 10% (Pearson, with p-values).

Outputs written to the working folder: `var_residual_tests.csv`, `tvar_multipliers.csv`, `tvar_residual_tests.csv` and `correlacoes_denominador_variando.xlsx`. The figures are drawn in the R graphics device.

## How to run

1. Install R (≥ 4.0) and the packages:

```r
install.packages(c("tvReg", "vars", "readxl", "seasonal", "showtext",
                   "tseries", "writexl", "rstudioapi"))
install.packages("tsDyn", repos = "https://matthieustigler.r-universe.dev")
```

2. Download or clone this repository.
3. Open `DynamicMultipliers.R` in RStudio and run it. The working directory is set to the script's folder automatically. Outside RStudio, set it to the repository folder first with `setwd()`.

Notes:

- The bootstrap of the monthly IRFs (1,000 draws) is slow. The result is cached in `ic_bootstrap_68_90.rds` after the first run.
- The figures in the paper use the Mena Grotesk typeface, which is not distributed here. To use it, place its four `.otf` files in a `fonts/` folder. Otherwise the default sans font is used.

## Data sources

- **Primary spending and social benefits**: Sanches et al. (2026), built on Orair and Gobetti (2017). Social benefits include Bolsa Família, BPC, unemployment insurance, salary bonus (abono salarial) and social security benefits.
- **GDP**: Central Bank of Brazil, series 4380.
- **Deflator**: IPCA, Central Bank of Brazil, series 433.
- **Real exchange rate**: Central Bank of Brazil, series 11753.
- **Inflation (control)**: IPCA, 12-month accumulated change (Ipeadata, source IBGE).
- **Poverty and extreme poverty**: share of people below the PPP$8.30 and PPP$3.00 lines (Ipeadata).
- **Gini**: Gini coefficient of per capita household income from all sources (Ipeadata).
- **Top 1% and Top 10% income shares**: World Inequality Database (WID).

Fiscal series and GDP are deflated, seasonally adjusted with X-13ARIMA-SEATS and used in log first differences. For the TVAR, the annual distributive variables are converted to quarterly frequency with cubic spline interpolation.

## Contact

Marina da Silva Sanches (corresponding author), UNESP, Araraquara, Brazil: ms.sanches@unesp.br
