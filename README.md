# INDENG 242A Project: Time-Series Application in Stocks

## Project Description
This project focuses on applying machine learning techniques to analyze time-series data for portfolio risk management. Specifically, we analyze the S&P 500 and Nikkei 225 indices to estimate the Value-at-Risk (VaR) of an equally weighted portfolio at 99% and 95% confidence levels.

### Key Highlights:
- **Data**: Weekly adjusted closing prices of the S&P 500 and Nikkei 225 indices from January 1, 2000, to December 7, 2024, sourced from Yahoo Finance.
- **Models**: ARMA-GARCH and GARCH models were used to capture volatility clustering, while Copula-based Monte Carlo simulations were employed for joint modeling and risk estimation.
- **Results**: The project successfully estimated VaR, with a 99% VaR of 0.07469 and a 95% VaR of 0.03707 for the portfolio.

## Files Included
- `main.tex`: LaTeX source file for the report.
- `main.pdf`: Final compiled report.
- `SP500_N225.csv`: Data file containing weekly adjusted closing prices for the S&P 500 and Nikkei 225 indices.
- `r code`: containing R code for data processing, modeling, and simulations.
- `plots/`: Folder containing plots and figures used in the report.

## How to Run
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/your-repo-name.git
