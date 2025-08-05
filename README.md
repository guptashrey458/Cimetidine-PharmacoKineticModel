# Cimetidine PBPK Model: Simulating Pharmacokinetics in Microgravity

## Project Overview

This project utilizes a Physiologically Based Pharmacokinetic (PBPK) model to simulate the behavior of the drug cimetidine in the human body. The primary goal is to predict how its concentration profile changes under microgravity conditions compared to a normal gravity baseline. This is achieved by integrating transcriptomic data (log-fold change values of gene expression) into the model's parameters.

The analysis workflow is as follows:

1. A baseline PBPK model is validated against observed experimental data.
2. Gene expression changes associated with microgravity are used to adjust key model parameters related to drug transport and metabolism.
3. Monte Carlo simulations are run for both baseline and microgravity conditions to generate 95% confidence intervals, representing population variability.
4. The results are visualized and compared to quantify the impact of microgravity on drug exposure (AUC) and peak concentration (Cmax).

## Key Features

- **Three-Compartment PBPK Model**: Simulates drug movement between Gut, Plasma, and Liver compartments using a system of ordinary differential equations (ODEs).
- **Transcriptomic Integration**: Modifies clearance parameters based on gene expression logFC values to simulate different physiological states.
- **Monte Carlo Analysis**: Assesses model uncertainty and predicts population variability to generate robust 95% confidence intervals.
- **Data-Driven Validation**: The baseline model is fitted and validated against real-world observed experimental data.
- **Comprehensive Visualization**: Generates clear plots comparing the baseline and microgravity scenarios, including changes in PK metrics and underlying parameters.

## Files in This Project

### `CimitadineModel.ipynb`
The main Jupyter Notebook that contains the complete, step-by-step analysis workflow. Run this file to reproduce the results.

### `CimetidineKineticModel.py`
A Python script containing the CimetidineKineticModel class. This class encapsulates the PBPK model, including its parameters, ODEs, and simulation functions.

### `cimetidine_obs.csv`
Input data file containing the observed experimental data (Time vs. Plasma Concentration) used to validate the baseline model.

### `CIMItadine (1).xlsx`
Input data file containing the gene expression data with logFC values for key drug transporters and metabolic enzymes.

### `logFC_average.xlsx`
Input data file containing the pre-averaged logFC values for key drug transporters and metabolic enzymes. These values are used to simulate the microgravity condition.

### `cimetidine_simcyp (1).csv`
Reference data file containing SimCYP simulation results for comparison.

### `README.md`
This documentation file.

## How to Use

### 1. Prerequisites
Ensure you have Python installed with the following libraries:
- pandas
- numpy
- scipy
- matplotlib
- openpyxl

You can install them using pip:
```bash
pip install pandas numpy scipy matplotlib openpyxl
```

### 2. Setup
Place all the project files (.ipynb, .py, .csv, .xlsx) in the same directory.

### 3. Running the Analysis
1. Open the `CimitadineModel.ipynb` file in a Jupyter Notebook environment.
2. Run the cells sequentially from top to bottom.
3. The notebook will load the data, run the simulations, and generate all the plots and summary tables automatically.

## Model Details

The PBPK model is based on three compartments. The effect of microgravity is simulated by adjusting the baseline clearance (CL) values of key transporters and enzymes using the following formula:

**CL_microgravity = CL_baseline × 2^logFC**

A positive logFC value indicates up-regulation (increased clearance), while a negative value indicates down-regulation (decreased clearance).

## Example Results

The primary output is a comparative plot showing the predicted plasma concentration curves for both conditions.

**Predicted Impact of Microgravity on Cimetidine PK**

This plot typically shows that microgravity leads to a higher peak concentration (Cmax) and greater overall drug exposure (AUC), indicating impaired drug clearance.

## Repository Structure

```
Cimetidine-PharmacoKineticModel/
├── README.md
├── CimitadineModel.ipynb          # Main analysis notebook
├── CimetidineKineticModel.py      # Core PBPK model class
├── cimetidine_obs.csv             # Observed experimental data
├── cimetidine_simcyp (1).csv      # SimCYP reference data
├── CIMItadine (1).xlsx            # Gene expression data
└── logFC_average.xlsx             # Processed logFC values
```

## Contributing

This project is part of ongoing research into pharmacokinetics under microgravity conditions. Contributions and suggestions are welcome.

## License

This project is available for academic and research purposes.