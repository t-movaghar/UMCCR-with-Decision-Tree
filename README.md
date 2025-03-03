# Integrating Decision Tree Learning into the Unified Mixture Cure Model

## Project Overview

This repository contains my attempt to modify Dr. Suvra Pal's Unified Mixture Cure Model (UMCCR) to integrate the predictive capabilities of decision trees in determining the cure rate. The UMCCR is a mixture cure model designed to handle multiple competing risks, providing a new framework for survival analysis in biomedical and clinical research.

My approach focuses on using decision trees to estimate the cure fraction, with the aim to improve the model's adaptability to complex, high-dimensional datasets. I seek to enhance the model’s predictive accuracy and interpretability when analyzing survival outcomes.

## Motivation

Survival analysis in biomedical research often involves patients who are "cured," meaning they are no longer at risk of the event of interest. Traditional parametric approaches may not fully capture the underlying complexities in determining cure rates, particularly when multiple competing risks are involved. Decision trees provide a data-driven approach to identifying patterns in survival data, allowing for improved estimation of the cure probability.

## Key Components

* Data Exploration: Visualization of key dataset characteristics, ensuring meaningful feature extraction.

* Model Integration: Modification of the UMCM to include decision tree-based predictions for the cure rate parameter.

* Performance Evaluation: Assessment of model accuracy and robustness using real and simulated datasets.

## Methodology

1. Exploratory Data Analysis

Data Visualizations:

![image](https://github.com/user-attachments/assets/071c8743-034a-4b9f-9c8c-d03bca013df8)
"Alive" patients are right-censored. The probability of belonging to a cured or uncured group must be estimated.

Initial estimations of cure rates using traditional methods:

Kaplan-Meier estimations by patient race:
![image](https://github.com/user-attachments/assets/ce4d284f-61e4-463d-8f8a-77bba4f347a6)

Kaplan-Meier survival estimations by cancer type:
![image](https://github.com/user-attachments/assets/f7bfdbf0-9ed6-41bb-b6a0-65a0c6739c29)

Survival probabilities estimated by the kaplan-meier curves do not appear to be approaching 0 as time increases, indicating the presence of a cured population.
Results of the KM curves at time = 90 months will serve as initial estimates for cure rate for the respective patient groups.

![image](https://github.com/user-attachments/assets/23ca3564-f6e7-4a93-b85c-64ebeffc184a)
Initial guess for the overall cure rate: 0.28904936549925275


2. Data
  * Type: TXT, CSV
    * Input: One TXT file (Brstdata.txt) containing information on patient race, cancer type, months survived, cause of death, and right censorship indicatior (D)
    * Input (Visualizations): One CSV file (Brstdata_Distant_5064.csv) containing the above information
  * Target Variable: Cause of Death = Breast Cancer
  * 
The data used in this project is a subset of the SEER breast cancer dataset:
![image](https://github.com/user-attachments/assets/fae46cc2-8acb-4319-a8a7-494cf1159d81)

Patient race = white if race_black and race_asian = 0
Cancer type = Luminal_A if all other cancer types = 0


3. Integration of Decision Trees in the M-Step (WIP)

Tree to model the cure rate within the EM algorithm framework.
![image](https://github.com/user-attachments/assets/51166f6f-1120-46bf-bc45-e0ddb6b632aa)
Using estimated cured probability in calculating the weight for the cured group.
Iterates until convergence.

Comparing traditional and tree-based estimation approaches. (WIP)

## Future Directions

Extending the model to incorporate ensemble methods (e.g., Random Forest, XGBoost) for further accuracy improvements.

Applying the approach to real-world clinical datasets with high-dimensional covariates.

Refining the decision tree integration to optimize interpretability and generalization.


This repository is a work in progress, and I welcome any feedback or suggestions. Feel free to open issues or contribute via pull requests!
