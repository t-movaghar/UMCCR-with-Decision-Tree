# Integrating Decision Tree Learning into the Unified Mixture Cure Model

## Project Overview

This repository contains my attempt to modify Dr. Suvra Pal's Unified Mixture Cure Model (UMCM) to integrate the predictive capabilities of decision trees in determining the cure rate. The UMCM is a mixture cure model designed to handle multiple competing risks, providing a robust framework for survival analysis in biomedical and clinical research.

My approach focuses on leveraging decision trees to estimate the cure fraction, aiming to improve the model's adaptability to complex, high-dimensional datasets. By incorporating machine learning techniques, I seek to enhance the model’s predictive accuracy and interpretability when analyzing survival outcomes.

## Motivation

Survival analysis in biomedical research often involves patients who are "cured," meaning they are no longer at risk of the event of interest. Traditional parametric approaches may not fully capture the underlying complexities in determining cure rates, particularly when multiple competing risks are involved. Decision trees provide a data-driven approach to identifying patterns in survival data, allowing for improved estimation of the cure probability.

## Key Components

Preprocessing Pipeline: Standardized methods for handling survival data, including feature selection and transformation.

Data Exploration: Visualization of key dataset characteristics, ensuring meaningful feature extraction.

Model Integration: Modification of the UMCM to include decision tree-based predictions for the cure rate parameter.

Performance Evaluation: Assessment of model accuracy and robustness using real and simulated datasets.

## Methodology

1. Exploratory Data Analysis

Visualization of survival distributions and covariate relationships.

Initial estimations of cure rates using traditional methods

1. Data Preprocessing

Transforming features for optimal model performance.

3. Integration of Decision Trees in the M-Step

Utilizing regression trees to model the cure rate within the EM algorithm framework.

Comparing traditional and tree-based estimation approaches.

## Future Directions

Extending the model to incorporate ensemble methods (e.g., Random Forest, XGBoost) for further accuracy improvements.

Applying the approach to real-world clinical datasets with high-dimensional covariates.

Refining the decision tree integration to optimize interpretability and generalization.


This repository is a work in progress, and I welcome any feedback or suggestions. Feel free to open issues or contribute via pull requests!
