Parental Emotional Expressivity and Children’s Emotional Development

This repository contains the analysis code associated with the manuscript:

Who Influences Children’s Emotional Development? Longitudinal Associations of Parental Emotional Expressivity with Children’s Emotional Intelligence and Depressive Symptoms

Overview

The study used a three-wave longitudinal design to examine associations among parental emotional expressivity, children’s emotional intelligence, and depressive symptoms.

Two complementary analytic approaches were used:

Latent profile analysis (LPA) to identify distinct profiles of paternal and maternal positive and negative emotional expressivity and to compare children’s emotional intelligence and depressive symptoms across these profiles.
Cross-lagged panel network (CLPN) analysis to examine cross-temporal predictive associations among specific dimensions of parental emotional expressivity, emotional intelligence, and depressive symptoms across T1–T2 and T2–T3.

Additional analyses included longitudinal measurement invariance, multiple imputation, network accuracy and stability analyses, node strength estimation, and cross-network similarity analyses.

Software

Analyses were conducted using:

SPSS 26.0
Mplus 8.3
R 4.5.2

The R analyses used packages including:

mice
bootnet
qgraph

Additional package dependencies are specified in the analysis scripts where applicable.

Repository Structure
R/

Contains R scripts for data preparation, multiple imputation, cross-lagged panel network estimation, network stability analyses, network visualization, and network similarity analyses.

Mplus/

Contains Mplus input files for longitudinal measurement invariance testing, latent profile analysis, and BCH distal outcome comparisons.

output/

Contains non-sensitive analysis outputs that can be shared publicly.

Analysis Workflow

The analyses were conducted in the following general order:

Data preparation and descriptive analyses
Longitudinal measurement invariance testing
Latent profile analysis
BCH comparisons of distal outcomes
Multiple imputation for the longitudinal network analyses
T1–T2 cross-lagged panel network estimation
T2–T3 cross-lagged panel network estimation
Network accuracy and centrality stability analyses
Network similarity analyses
Data Availability

The participant-level data are not publicly available at this time because the study involves longitudinal data collected from children and includes potentially sensitive psychological information.

The analysis code used to produce the main results reported in the manuscript is publicly available in this repository.

Because the raw data are not included, the analysis scripts cannot be reproduced directly without access to the original dataset.

Citation

If you use or refer to this code, please cite the associated article once it is published.

Manuscript:
Who Influences Children’s Emotional Development? Longitudinal Associations of Parental Emotional Expressivity with Children’s Emotional Intelligence and Depressive Symptoms

Contact

For questions regarding the analysis code or research materials, please contact the corresponding author.
