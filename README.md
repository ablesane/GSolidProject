# GSolid Project - Solid Solubility of Oxygen in Metal Alloys

## Overview

The development of high entropy alloys (HEA) has expanded the field of materials research, providing enhanced material properties in structural metal alloys. HEAs are different combinations of elements with equi-atomic or nearly equi-atomic compositions to fulfill property requirements for different fields of applications. HEAs have unique behavior due to the high mixing entropy, sluggish diffusion, severe lattice distortion, and multi-metallic cocktail effects, also known as core effects. These core effects lead to the forming of a simple solid solution instead of complex intermetallic phases and help acquire a unique combination of properties, such as great strength, large ductility, corrosion resistance, and high thermal stability.

Oxygen, as an alloying element, plays a profound role in influencing metal properties and can significantly alter the material's characteristics. This understanding is not just crucial but also highly relevant for achieving the desired material properties in practical applications, underscoring the importance of our research on the role of oxygen in high entropy alloys.

The ability of a HEA to form and sustain a protective oxide layer is directly influenced by the concentration of dissolved oxygen on its surface. Therefore, understanding the thermodynamics of oxygen dissolution is critical in designing HEAs. However, the lack of available thermodynamic data has been a significant obstacle. By addressing this gap, this research provides invaluable insights for HEA design, significantly advancing our understanding of these alloys. 


## Table of Contents

- [GSolid Project - Solid Solubility of Oxygen in Metal Alloys](#gsolid-project---solid-solubility-of-oxygen-in-metal-alloys)
  - [Overview](#overview)
  - [Table of Contents](#table-of-contents)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
  - [Project Structure](#project-structure)
  - [Data Overview](#data-overview)
    - [Data Sources](#data-sources)
    - [Data Preparation](#data-preparation)
  - [Analysis](#analysis)
    - [Analytical Questions](#analytical-questions)
    - [Methods and Tools](#methods-and-tools)
  - [Results](#results)
  - [References & Acknowledgments](#references-&-acknowledgments)

## Getting Started

### Prerequisites

Ensure that you have the Python package dependencies installed before running this code. The required packages are listed in the `matprojenv.yml ` file of this repository.

### Project Structure

`data/`: Folder contains open source data sets from The Materials Project and Materials Platform for Data Science.

`notebooks/`: Jupyter notebooks with the data processing procedure.

`src/`: utility scripts used in this project.

`workflows/`: Machine Learning & Data Visualization workflows utilized in Orange Data Mining.

## Data Access

### Data Sources

#### Materials Platform for Data Science 
[The Materials Platform for Data Science](https://mpds.io/#start) (MPDS) is a digitalized version of the PAULING FILE materials database. The databases contains experimental observations of crystal structure, phase diagrams, and physical properties obtained from original scientific publications from 1891 to present. Data can be accessed via API and is returned in JSON format. This project will be extracting liquidus curve data from experimental phase diagrams of binary and ternary systems. Phase diagrams are predomiantely sourced from journal publications, with the majority sourced from the Russian Journal of Inorganic Chemistry, Journal of Alloys and Compounds, and CALPHAD. Open source data specified by the API license is released under the Creative Commons Attribution 4.0 International license. Commercial proprietary data is closed, confidential and not available for distribution.


#### The Materials Project
[The Materials Project](https://next-gen.materialsproject.org/) is a multi-institution, multi-national effort to compute the properties of all inorganic materials and provide the data and associated analysis algorithms for every materials researcher free of charge. The Materials Project is powered by the Python Materials Genomics package (pymatgen), an open source python library for materials analysis. Pymatgen has analytical tools for phase diagram generation, reaction balancing & calculation, diffusion analyses, and more. Pymatgen is integated with the Materials Project API to easily access materials data. Users of The Materials Project agree to accept the Creative Commons Attribution 4.0 License.

### Data Preparation
The data sourced from MPDS and The Materials Project will be processed to calculate the regular solution mixing enthalpy of oxygen in various metal oxide systems.
#### Hume-Rothery Rules
The Hume-Rothery rules are a helpful tool for predicting the formation of a solid solution. These rules specify that for a solid solution to form, the two elements involved must have similar atomic sizes, crystal structures, electronegativities, and valencies.

For a given material, atomic size, electronegativity,crystal structure, and valence are strongly related to one another. This project will only study the electronegativity of the metal oxide systems explored to supplement the regular solution mixing enthalpy calculations. 

## Analysis

### Analytical Questions
-	How can we characterize oxygen solubility limits in multi-metal alloy solid solutions?
-	How can we model a general thermodynamic framework that can characterize solubility limits in HEAs of various compositions?
-	What other material structure-property relationships can be utilized for systems with limited experimental data to improve the predictive model?

### Methods and Tools

The CALPHAD (Calculation of Phase Diagrams) method is the standard technique for modeling thermodynamic processes. The method transforms the available experimental data of material systems into physical-based mathematical models.3 Because it uses verified experimental data, the method can extrapolate material properties within the known domain with great accuracy but cannot scale to predict new alloy systems for which experimental data cannot be obtained.

Density functional theory (DFT) is a method employed to predict the properties of novel alloy systems. DFT calculations predict material behavior based on the behavior of materials at and below the atomic scale. It's worth noting that DFT-based calculations have lower thermodynamic accuracy compared to properties calculated by the CALPHAD method. 

Given the known limitations of CALPHAD and DFT-based calculations, the strategy is to develop a computational method that combines experimental phase boundaries and thermodynamic calculations from MPDS with density functional theory (DFT)-based solid-state energy calculations from Materials Project. It is through this approach that the oxidation resistance of novel structural alloys can be characterized.

#### Orange Data Mining
[Orange Data Mining](https://orangedatamining.com/) is an open source software tool used for creating machine learning pipelines and interactive data visualizations.  This project will use Orange to creat machine learning pipeline to assess and predict the calculated regular mixing enthalpy of various metal systems.

## References & Acknowledgement

1. Approach based on work conducted by S. Tan, J. Willwerth, A. Rauf & W. Sun, “Rapid estimation of the liquidus curve for unexplored alloy chemistries” [unpublished work] (2024)
2.	Poresh Kumar, Tu-Ngoc Lam, Pawan Kumar Tripathi, Sudhanshu Shekhar Singh, Peter K. Liaw, E-Wen Huang; Recent progress in oxidation behavior of high-entropy alloys: A review. APL Mater (December 2022)
3.	Domingo Jullian, Jianqiang Zhang, D. Brynn Hibbert, David J. Young, Oxygen solubility in austenitic Fe-Ni alloys at high temperatures, Journal of Alloys and Compounds, Volume 732, (2018)
4.	Ohtani, H. The CALPHAD Method. In: Czichos, H., Saito, T., Smith, L. (eds) Springer Handbook of Materials Measurement Methods. Springer Handbooks. Springer, Berlin, Heidelberg. (2006)

I would like to thank Dr. Wenhao Sun, Joshua Willwerth, and Shibo Tan and the rest of the Sun Research Group at the University of Michigan - Ann Arbor for their valuable guidance and support.

Project collaboration was made possible through the Michigan Insitiute For Data & AI In Society (MIDAS) Faculty Student Reasearch Connections Program.
