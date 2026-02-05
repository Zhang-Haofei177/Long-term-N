# Long-term-N
Overview
This repository contains all data and reproducible analysis code for the paper "Long-term nitrogen input reshapes symbiotic communities and drives functional trade-offs in legumes".

Project Structure
Fig.1_and_Fig.S4/ - Microbial diversity, community structure, and Mantel test

Fig.2_and_Fig.3/ - Community assembly processes of AMF and rhizobia

Fig.4/ - Network analysis and topological parameter line plots (network visualization requires edge and node files for Gephi)

Fig.5/ - Spearman correlation and LMM analysis of nitrogen cycling gene abundance in response to nitrogen application and host selection, plus PERMANOVA (nitrogen cycling figure generated from Niche_LMM_result.csv and Nitrogen_LMM_result.csv using Adobe Illustrator)

Fig.6/ - Effects of symbiotic systems on soil factors and SEM analysis (SEM figure based on SEM_summary.txt; model fit from SEM_fit_measures.txt)

Fig.S2/ - Community composition of bacteria, rhizobia, and AMF

Fig.S3/ - Rarefaction curves

Fig.S5/ - Source tracking analysis (figure requires AMF_source_track_result.csv and Rhizobia_source_track_result.csv for Adobe Illustrator)

Fig.S6/ - Neutral Community Model (NCM) analysis

Fig.S7/ - Boxplots of soil factors

Usage
Each folder contains:

Fully executable R scripts

Required data files

Direct script execution generates corresponding figures

Requirements
Ensure the following R packages are installed:

Core: tidyverse, vegan, ggplot2, ComplexHeatmap, circlize, linkET

Additional: patchwork, reshape2, dplyr, tidyr, ggsci

Run scripts from within each figure's directory after adjusting file paths if necessary.

Notes
Some visualizations require final editing in Adobe Illustrator for layout and aesthetics

File paths in scripts may need adjustment based on local directory structure

Network visualizations should be created using edge and node files in Gephi (v0.9.2 or later recommended)

License
MIT License

Contact
For questions or issues, please open an Issue in this repository or send to Zhanghaofei101@163.com.

