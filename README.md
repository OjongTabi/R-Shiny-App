# R Shiny App: Multi-function Data Analysis Tool

This repository contains an R Shiny application designed to perform a variety of data visualization and analysis tasks on user-provided data. This multi-function tool consists of four distinct modules, each performing different analyses and generating interactive visualizations.

## Modules

1. **Data Summary and Histogram Module (app1_server)**: This module allows users to upload a CSV file, displays a summary of the data, and plots a histogram for a selected numerical column. The summary includes information such as column names, types, mean, standard deviation, and distinct values.

2. **Data Filtering and Visualization Module (app2_server)**: This module allows the user to filter gene expression data based on variance and zero-count thresholds. It provides summary statistics of the filtered data and displays two scatter plots comparing various gene expression metrics and a heatmap of the filtered data. A Principal Component Analysis (PCA) scatter plot is also included.

3. **Volcano Plot Module (app3_server)**: This module enables users to create a volcano plot from gene expression data. The user can control the x and y-axes variables and set a p-value threshold for highlighting significant genes.

4. **Gene Expression Visualization Module (app4_server)**: This module takes count data and sample information as input, then generates a boxplot, violin plot, or bar plot of gene expression for a selected gene across different sample groups.

## Installation and Usage

### Prerequisites
- R version 4.0 or higher
- Shiny package installed in R
- Additional R packages: ggplot2, dplyr, DT, gridExtra, pheatmap, shinyWidgets, shinyFiles

### Steps

1. Clone this repository to your local machine.
2. Open RStudio and set your working directory to the location of the cloned repository.
3. Run the app using `shiny::runApp()` command.
4. Upload your data using the file upload buttons in the respective modules.
5. Adjust settings and parameters as per your requirements.


## License

This project is licensed under the MIT License 


Please feel free to raise any issues or to contribute to this project. If you need further assistance, contact the project maintainers.
