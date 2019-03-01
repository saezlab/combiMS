# combiMS

Data, code and analysis results of the CombiMS project.

For more information, please visit the [website](https://cordis.europa.eu/project/rcn/106479/factsheet/en) of the project. This project was supported by the European Union Seventh Framework Programme (FP7/2007-2) under grant agreement No. 305397 and the European Sys4MS project (Horizon2020: Eracosysmed: ID-43).

The analysis results compiled here are presented in the following publication:
Bernardo-Faura, M. et al, Prediction of combination therapies based on topological modeling of the immune signaling network in Multiple Sclerosis, [bioRxiv 541458](https://www.biorxiv.org/content/10.1101/541458v1) and under submission, 2019


## Workflow of the Project

1. Normalization of the raw data with [normalization_pipeline.R](https://github.com/saezlab/combiMS/blob/master/code/data_processing_and_normalization/normalization_pipeline.R)
2. Patient-specific modeling with CellNOptR, see [single_model_optimization](https://github.com/saezlab/combiMS/tree/master/code/single_model_optimization)
3. Model merging by subgroups, see [model_merging](https://github.com/saezlab/combiMS/tree/master/code/model_merging)
4. Analysis of model similarities after merging, see [similarity](https://github.com/saezlab/combiMS/tree/master/code/similarity)
5. Prediction of combination therapies, see [prediction_of_combination_therapies](https://github.com/saezlab/combiMS/tree/master/code/prediction_of_combination_therapies)


## License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/combiMS/blob/master/LICENSE.txt) or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).

## Requirements

The scripts collected in this repository are written in [R 3.4.0 (2017-04-21)](https://cran.r-project.org/) and were run in [RStudio Version 0.99.893](https://www.rstudio.com) and on the Cluster of the [Rheinisch-Westfälische Technische Hochschule Aachen](https://www.rwth-aachen.de/), which uses Centos 7.3 and the LSF job scheduling system version 9.1.3.0.

The [CellNOptR](https://saezlab.github.io/cellnopt/) packages to fit the signaling models for each individual patient was used in version 1.22.0.
