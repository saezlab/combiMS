# combiMS

Code and analysis results for the combiMS project.

For more information, please visit the [website](http://combims.eu/) of the project. This project was supported by the European Union Seventh Framework Programme (FP7/2007-2) under grant agreement No. 305397.

The analysis results compiled here are presented in the following publication:
> Bernardo-Faura, M. et al, Prediction of combination therapy based on perturbation modeling of the multiple sclerosis signaling network, submitted 2017.

## Workflow of the Project

1. Normalization of the raw data with [normalization_pipeline.R](https://github.com/saezlab/combiMS/blob/master/code/data_processing_and_normalization/normalization_pipeline.R)
2. Patient-specific modeling with CellNOptR, see [single_model_optimization](https://github.com/saezlab/combiMS/tree/master/code/single_model_optimization)
3. Model merging by subgroups, see [model_merging](https://github.com/saezlab/combiMS/tree/master/code/model_merging)
4. Analysis of model similarities after merging, see [similarity](https://github.com/saezlab/combiMS/tree/master/code/similarity)
5. Prediction of combination therapies, see [prediction_of_combination_therapies](https://github.com/saezlab/combiMS/tree/master/code/prediction_of_combination_therapies)


## License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/combiMS/blob/master/LICENSE.txt) or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).
