# combiMS

Code and analysis results for the combiMS project.

For more information, please visit the [website](http://combims.eu/) of the project. This project was supported by the European Union Seventh Framework Programme (FP7/2007-2) under grant agreement No. 305397.

The analysis results compiled here are presented in the following publication:
> Bernardo-Faura, M. et al, Prediction of combination therapy based on perturbation modeling of the multiple sclerosis signaling network, submitted 2017.

## Workflow of the Project

1. Normalization of the raw data with [normalize_data.R](https://github.com/saezlab/combiMS/tree/master/code/normalization/normalize_data.R)
2. Patient-specific modeling on the Cluster
3. Model merging by subgroups
4. Analysis of model similarities after merging
5. Prediction of combination therapies


## License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/combiMS/blob/master/LICENSE.txt) or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).
