# Nicholas_Cowie_Masters_Thesis
This repository stores the data and scripts used to generate the results to this thesis

## Instructions
By running these instructions you should be able to:
1. Sample the model used in this thesis
2. Reproduce the plots used in the thesis

### Installing Maud
The first step requires installing maud in order to conduct the Bayesian sampling.
This process is simplified with virtualenvironments and a bash environment.

```
pip install virtualenv
```
Go to the directory Maud in Nicholas_Cowie_Masters_Thesis and create a virtual
environment there

```
virtualenv maud_venv
```

activate the virtual environment

```
source maud_venv/activate/bin
```

After you've created and activated the virtual environment you can install all of the
dependencies for Maud

```
pip install -e .
```

This will install the directories related to cmdstanpy. A library which interfaces
with Stan, the probabilistic programming language.

To install the stan library you need to run the following command:

```
install_cmdstan
```

Further information can be found on the Stan Website: https://mc-stan.org

### Running a Maud model

The first step in generating data is to run the Maud model.
This model contains all of the priors with uncertainties, measurements for conditions,
and model structure and interactions.

The model used in this thesis is found in: "Models/data/july-2020/yeast_ethanol_fx_no_split.toml"

To run the model in Maud, use the following command:

```
 maud sample ../Models/data/july-2020/yeast_ethanol_fx_no_split_compartment.toml --output_dir="../Models/results/july-2020/yeast_test" --n_samples=250 --n_warmup=250 --timepoint=1000 --n_chains=4
```

This will simulate the model and generate the results in "Models/results/july-2020/"

### Validating model

From the Maud directory, with the virtual environment still active. to simulate
Nitrogen limited conditions using the samples, run the following command, whilst changing 
the csv_filename in "Maud/src/maud/validate_samples.py" to the string of numbers
following the csv file names:

```
python src/maud/validate_samples.py
```

In this instance, the numbers should already be set up
to run the pre-generated samples.

This will generate the .csv files used in the analysis in the "Models/results/july-2020"
folder.

### Simulating trained model

Similar to model validation, runnning the command

```
python src/maud/simulate_samples.py
```

Will generate the csv files that define the metabolite concentrations,
free enzyme amount and enzyme fluxes for the samples parameters.

### Plotting 

To plot the results. after running all of the simulations, you should deactivate the
environment and head to the Models directory. To run the jupyter notebook, it's advised
to create an enviroment using the `pipenv install` command. And adding that to jupyter using
iPykernel:

```
source activate Models
python -m ipykernel install --user --name Models --display-name "Models Plotting"
```

You can then open the jupyter notebook in "Models/src/analysis/july-2020/"
and run each of the commands there.

This should plot and save all of the figures into the results section:
"Models/results/july-2020/yeast_ma"




