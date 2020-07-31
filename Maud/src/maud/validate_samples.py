import os
import json

import arviz as az
import pandas as pd
import numpy as np
import cmdstanpy
import toml

from maud import io, utils, code_generation

here = os.path.dirname(os.path.realpath(__file__))
home = os.path.join(here, "../../")

PATHS = {
    'DATA': os.path.join(home, "../Models/data/july-2020/"),
    'RESULTS': os.path.join(home, "../Models/results/july-2020/yeast_ma"),
    'DOC': os.path.join(home, "../Models/doc")
}

## Input Section and Configuration
model_name = "yeast_ethanol_fx_no_split_compartment"
validation_name = "yeast_validation"
toml_input = f'{model_name}_wo_allo.toml'
validation_input = f'{validation_name}.toml'
# validation_file = os.path.join(PATHS['RESULTS'], f"yeast_validation_model.toml")
csv_filename = f"inference_model_{model_name}-202007241851"
inference_json = f"input_data_{model_name}.json"
stan_filename = f"validate_model_{model_name}.stan"

timepoints = [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000]

# Sensitivity analysis input
""""
where each key represents an unbalanced metabolite,
and the value represents the [amount]*[sampled_value].
For example:
    {'g6p_c': 2} means that the g6p_c concentration
    will be doubled
"""
# def get_validation_input(input_file):

#     input_dict = toml.load(input_file)

#     enzyme_dict = input_dict["enzyme"]
#     unbalanced_metabolite_dict = input_dict["unbalanced_metabolites"]

#     enzyme_df = pd.DataFrame.from_dict(enzyme_dict)
#     unbalanced_metabolite_df = pd.DataFrame.from_dict(unbalanced_metabolite_dict)

#     return enzyme_df, unbalanced_metabolite_df



## Functions
def get_keys(stan_codes):
    return [key for key in stan_codes.keys()]


def get_input(gibbs_energy,
              kinetic_parameters,
              conc_sim,
              enzyme_concentration,
              timepoints,
              ii,
              mi_validation):
    
    return {
    "N_mic": ii["N_mic"],
    "N_unbalanced": ii["N_unbalanced"],
    "N_kinetic_parameters": ii["N_kinetic_parameters"],
    "N_reaction": ii["N_reaction"],
    "N_enzyme": ii["N_enzyme"],
    "N_experiment": len(mi_validation.experiments),
    "N_flux_measurement": 1,
    "N_enzyme_measurement": len(enzyme_concentration),
    "N_conc_measurement": len(conc_sim),
    "N_metabolite": ii["N_metabolite"],
    "N_timepoint": len(timepoints),
    "unbalanced_mic_ix": ii["unbalanced_mic_ix"],
    "balanced_mic_ix": ii["balanced_mic_ix"],
    "kinetic_parameters": kinetic_parameters,
    "enzyme_concentration": enzyme_concentration,
    "conc_sim": conc_sim,
    "delta_g": gibbs_energy,
    "timepoints": timepoints,
    "adjustment": ii["adjustment"]
    }
    


## Script
mi = io.load_maud_input_from_toml(os.path.join(PATHS['DATA'], toml_input))
mi_validation = io.load_maud_input_from_toml(os.path.join(PATHS['DATA'], validation_input))

# Training Dataset
experiment_codes = get_keys(mi.stan_codes["experiment"])

# Validation Dataset
experiment_codes_validation = get_keys(mi_validation.stan_codes["experiment"])


# General Model Parameters
reaction_codes = get_keys(mi.stan_codes["reaction"])
reaction_index = mi.stan_codes["reaction"]
enzyme_codes = get_keys(mi.stan_codes["enzyme"])
enzyme_index = mi.stan_codes["enzyme"]
metabolite_codes = get_keys(mi.stan_codes["metabolite"])
kinetic_parameter_codes = get_keys(mi.stan_codes["kinetic_parameter"])
mic_codes = get_keys(mi.stan_codes["metabolite_in_compartment"])
mic_index = mi.stan_codes["metabolite_in_compartment"]
balanced_mic_codes = utils.codify(get_keys(mi.stan_codes["balanced_mic"]))
unbalanced_mic_codes = utils.codify(get_keys(mi.stan_codes["unbalanced_mic"]))

enzyme_list = pd.DataFrame.from_dict(enzyme_index, orient='index').reset_index()
metabolite_list = pd.DataFrame.from_dict(mic_index, orient='index').reset_index()


mic_measurements, reaction_measurements, enzyme_measurements = (
        pd.DataFrame(
            [
                [exp.id, meas.target_id, meas.value, meas.uncertainty]
                for exp in mi_validation.experiments.values()
                for meas in exp.measurements[measurement_type].values()
            ],
            columns=["experiment_id", "target_id", "value", "uncertainty"],
        )
        for measurement_type in ["metabolite", "reaction", "enzyme"]
    )

conc_init = pd.DataFrame(
        0.01, index=mi_validation.stan_codes["experiment"].values(), columns=mi.stan_codes["metabolite_in_compartment"].values()
    )
for _i, row in mic_measurements.iterrows():
    if row["target_id"] in utils.codify(mic_codes).keys():
        row_ix = mi_validation.stan_codes["experiment"][row["experiment_id"]]
        column_ix = mi.stan_codes["metabolite_in_compartment"][row["target_id"]]
        conc_init.loc[row_ix, column_ix] = row["value"]


files_trained = [os.path.join(PATHS['RESULTS'], f"{csv_filename}-{i}.csv") for i in range(1,5)]
data_trained = az.from_cmdstan(
        posterior=files_trained,
        coords={
            'reactions': reaction_codes,
            'enzymes': enzyme_codes,
            'mic': mic_codes,
            'metabolites': metabolite_codes,
            'experiments': experiment_codes,
            'parameter_names': kinetic_parameter_codes,
        },
        dims={
            'conc': ['experiments', 'mic'],
            'flux': ['experiments', 'reactions'],
            'kinetic_parameters': ['parameter_names'],
            'delta_g': ['enzymes'],
            'formation_energy': ['metabolites'],
        }
    )

with open(os.path.join(PATHS["RESULTS"], inference_json), 'r') as json_file:
    inference_input = json.load(json_file)

model = cmdstanpy.CmdStanModel(
    stan_file=os.path.join(PATHS['DATA'], stan_filename)
    )

## Initialising empty dataframe for ODE results
ode_results = pd.DataFrame(columns = ['experiment_code', 'metabolite_code'] + [f'{time}' for time in timepoints])
flux_results = pd.DataFrame(columns = ['experiment_code', 'reaction_code'] + [f'{time}' for time in timepoints])

idm = 0
idr = 0

enzyme_data =  pd.DataFrame(
        0.01, index=mi_validation.stan_codes["experiment"].values(), columns=mi.stan_codes["enzyme"].values()
    )

for _i, row in enzyme_measurements.iterrows():
    if row["target_id"] in utils.codify(enzyme_codes).keys():
        row_ix = mi_validation.stan_codes["experiment"][row["experiment_id"]]
        column_ix = mi.stan_codes["enzyme"][row["target_id"]]
        enzyme_data.loc[row_ix, column_ix] = row["value"]


unbalanced_metabolite_data = pd.DataFrame(
        0.05, index=mi_validation.stan_codes["experiment"].values(), columns=mi.stan_codes["metabolite_in_compartment"].values()
    ) # Median of nadph

for _i, row in mic_measurements.iterrows():
    if row["target_id"] in utils.codify(mic_codes).keys():
        row_ix = mi_validation.stan_codes["experiment"][row["experiment_id"]]
        column_ix = mi.stan_codes["metabolite_in_compartment"][row["target_id"]]
        unbalanced_metabolite_data.loc[row_ix, column_ix] = row["value"]


for chain in range(0, 4): # data_trained.posterior['chain']:
    for draw in range(0, 50): # data_trained.posterior['draw']:
        tmp_gibbs_energy = data_trained.posterior.delta_g[chain][draw].values
        tmp_kinetic_parameters = data_trained.posterior.kinetic_parameters[chain][draw].values

        tmp_input = get_input(gibbs_energy = tmp_gibbs_energy,
                              kinetic_parameters = tmp_kinetic_parameters, 
                              conc_sim = unbalanced_metabolite_data.values,
                              enzyme_concentration = enzyme_data.values,
                              timepoints=timepoints,
                              ii = inference_input,
                              mi_validation = mi_validation)

        fit = model.sample(
            data=tmp_input,
            chains=1,
            save_warmup=False,
            fixed_param=True
            )

        tmp_ode_results = {key[4:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'conc' in key}
        tmp_flux_results = {key[4:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'flux' in key}

        for exp_id, exp in utils.codify(experiment_codes_validation).items():
            for met_id, met in mic_index.items():
                for time_step, time in enumerate(timepoints):
                    ode_results.loc[idm, 'experiment_code'] = exp_id
                    ode_results.loc[idm, 'metabolite_code'] = met_id
                    ode_results.loc[idm, f'{time}'] = tmp_ode_results[f'[{exp},{time_step+1},{met}]']

                idm += 1


        for exp_id, exp in utils.codify(experiment_codes_validation).items():
            for rxn_id, rxn in reaction_index.items():
                for time_step, time in enumerate(timepoints):
                    flux_results.loc[idr, 'experiment_code'] = exp_id
                    flux_results.loc[idr, 'reaction_code'] = rxn_id
                    flux_results.loc[idr, f'{time}'] = tmp_flux_results[f'[{exp},{time_step+1},{rxn}]']

                idr += 1


# Exporting ODE Results
ode_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_ode_results_validation.csv"))
flux_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_flux_results_validation.csv"))
## Analysing dG




