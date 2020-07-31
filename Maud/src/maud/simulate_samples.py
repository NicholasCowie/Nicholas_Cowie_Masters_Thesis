import os
import json

import arviz as az
import pandas as pd
import numpy as np
import cmdstanpy

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
toml_input = f'{model_name}_wo_allo.toml'
csv_filename = f"inference_model_{model_name}-202007241852"
inference_json = f"input_data_{model_name}.json"
stan_filename = f"simulate_model_{model_name}.stan"

timepoints = [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000]

# Sensitivity analysis input
""""
where each key represents an unbalanced metabolite,
and the value represents the [amount]*[sampled_value].
For example:
    {'g6p_c': 2} means that the g6p_c concentration
    will be doubled
"""
modifications = {'g6p_c': 1}


## Functions
def get_keys(stan_codes):
    return [key for key in stan_codes.keys()]


def get_input(gibbs_energy,
              kinetic_parameters,
              conc_sim,
              enzyme_concentration,
              timepoints,
              ii):
    
    return {
    "N_mic": ii["N_mic"],
    "N_unbalanced": ii["N_unbalanced"],
    "N_kinetic_parameters": ii["N_kinetic_parameters"],
    "N_reaction": ii["N_reaction"],
    "N_enzyme": ii["N_enzyme"],
    "N_experiment": ii["N_experiment"],
    "N_flux_measurement": ii["N_flux_measurement"],
    "N_enzyme_measurement": ii["N_enzyme_measurement"],
    "N_conc_measurement": ii["N_conc_measurement"],
    "N_metabolite": ii["N_metabolite"],
    "N_timepoint": len(timepoints),
    "unbalanced_mic_ix": ii["unbalanced_mic_ix"],
    "balanced_mic_ix": ii["balanced_mic_ix"],
    "experiment_yconc": ii["experiment_yconc"],
    "mic_ix_yconc": ii["mic_ix_yconc"],
    "reaction_yflux": ii["reaction_yflux"],
    "experiment_yenz": ii["experiment_yenz"],
    "enzyme_yenz": ii["enzyme_yenz"],
    "kinetic_parameters": kinetic_parameters,
    "enzyme_concentration": enzyme_concentration,
    "conc_sim": conc_sim,
    "delta_g": gibbs_energy,
    "timepoints": timepoints,
    "adjustment": ii["adjustment"]
    }
    


## Script
mi = io.load_maud_input_from_toml(os.path.join(PATHS['DATA'], toml_input))

experiment_codes = get_keys(mi.stan_codes["experiment"])
reaction_codes = get_keys(mi.stan_codes["reaction"])
reaction_index = mi.stan_codes["reaction"]
enzyme_codes = get_keys(mi.stan_codes["enzyme"])
enzyme_index = mi.stan_codes["enzyme"]
metabolite_codes = get_keys(mi.stan_codes["metabolite"])
kinetic_parameter_codes = get_keys(mi.stan_codes["kinetic_parameter"])
mic_codes = get_keys(mi.stan_codes["metabolite_in_compartment"])
mic_index = mi.stan_codes["metabolite_in_compartment"]
balanced_mic_codes = utils.codify(get_keys(mi.stan_codes["balanced_mic"]))

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
    stan_file=os.path.join(PATHS['RESULTS'], stan_filename)
    )

## Initialising empty dataframe for ODE results
# ode_results = pd.DataFrame(columns = ['experiment_code', 'metabolite_code'] + [f'{time}' for time in timepoints])
enzyme_flux_results = pd.DataFrame(columns = ['experiment_code', 'enzyme_code', 'flux', 'draw'])
enzyme_denom_results = pd.DataFrame(columns = ['experiment_code', 'enzyme_code', 'denom', 'draw'])
# flux_results = pd.DataFrame(columns = ['experiment_code', 'reaction_code'] + [f'{time}' for time in timepoints])
# dG_results = pd.DataFrame(columns = ['experiment_code', 'enzyme_code', 'dG_0'] + [f'{time}' for time in timepoints])
# reg_component_results = pd.DataFrame(columns = ['experiment_code', 'enzyme_code'] + [f'{time}' for time in timepoints])
# idm = 0
idef = 0
# idr = 0
# idf = 0
# ida = 0

stoichiometric_matrix_T = np.matrix(inference_input['stoichiometric_matrix']).T

for chain in range(0, 4): # data_trained.posterior['chain']:
    for draw in range(0, 500): # data_trained.posterior['draw']:
        tmp_gibbs_energy = data_trained.posterior.delta_g[chain][draw].values
        tmp_kinetic_parameters = data_trained.posterior.kinetic_parameters[chain][draw].values
        tmp_conc_sim = data_trained.posterior.conc[chain][draw].values
        if len(modifications) > 0:
            for mod_met, mod_val in modifications.items():
                tmp_conc_sim[:, mic_index[mod_met]-1] = np.multiply(mod_val, tmp_conc_sim[:, mic_index[mod_met]-1])

        tmp_enzyme_concentration = data_trained.posterior.enzyme_concentration[chain][draw].values
        tmp_input = get_input(gibbs_energy = tmp_gibbs_energy,
                              kinetic_parameters = tmp_kinetic_parameters,
                              conc_sim = tmp_conc_sim,
                              enzyme_concentration = tmp_enzyme_concentration,
                              timepoints=timepoints,
                              ii = inference_input)

        fit = model.sample(
            data=tmp_input,
            chains=1,
            save_warmup=False,
            fixed_param=True
            )
        
        tmp_enzyme_flux_results = {key[11:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'enzyme_flux' in key}
        tmp_enzyme_denom_results = {key[12:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'enzyme_denom' in key}
        # tmp_ode_results = {key[4:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'conc' in key}
        # tmp_flux_results = {key[4:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'flux' in key}
        # tmp_dG_0_results = tmp_gibbs_energy
        # tmp_reg_results = {key[4:]: val for key, val in fit.summary()['Mean'].to_dict().items() if 'allo' in key}
        # for exp_id, exp in utils.codify(experiment_codes).items():
        #     for met_id, met in mic_index.items():
        #         for time_step, time in enumerate(timepoints):
        #             ode_results.loc[idm, 'experiment_code'] = exp_id
        #             ode_results.loc[idm, 'metabolite_code'] = met_id
        #             ode_results.loc[idm, f'{time}'] = tmp_ode_results[f'[{exp},{time_step+1},{met}]']


        #         idm += 1

        # for exp_id, exp in utils.codify(experiment_codes).items():
        #     for rxn_id, rxn in reaction_index.items():
        #         for time_step, time in enumerate(timepoints):
        #             flux_results.loc[idr, 'experiment_code'] = exp_id
        #             flux_results.loc[idr, 'reaction_code'] = rxn_id
        #             flux_results.loc[idr, f'{time}'] = tmp_flux_results[f'[{exp},{time_step+1},{rxn}]']

        #         idr += 1

        for exp_id, exp in utils.codify(experiment_codes).items():
            for enz_id, enz in enzyme_index.items():
                enzyme_flux_results.loc[idef, 'experiment_code'] = exp_id
                enzyme_flux_results.loc[idef, 'enzyme_code'] = enz_id
                enzyme_flux_results.loc[idef, 'draw'] = draw+500*chain
                enzyme_flux_results.loc[idef, 'flux'] = tmp_enzyme_flux_results[f'[{exp},{enz}]']

                idef += 1

        for exp_id, exp in utils.codify(experiment_codes).items():
            for enz_id, enz in enzyme_index.items():
                enzyme_denom_results.loc[idef, 'experiment_code'] = exp_id
                enzyme_denom_results.loc[idef, 'enzyme_code'] = enz_id
                enzyme_denom_results.loc[idef, 'draw'] = draw+500*chain
                enzyme_denom_results.loc[idef, 'denom'] = tmp_enzyme_denom_results[f'[{exp},{enz}]']

                idef += 1

        # for exp_id, exp in utils.codify(experiment_codes).items():
        #     for rxn_id, rxn in reaction_index.items():
        #         for time_step, time in enumerate(timepoints):
        #             reg_component_results.loc[ida, 'experiment_code'] = exp_id
        #             reg_component_results.loc[ida, 'enzyme_code'] = rxn_id
        #             reg_component_results.loc[ida, f'{time}'] = tmp_reg_results[f'[{exp},{time_step+1},{rxn}]']

        #         ida += 1

        # for exp_id, exp in utils.codify(experiment_codes).items():
        #     tmp_conc_vector = [conc for key, conc in tmp_ode_results.items() if f'[{exp}' in key]
        #     tmp_conc_matrix = np.reshape(tmp_conc_vector, (len(mic_codes), len(timepoints)), order='F')
        #     tmp_Q = np.multiply(298.15*0.008314, np.matmul(stoichiometric_matrix_T, np.log(tmp_conc_matrix)))
        #     for rxn_id, rxn in reaction_index.items():
        #         for time_step, time in enumerate(timepoints):
        #             dG_results.loc[idr, 'experiment_code'] = exp_id
        #             dG_results.loc[idr, 'enzyme_code'] = rxn_id
        #             dG_results.loc[idr, 'dG_0'] = tmp_gibbs_energy[rxn-1]
        #             dG_results.loc[idr, f'{time}'] = tmp_gibbs_energy[rxn-1] + tmp_Q[rxn-1, time_step]

        #         idf += 1


# Exporting ODE Results
# ode_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_ode_results.csv"))
# dG_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_dG_results.csv"))
# flux_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_flux_results.csv"))
enzyme_flux_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_enzyme_flux_results.csv"))
enzyme_denom_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_enzyme_denom_results.csv"))
# reg_component_results.to_csv(os.path.join(PATHS['RESULTS'], f"{model_name}_reg_results.csv"))
## Analysing dG




