functions{
#include big_k_rate_equations.stan
#include haldane_relationships.stan
#include allostery.stan
  vector get_fluxes(real[] m, real[] p){
  real empty_array[0];
  real Tr_GLCt1 = p[9]*p[27]
                * (p[1]/p[25])^(-1*-1) 
                - p[9]*p[27]/p[17]
                * (p[25])^-1 
                * (p[26])^1 
                * (m[1]/p[26])^1 ;

real Dr_GLCt1 = (1 + p[1]/p[25])^(-1*-1) 
                + (1 + m[1]/p[26])^1 
                - 1;


real Dr_reg_GLCt1 = 0;
  real Tr_HEX1 = p[10]*p[32]
                * (m[1]/p[28])^(-1*-1) * (p[2]/p[29])^(-1*-1) 
                - p[10]*p[32]/p[18]
                * (p[28])^-1 * (p[29])^-1 
                * (p[30])^1 * (p[31])^1 
                * (m[2]/p[30])^1 * (p[3]/p[31])^1 ;

real Dr_HEX1 = (1 + m[1]/p[28])^(-1*-1) * (1 + p[2]/p[29])^(-1*-1) 
                + (1 + m[2]/p[30])^1 * (1 + p[3]/p[31])^1 
                - 1;


real Dr_reg_HEX1 = 0;
  real Tr_PGI = p[11]*p[35]
                * (m[2]/p[33])^(-1*-1) 
                - p[11]*p[35]/p[19]
                * (p[33])^-1 
                * (p[34])^1 
                * (m[3]/p[34])^1 ;

real Dr_PGI = (1 + m[2]/p[33])^(-1*-1) 
                + (1 + m[3]/p[34])^1 
                - 1;


real Dr_reg_PGI = 0;
  real Tr_PFK = p[12]*p[40]
                * (m[3]/p[36])^(-1*-1) * (p[2]/p[37])^(-1*-1) 
                - p[12]*p[40]/p[20]
                * (p[36])^-1 * (p[37])^-1 
                * (p[38])^1 * (p[39])^1 
                * (m[4]/p[38])^1 * (p[3]/p[39])^1 ;

real Dr_PFK = (1 + m[3]/p[36])^(-1*-1) * (1 + p[2]/p[37])^(-1*-1) 
                + (1 + m[4]/p[38])^1 * (1 + p[3]/p[39])^1 
                - 1;


real Dr_reg_PFK = 0;
  real Tr_FBA = p[13]*p[44]
                * (m[4]/p[41])^(-1*-1) 
                - p[13]*p[44]/p[21]
                * (p[41])^-1 
                * (p[42])^1 * (p[43])^1 
                * (p[8]/p[42])^1 * (m[5]/p[43])^1 ;

real Dr_FBA = (1 + m[4]/p[41])^(-1*-1) 
                + (1 + p[8]/p[42])^1 * (1 + m[5]/p[43])^1 
                - 1;


real Dr_reg_FBA = 0;
  real Tr_TPI = p[14]*p[47]
                * (p[8]/p[45])^(-1*-1) 
                - p[14]*p[47]/p[22]
                * (p[45])^-1 
                * (p[46])^1 
                * (m[5]/p[46])^1 ;

real Dr_TPI = (1 + p[8]/p[45])^(-1*-1) 
                + (1 + m[5]/p[46])^1 
                - 1;


real Dr_reg_TPI = 0;
  real Tr_GAPD = p[15]*p[53]
                * (m[5]/p[48])^(-1*-1) * (p[4]/p[49])^(-1*-1) * (p[7]/p[50])^(-1*-1) 
                - p[15]*p[53]/p[23]
                * (p[48])^-1 * (p[49])^-1 * (p[50])^-1 
                * (p[51])^1 * (p[52])^1 
                * (p[5]/p[51])^1 * (m[6]/p[52])^1 ;

real Dr_GAPD = (1 + m[5]/p[48])^(-1*-1) * (1 + p[4]/p[49])^(-1*-1) * (1 + p[7]/p[50])^(-1*-1) 
                + (1 + p[5]/p[51])^1 * (1 + m[6]/p[52])^1 
                - 1;


real Dr_reg_GAPD = 0;
  real Tr_PGK = p[16]*p[58]
                * (p[3]/p[54])^(-1*-1) * (m[6]/p[55])^(-1*-1) 
                - p[16]*p[58]/p[24]
                * (p[54])^-1 * (p[55])^-1 
                * (p[56])^1 * (p[57])^1 
                * (p[6]/p[56])^1 * (p[2]/p[57])^1 ;

real Dr_PGK = (1 + p[3]/p[54])^(-1*-1) * (1 + m[6]/p[55])^(-1*-1) 
                + (1 + p[6]/p[56])^1 * (1 + p[2]/p[57])^1 
                - 1;


real Dr_reg_PGK = 0;
  return [
    modular_rate_law(Tr_GLCt1, Dr_GLCt1, Dr_reg_GLCt1),
    modular_rate_law(Tr_HEX1, Dr_HEX1, Dr_reg_HEX1),
    modular_rate_law(Tr_PGI, Dr_PGI, Dr_reg_PGI),
    modular_rate_law(Tr_PFK, Dr_PFK, Dr_reg_PFK),
    modular_rate_law(Tr_FBA, Dr_FBA, Dr_reg_FBA),
    modular_rate_law(Tr_TPI, Dr_TPI, Dr_reg_TPI),
    modular_rate_law(Tr_GAPD, Dr_GAPD, Dr_reg_GAPD),
    modular_rate_law(Tr_PGK, Dr_PGK, Dr_reg_PGK)
  ]';
}
  real[] ode_func(real t, real[] m, real[] p, real[] xr, int[] xi){
  vector[8] fluxes = get_fluxes(m, p);
  return {
    1*fluxes[1]-1*fluxes[2],
    1*fluxes[2]-1*fluxes[3],
    1*fluxes[3]-1*fluxes[4],
    1*fluxes[4]-1*fluxes[5],
    1*fluxes[5]+1*fluxes[6]-1*fluxes[7],
    1*fluxes[7]-1*fluxes[8]
  };
}
}
data {
  // dimensions
  int<lower=1> N_mic;         // Total number of metabolites in compartments
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameters;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_enzyme_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=1> N_metabolite;  // NB metabolites in multiple compartments only count once here
  // measurements
  int<lower=1,upper=N_mic> unbalanced_mic_ix[N_unbalanced];
  int<lower=1,upper=N_mic> balanced_mic_ix[N_mic-N_unbalanced];
  int<lower=1,upper=N_experiment> experiment_yconc[N_conc_measurement];
  int<lower=1,upper=N_mic> mic_ix_yconc[N_conc_measurement];
  vector[N_conc_measurement] yconc;
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  vector[N_flux_measurement] yflux;
  vector<lower=0>[N_flux_measurement] sigma_flux;
  int<lower=1,upper=N_experiment> experiment_yenz[N_enzyme_measurement];
  int<lower=1,upper=N_enzyme> enzyme_yenz[N_enzyme_measurement];
  vector[N_enzyme_measurement] yenz;
  vector<lower=0>[N_enzyme_measurement] sigma_enz;
  // hardcoded priors
  vector[N_metabolite] prior_loc_formation_energy;
  vector<lower=0>[N_metabolite] prior_scale_formation_energy;
  vector[N_kinetic_parameters] prior_loc_kinetic_parameters;
  vector<lower=0>[N_kinetic_parameters] prior_scale_kinetic_parameters;
  real prior_loc_unbalanced[N_experiment, N_unbalanced];
  real<lower=0> prior_scale_unbalanced[N_experiment, N_unbalanced];
  real prior_loc_enzyme[N_experiment, N_enzyme];
  real<lower=0> prior_scale_enzyme[N_experiment, N_enzyme];
  // network properties
  matrix[N_mic, N_enzyme] stoichiometric_matrix;
  int<lower=1,upper=N_metabolite> metabolite_ix_stoichiometric_matrix[N_mic];
  // configuration
  real<lower=0> conc_init[N_experiment, N_mic-N_unbalanced];
  real rtol;
  real ftol;
  int steps;
  int<lower=0,upper=1> LIKELIHOOD;  // set to 0 for priors-only mode
  real<lower=0> timepoint;
}
transformed data {
  real xr[0];
  int xi[0];
  real minus_RT = - 0.008314 * 298.15;

}
parameters {
  vector[N_metabolite] formation_energy;
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  vector<lower=0>[N_enzyme] enzyme_concentration[N_experiment];
  vector<lower=0>[N_unbalanced] conc_unbalanced[N_experiment];
}
transformed parameters {
  real initial_time = 0;
  vector<lower=0>[N_mic] conc[N_experiment];
  vector[N_reaction] flux[N_experiment];
  vector[N_enzyme] delta_g = stoichiometric_matrix' * formation_energy[metabolite_ix_stoichiometric_matrix];
  for (e in 1:N_experiment){
    vector[N_enzyme] keq = exp(delta_g / minus_RT);
    vector[N_unbalanced+N_enzyme+N_enzyme+N_kinetic_parameters] theta = append_row(append_row(append_row(
      conc_unbalanced[e], enzyme_concentration[e]), keq), kinetic_parameters);
    conc[e, balanced_mic_ix] = to_vector(integrate_ode_bdf(
                                    ode_func,
                                    conc_init[e,],
                                    initial_time,
                                    rep_array(timepoint, 1),
                                    to_array_1d(theta),
                                    xr,
                                    rep_array(0, 1),
                                    1e-8, 1e-12, 1e5
                                  )[1, ]); 
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e];
    flux[e] = get_fluxes(to_array_1d(conc[e, balanced_mic_ix]), to_array_1d(theta));
  }
}
model {
  kinetic_parameters ~ lognormal(log(prior_loc_kinetic_parameters), prior_scale_kinetic_parameters);
  formation_energy ~ normal(prior_loc_formation_energy, prior_scale_formation_energy);
  for (e in 1:N_experiment){
    conc_unbalanced[e] ~ lognormal(log(prior_loc_unbalanced[e]), prior_scale_unbalanced[e]);
    enzyme_concentration[e] ~ lognormal(log(prior_loc_enzyme[e]), prior_scale_enzyme[e]);
  }
  if (LIKELIHOOD == 1){
    for (c in 1:N_conc_measurement){
      target += lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
    }
    for (ec in 1:N_enzyme_measurement){
      target += lognormal_lpdf(yenz[ec] | log(enzyme_concentration[experiment_yenz[ec], enzyme_yenz[ec]]), sigma_enz[ec]);
    }
    for (f in 1:N_flux_measurement){
      target += normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
    }
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_enzyme_measurement] yenz_sim;
  vector[N_flux_measurement] yflux_sim;
  vector[N_flux_measurement+N_conc_measurement] log_like;

  for (c in 1:N_conc_measurement){
    log_like[N_flux_measurement+c] = lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }

  for (f in 1:N_flux_measurement){
    log_like[f] = normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }


  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (ec in 1:N_enzyme_measurement){
    yenz_sim[ec] = lognormal_rng(log(enzyme_concentration[experiment_yenz[ec], enzyme_yenz[ec]]), sigma_enz[ec]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}