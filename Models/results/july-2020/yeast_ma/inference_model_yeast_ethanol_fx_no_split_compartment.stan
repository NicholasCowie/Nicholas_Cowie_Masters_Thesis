functions{
#include big_k_rate_equations.stan
#include partial_sums.stan
#include haldane_relationships.stan
#include allostery.stan
  vector get_fluxes(real[] m, real[] p){
  real empty_array[0];
  real Tr_PGI = p[9]*p[55]
                * (p[7]/p[53])^(-1*-1) 
                - p[9]*p[55]/p[31]
                * (p[53])^-1 
                * (p[54])^1 
                * (m[3]/p[54])^1 ;

real Dr_PGI = (1 + p[7]/p[53])^(-1*-1) 
                + (1 + m[3]/p[54])^1 
                - 1;


real Dr_reg_PGI = 0;
  real Tr_PFK = p[10]*p[60]
                * (m[3]/p[56])^(-1*-1) * (p[1]/p[57])^(-1*-1) 
                - p[10]*p[60]/p[32]
                * (p[56])^-1 * (p[57])^-1 
                * (p[58])^1 * (p[59])^1 
                * (m[4]/p[58])^1 * (p[2]/p[59])^1 ;

real Dr_PFK = (1 + m[3]/p[56])^(-1*-1) * (1 + p[1]/p[57])^(-1*-1) 
                + (1 + m[4]/p[58])^1 * (1 + p[2]/p[59])^1 
                - 1;


real Dr_reg_PFK = 0;
  real Tr_FBA = p[11]*p[67]
                * (m[4]/p[64])^(-1*-1) 
                - p[11]*p[67]/p[33]
                * (p[64])^-1 
                * (p[65])^1 * (p[66])^1 
                * (m[7]/p[65])^1 * (m[5]/p[66])^1 ;

real Dr_FBA = (1 + m[4]/p[64])^(-1*-1) 
                + (1 + m[7]/p[65])^1 * (1 + m[5]/p[66])^1 
                - 1;


real Dr_reg_FBA = (m[5]/p[68]) ;


  real Tr_TPI = p[12]*p[71]
                * (m[7]/p[69])^(-1*-1) 
                - p[12]*p[71]/p[34]
                * (p[69])^-1 
                * (p[70])^1 
                * (m[5]/p[70])^1 ;

real Dr_TPI = (1 + m[7]/p[69])^(-1*-1) 
                + (1 + m[5]/p[70])^1 
                - 1;


real Dr_reg_TPI = (m[5]/p[72]) ;


  real Tr_TDH1 = p[13]*p[78]
                * (m[5]/p[73])^(-1*-1) * (p[3]/p[74])^(-1*-1) * (p[8]/p[75])^(-1*-1) 
                - p[13]*p[78]/p[35]
                * (p[73])^-1 * (p[74])^-1 * (p[75])^-1 
                * (p[76])^1 * (p[77])^1 
                * (p[4]/p[76])^1 * (m[6]/p[77])^1 ;

real Dr_TDH1 = (1 + m[5]/p[73])^(-1*-1) * (1 + p[3]/p[74])^(-1*-1) * (1 + p[8]/p[75])^(-1*-1) 
                + (1 + p[4]/p[76])^1 * (1 + m[6]/p[77])^1 
                - 1;


real Dr_reg_TDH1 = 0;
  real Tr_TDH3 = p[14]*p[84]
                * (m[5]/p[79])^(-1*-1) * (p[3]/p[80])^(-1*-1) * (p[8]/p[81])^(-1*-1) 
                - p[14]*p[84]/p[36]
                * (p[79])^-1 * (p[80])^-1 * (p[81])^-1 
                * (p[82])^1 * (p[83])^1 
                * (p[4]/p[82])^1 * (m[6]/p[83])^1 ;

real Dr_TDH3 = (1 + m[5]/p[79])^(-1*-1) * (1 + p[3]/p[80])^(-1*-1) * (1 + p[8]/p[81])^(-1*-1) 
                + (1 + p[4]/p[82])^1 * (1 + m[6]/p[83])^1 
                - 1;


real Dr_reg_TDH3 = 0;
  real Tr_PGK = p[15]*p[89]
                * (p[2]/p[85])^(-1*-1) * (m[6]/p[86])^(-1*-1) 
                - p[15]*p[89]/p[37]
                * (p[85])^-1 * (p[86])^-1 
                * (p[87])^1 * (p[88])^1 
                * (m[1]/p[87])^1 * (p[1]/p[88])^1 ;

real Dr_PGK = (1 + p[2]/p[85])^(-1*-1) * (1 + m[6]/p[86])^(-1*-1) 
                + (1 + m[1]/p[87])^1 * (1 + p[1]/p[88])^1 
                - 1;


real Dr_reg_PGK = 0;
  real Tr_GPM1 = p[16]*p[94]
                * (m[1]/p[92])^(-1*-1) 
                - p[16]*p[94]/p[38]
                * (p[92])^-1 
                * (p[93])^1 
                * (m[2]/p[93])^1 ;

real Dr_GPM1 = (1 + m[1]/p[92])^(-1*-1) 
                + (1 + m[2]/p[93])^1 
                - 1;


real Dr_reg_GPM1 = 0;
  real Tr_ENO1 = p[17]*p[97]
                * (m[2]/p[95])^(-1*-1) 
                - p[17]*p[97]/p[39]
                * (p[95])^-1 
                * (p[96])^1 
                * (m[8]/p[96])^1 ;

real Dr_ENO1 = (1 + m[2]/p[95])^(-1*-1) 
                + (1 + m[8]/p[96])^1 
                - 1;


real Dr_reg_ENO1 = 0;
  real Tr_ENO2 = p[18]*p[100]
                * (m[2]/p[98])^(-1*-1) 
                - p[18]*p[100]/p[40]
                * (p[98])^-1 
                * (p[99])^1 
                * (m[8]/p[99])^1 ;

real Dr_ENO2 = (1 + m[2]/p[98])^(-1*-1) 
                + (1 + m[8]/p[99])^1 
                - 1;


real Dr_reg_ENO2 = 0;
  real Tr_CDC19 = p[19]*p[105]
                * (m[8]/p[101])^(-1*-1) * (p[2]/p[102])^(-1*-1) 
                - p[19]*p[105]/p[41]
                * (p[101])^-1 * (p[102])^-1 
                * (p[103])^1 * (p[104])^1 
                * (p[1]/p[103])^1 * (m[9]/p[104])^1 ;

real Dr_CDC19 = (1 + m[8]/p[101])^(-1*-1) * (1 + p[2]/p[102])^(-1*-1) 
                + (1 + p[1]/p[103])^1 * (1 + m[9]/p[104])^1 
                - 1;


real Dr_reg_CDC19 = 0;
  real Tr_ZWF1 = p[20]*p[112]
                * (p[5]/p[108])^(-1*-1) * (p[7]/p[109])^(-1*-1) 
                - p[20]*p[112]/p[42]
                * (p[108])^-1 * (p[109])^-1 
                * (p[110])^1 * (p[111])^1 
                * (m[10]/p[110])^1 * (p[6]/p[111])^1 ;

real Dr_ZWF1 = (1 + p[5]/p[108])^(-1*-1) * (1 + p[7]/p[109])^(-1*-1) 
                + (1 + m[10]/p[110])^1 * (1 + p[6]/p[111])^1 
                - 1;


real Dr_reg_ZWF1 = 0;
  real Tr_SOL3 = p[21]*p[115]
                * (m[10]/p[113])^(-1*-1) 
                - p[21]*p[115]/p[43]
                * (p[113])^-1 
                * (p[114])^1 
                * (m[11]/p[114])^1 ;

real Dr_SOL3 = (1 + m[10]/p[113])^(-1*-1) 
                + (1 + m[11]/p[114])^1 
                - 1;


real Dr_reg_SOL3 = 0;
  real Tr_GND1 = p[22]*p[120]
                * (m[11]/p[116])^(-1*-1) * (p[5]/p[117])^(-1*-1) 
                - p[22]*p[120]/p[44]
                * (p[116])^-1 * (p[117])^-1 
                * (p[118])^1 * (p[119])^1 
                * (p[6]/p[118])^1 * (m[12]/p[119])^1 ;

real Dr_GND1 = (1 + m[11]/p[116])^(-1*-1) * (1 + p[5]/p[117])^(-1*-1) 
                + (1 + p[6]/p[118])^1 * (1 + m[12]/p[119])^1 
                - 1;


real Dr_reg_GND1 = 0;
  real Tr_RPE1 = p[23]*p[123]
                * (m[12]/p[121])^(-1*-1) 
                - p[23]*p[123]/p[45]
                * (p[121])^-1 
                * (p[122])^1 
                * (m[15]/p[122])^1 ;

real Dr_RPE1 = (1 + m[12]/p[121])^(-1*-1) 
                + (1 + m[15]/p[122])^1 
                - 1;


real Dr_reg_RPE1 = 0;
  real Tr_RKI1 = p[24]*p[126]
                * (m[12]/p[124])^(-1*-1) 
                - p[24]*p[126]/p[46]
                * (p[124])^-1 
                * (p[125])^1 
                * (m[13]/p[125])^1 ;

real Dr_RKI1 = (1 + m[12]/p[124])^(-1*-1) 
                + (1 + m[13]/p[125])^1 
                - 1;


real Dr_reg_RKI1 = 0;
  real Tr_TKL1 = p[25]*p[131]
                * (m[15]/p[127])^(-1*-1) * (m[13]/p[128])^(-1*-1) 
                - p[25]*p[131]/p[47]
                * (p[127])^-1 * (p[128])^-1 
                * (p[129])^1 * (p[130])^1 
                * (m[16]/p[129])^1 * (m[5]/p[130])^1 ;

real Dr_TKL1 = (1 + m[15]/p[127])^(-1*-1) * (1 + m[13]/p[128])^(-1*-1) 
                + (1 + m[16]/p[129])^1 * (1 + m[5]/p[130])^1 
                - 1;


real Dr_reg_TKL1 = 0;
  real Tr_TKL2 = p[26]*p[136]
                * (m[14]/p[132])^(-1*-1) * (m[15]/p[133])^(-1*-1) 
                - p[26]*p[136]/p[48]
                * (p[132])^-1 * (p[133])^-1 
                * (p[134])^1 * (p[135])^1 
                * (m[3]/p[134])^1 * (m[5]/p[135])^1 ;

real Dr_TKL2 = (1 + m[14]/p[132])^(-1*-1) * (1 + m[15]/p[133])^(-1*-1) 
                + (1 + m[3]/p[134])^1 * (1 + m[5]/p[135])^1 
                - 1;


real Dr_reg_TKL2 = 0;
  real Tr_TAL1 = p[27]*p[141]
                * (m[5]/p[137])^(-1*-1) * (m[16]/p[138])^(-1*-1) 
                - p[27]*p[141]/p[49]
                * (p[137])^-1 * (p[138])^-1 
                * (p[139])^1 * (p[140])^1 
                * (m[3]/p[139])^1 * (m[14]/p[140])^1 ;

real Dr_TAL1 = (1 + m[5]/p[137])^(-1*-1) * (1 + m[16]/p[138])^(-1*-1) 
                + (1 + m[3]/p[139])^1 * (1 + m[14]/p[140])^1 
                - 1;


real Dr_reg_TAL1 = 0;
  real free_enzyme_ratio_PFK = get_free_enzyme_ratio_modular_rate_law(Tr_PFK, Dr_PFK, Dr_reg_PFK);
  real free_enzyme_ratio_PGK = get_free_enzyme_ratio_modular_rate_law(Tr_PGK, Dr_PGK, Dr_reg_PGK);
  real free_enzyme_ratio_CDC19 = get_free_enzyme_ratio_modular_rate_law(Tr_CDC19, Dr_CDC19, Dr_reg_CDC19);
  return [
    modular_rate_law(Tr_PGI, Dr_PGI, Dr_reg_PGI),
    modular_rate_law(Tr_PFK, Dr_PFK, Dr_reg_PFK)*get_regulatory_effect({m[4]},{p[1]},free_enzyme_ratio_PFK,{p[61]},{p[63]},p[62],8),
    modular_rate_law(Tr_FBA, Dr_FBA, Dr_reg_FBA),
    modular_rate_law(Tr_TPI, Dr_TPI, Dr_reg_TPI),
    modular_rate_law(Tr_TDH1, Dr_TDH1, Dr_reg_TDH1)+modular_rate_law(Tr_TDH3, Dr_TDH3, Dr_reg_TDH3),
    modular_rate_law(Tr_PGK, Dr_PGK, Dr_reg_PGK)*get_regulatory_effect({p[2]},empty_array,free_enzyme_ratio_PGK,{p[90]},empty_array,p[91],2),
    modular_rate_law(Tr_GPM1, Dr_GPM1, Dr_reg_GPM1),
    modular_rate_law(Tr_ENO1, Dr_ENO1, Dr_reg_ENO1)+modular_rate_law(Tr_ENO2, Dr_ENO2, Dr_reg_ENO2),
    modular_rate_law(Tr_CDC19, Dr_CDC19, Dr_reg_CDC19)*get_regulatory_effect({m[4]},empty_array,free_enzyme_ratio_CDC19,{p[106]},empty_array,p[107],4),
    modular_rate_law(Tr_ZWF1, Dr_ZWF1, Dr_reg_ZWF1),
    modular_rate_law(Tr_SOL3, Dr_SOL3, Dr_reg_SOL3),
    modular_rate_law(Tr_GND1, Dr_GND1, Dr_reg_GND1),
    modular_rate_law(Tr_RPE1, Dr_RPE1, Dr_reg_RPE1),
    modular_rate_law(Tr_RKI1, Dr_RKI1, Dr_reg_RKI1),
    modular_rate_law(Tr_TKL1, Dr_TKL1, Dr_reg_TKL1),
    modular_rate_law(Tr_TKL2, Dr_TKL2, Dr_reg_TKL2),
    modular_rate_law(Tr_TAL1, Dr_TAL1, Dr_reg_TAL1),
    irr_mass_action(m[9],p[28]),
    irr_mass_action(m[7],p[29]),
    irr_mass_action(m[13],p[30])
  ]';
}
  real[] ode_func(real t, real[] m, real[] p, real[] xr, int[] xi){
  vector[20] fluxes = get_fluxes(m, p);
  return {
    1*fluxes[6]-1*fluxes[7],
    1*fluxes[7]-1*fluxes[8],
    1*fluxes[1]-1*fluxes[2]+1*fluxes[16]+1*fluxes[17],
    1*fluxes[2]-1*fluxes[3],
    1*fluxes[3]+1*fluxes[4]-1*fluxes[5]+1*fluxes[15]+1*fluxes[16]-1*fluxes[17],
    1*fluxes[5]-1*fluxes[6],
    1*fluxes[3]-1*fluxes[4]-1*fluxes[19],
    1*fluxes[8]-1*fluxes[9],
    1*fluxes[9]-1*fluxes[18],
    1*fluxes[10]-1*fluxes[11],
    1*fluxes[11]-1*fluxes[12],
    1*fluxes[12]-1*fluxes[13]-1*fluxes[14],
    1*fluxes[14]-1*fluxes[15]-1*fluxes[20],
    -1*fluxes[16]+1*fluxes[17],
    1*fluxes[13]-1*fluxes[15]-1*fluxes[16],
    1*fluxes[15]-1*fluxes[17]
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
  real yconc[N_conc_measurement];
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  real yflux[N_flux_measurement];
  vector<lower=0>[N_flux_measurement] sigma_flux;
  int<lower=1,upper=N_experiment> experiment_yenz[N_enzyme_measurement];
  int<lower=1,upper=N_enzyme> enzyme_yenz[N_enzyme_measurement];
  real yenz[N_enzyme_measurement];
  vector<lower=0>[N_enzyme_measurement] sigma_enz;
  // hardcoded priors
  vector[N_metabolite] formation_energy;
  vector[N_kinetic_parameters] prior_loc_kinetic_parameters;
  vector<lower=0>[N_kinetic_parameters] prior_scale_kinetic_parameters;
  real prior_loc_unbalanced[N_experiment, N_unbalanced];
  real<lower=0> prior_scale_unbalanced[N_experiment, N_unbalanced];
  real prior_loc_enzyme[N_experiment, N_enzyme];
  real<lower=0> prior_scale_enzyme[N_experiment, N_enzyme];
  // network properties
  matrix[N_mic, N_enzyme] stoichiometric_matrix;
  int<lower=1,upper=N_metabolite> metabolite_ix_stoichiometric_matrix[N_mic];
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] knockout_enzymes;
  vector[N_enzyme] adjustment;
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
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  matrix<lower=0>[N_experiment, N_enzyme] enzyme_concentration;
  matrix<lower=0>[N_experiment, N_unbalanced] conc_unbalanced;
}
transformed parameters {
  real initial_time = 0;
  vector<lower=0>[N_mic] conc[N_experiment];
  matrix[N_experiment, N_reaction] flux;
  vector[N_enzyme] delta_g = stoichiometric_matrix' * formation_energy[metabolite_ix_stoichiometric_matrix] + adjustment;
  vector[N_enzyme] keq = exp(delta_g / minus_RT);

  for (e in 1:N_experiment){
    vector[N_enzyme] temp_enzyme_concentration = (enzyme_concentration[e, ] .* knockout_enzymes[e, ])';
    vector[N_unbalanced+N_enzyme+N_enzyme+N_kinetic_parameters] theta = append_row(append_row(append_row(
      conc_unbalanced[e, ]', temp_enzyme_concentration), keq), kinetic_parameters);
    conc[e, balanced_mic_ix] = to_vector(integrate_ode_bdf(
                                    ode_func,
                                    conc_init[e,],
                                    initial_time,
                                    rep_array(timepoint, 1),
                                    to_array_1d(theta),
                                    xr,
                                    rep_array(0, 1),
                                    1e-10, 1e-12, 1e5
                                  )[1, ]); 
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e, ]';
    flux[e, ] = get_fluxes(to_array_1d(conc[e, balanced_mic_ix]), to_array_1d(theta))';
  }
}
model {
  kinetic_parameters ~ lognormal(log(prior_loc_kinetic_parameters), prior_scale_kinetic_parameters);
  for (e in 1:N_experiment){
    conc_unbalanced[e] ~ lognormal(log(prior_loc_unbalanced[e]), prior_scale_unbalanced[e]);
    enzyme_concentration[e] ~ lognormal(log(prior_loc_enzyme[e]), prior_scale_enzyme[e]);
  }
  if (LIKELIHOOD == 1){
    target += reduce_sum(partial_sum_conc, yconc, 1,
                         conc, experiment_yconc, mic_ix_yconc, sigma_conc);
    target += reduce_sum(partial_sum_enz, yenz, 1,
                         enzyme_concentration, experiment_yenz, enzyme_yenz, sigma_enz);
    target += reduce_sum(partial_sum_flux, yflux, 1,
                         flux, experiment_yflux, reaction_yflux, sigma_flux);
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_enzyme_measurement] yenz_sim;
  vector[N_flux_measurement] yflux_sim;
  vector[N_flux_measurement+N_conc_measurement] log_like;


  for (f in 1:N_flux_measurement){
    log_like[f] = normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }

  for (c in 1:N_conc_measurement){
    log_like[N_flux_measurement+c] = lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
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