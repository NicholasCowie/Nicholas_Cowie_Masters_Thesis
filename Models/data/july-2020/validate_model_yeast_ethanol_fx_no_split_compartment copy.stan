functions{
real modular_rate_law(real Tr, real Dr, real Dr_reg){
  return(Tr/(Dr + Dr_reg));
}

real get_regulatory_effect(real[] activator_concentration,    // metabolite
                           real[] inhibitor_concentration,    // metabolite
                           real free_enzyme_ratio,            // derived from rate equation
                           real[] dissociation_constant_r,    // parameter
                           real[] dissociation_constant_t,    // parameter
                           real transfer_constant){           // parameter
  real Q_num = size(inhibitor_concentration) == 0 ? 1 :
    1 + sum(to_vector(inhibitor_concentration) ./ to_vector(dissociation_constant_t));
  real Q_denom = size(activator_concentration) == 0 ? 1
    : 1 + sum(to_vector(activator_concentration) ./ to_vector(dissociation_constant_r));
  real Q = transfer_constant * free_enzyme_ratio * Q_num / Q_denom;
  return inv(1 + Q); 
}

real get_free_enzyme_ratio_modular_rate_law(real Tr, real Dr, real Dr_reg){
  return 1 / (Dr + Dr_reg);
}

vector get_fluxes(real[] m, real[] p){
real empty_array[0];
real Tr_PGI = p[10]*p[60]
              * (p[7]/p[58])^(-1*-1) 
              - p[10]*p[60]/p[34]
              * (p[58])^-1 
              * (p[59])^1 
              * (m[3]/p[59])^1 ;

real Dr_PGI = (1 + p[7]/p[58])^(-1*-1) 
                + (1 + m[3]/p[59])^1 
                - 1;


real Dr_reg_PGI = 0;
  real Tr_PFK = p[11]*p[65]
                * (m[3]/p[61])^(-1*-1) * (p[1]/p[62])^(-1*-1) 
                - p[11]*p[65]/p[35]
                * (p[61])^-1 * (p[62])^-1 
                * (p[63])^1 * (p[64])^1 
                * (m[4]/p[63])^1 * (p[2]/p[64])^1 ;

real Dr_PFK = (1 + m[3]/p[61])^(-1*-1) * (1 + p[1]/p[62])^(-1*-1) 
                + (1 + m[4]/p[63])^1 * (1 + p[2]/p[64])^1 
                - 1;


real Dr_reg_PFK = 0;
  real Tr_FBA = p[12]*p[72]
                * (m[4]/p[69])^(-1*-1) 
                - p[12]*p[72]/p[36]
                * (p[69])^-1 
                * (p[70])^1 * (p[71])^1 
                * (m[7]/p[70])^1 * (m[5]/p[71])^1 ;

real Dr_FBA = (1 + m[4]/p[69])^(-1*-1) 
                + (1 + m[7]/p[70])^1 * (1 + m[5]/p[71])^1 
                - 1;


real Dr_reg_FBA = (m[5]/p[73]) ;


  real Tr_TPI = p[13]*p[76]
                * (m[7]/p[74])^(-1*-1) 
                - p[13]*p[76]/p[37]
                * (p[74])^-1 
                * (p[75])^1 
                * (m[5]/p[75])^1 ;

real Dr_TPI = (1 + m[7]/p[74])^(-1*-1) 
                + (1 + m[5]/p[75])^1 
                - 1;


real Dr_reg_TPI = (m[5]/p[77]) ;


  real Tr_TDH1 = p[14]*p[83]
                * (m[5]/p[78])^(-1*-1) * (p[3]/p[79])^(-1*-1) * (p[8]/p[80])^(-1*-1) 
                - p[14]*p[83]/p[38]
                * (p[78])^-1 * (p[79])^-1 * (p[80])^-1 
                * (p[81])^1 * (p[82])^1 
                * (p[4]/p[81])^1 * (m[6]/p[82])^1 ;

real Dr_TDH1 = (1 + m[5]/p[78])^(-1*-1) * (1 + p[3]/p[79])^(-1*-1) * (1 + p[8]/p[80])^(-1*-1) 
                + (1 + p[4]/p[81])^1 * (1 + m[6]/p[82])^1 
                - 1;


real Dr_reg_TDH1 = 0;
  real Tr_TDH3 = p[15]*p[89]
                * (m[5]/p[84])^(-1*-1) * (p[3]/p[85])^(-1*-1) * (p[8]/p[86])^(-1*-1) 
                - p[15]*p[89]/p[39]
                * (p[84])^-1 * (p[85])^-1 * (p[86])^-1 
                * (p[87])^1 * (p[88])^1 
                * (p[4]/p[87])^1 * (m[6]/p[88])^1 ;

real Dr_TDH3 = (1 + m[5]/p[84])^(-1*-1) * (1 + p[3]/p[85])^(-1*-1) * (1 + p[8]/p[86])^(-1*-1) 
                + (1 + p[4]/p[87])^1 * (1 + m[6]/p[88])^1 
                - 1;


real Dr_reg_TDH3 = 0;
  real Tr_PGK = p[16]*p[94]
                * (p[2]/p[90])^(-1*-1) * (m[6]/p[91])^(-1*-1) 
                - p[16]*p[94]/p[40]
                * (p[90])^-1 * (p[91])^-1 
                * (p[92])^1 * (p[93])^1 
                * (m[1]/p[92])^1 * (p[1]/p[93])^1 ;

real Dr_PGK = (1 + p[2]/p[90])^(-1*-1) * (1 + m[6]/p[91])^(-1*-1) 
                + (1 + m[1]/p[92])^1 * (1 + p[1]/p[93])^1 
                - 1;


real Dr_reg_PGK = 0;
  real Tr_GPM1 = p[17]*p[99]
                * (m[1]/p[97])^(-1*-1) 
                - p[17]*p[99]/p[41]
                * (p[97])^-1 
                * (p[98])^1 
                * (m[2]/p[98])^1 ;

real Dr_GPM1 = (1 + m[1]/p[97])^(-1*-1) 
                + (1 + m[2]/p[98])^1 
                - 1;


real Dr_reg_GPM1 = 0;
  real Tr_ENO1 = p[18]*p[102]
                * (m[2]/p[100])^(-1*-1) 
                - p[18]*p[102]/p[42]
                * (p[100])^-1 
                * (p[101])^1 
                * (m[8]/p[101])^1 ;

real Dr_ENO1 = (1 + m[2]/p[100])^(-1*-1) 
                + (1 + m[8]/p[101])^1 
                - 1;


real Dr_reg_ENO1 = 0;
  real Tr_ENO2 = p[19]*p[105]
                * (m[2]/p[103])^(-1*-1) 
                - p[19]*p[105]/p[43]
                * (p[103])^-1 
                * (p[104])^1 
                * (m[8]/p[104])^1 ;

real Dr_ENO2 = (1 + m[2]/p[103])^(-1*-1) 
                + (1 + m[8]/p[104])^1 
                - 1;


real Dr_reg_ENO2 = 0;
  real Tr_CDC19 = p[20]*p[110]
                * (m[8]/p[106])^(-1*-1) * (p[2]/p[107])^(-1*-1) 
                - p[20]*p[110]/p[44]
                * (p[106])^-1 * (p[107])^-1 
                * (p[108])^1 * (p[109])^1 
                * (p[1]/p[108])^1 * (m[9]/p[109])^1 ;

real Dr_CDC19 = (1 + m[8]/p[106])^(-1*-1) * (1 + p[2]/p[107])^(-1*-1) 
                + (1 + p[1]/p[108])^1 * (1 + m[9]/p[109])^1 
                - 1;


real Dr_reg_CDC19 = 0;
  real Tr_ZWF1 = p[21]*p[117]
                * (p[5]/p[113])^(-1*-1) * (p[7]/p[114])^(-1*-1) 
                - p[21]*p[117]/p[45]
                * (p[113])^-1 * (p[114])^-1 
                * (p[115])^1 * (p[116])^1 
                * (m[10]/p[115])^1 * (p[6]/p[116])^1 ;

real Dr_ZWF1 = (1 + p[5]/p[113])^(-1*-1) * (1 + p[7]/p[114])^(-1*-1) 
                + (1 + m[10]/p[115])^1 * (1 + p[6]/p[116])^1 
                - 1;


real Dr_reg_ZWF1 = 0;
  real Tr_SOL3 = p[22]*p[120]
                * (m[10]/p[118])^(-1*-1) 
                - p[22]*p[120]/p[46]
                * (p[118])^-1 
                * (p[119])^1 
                * (m[11]/p[119])^1 ;

real Dr_SOL3 = (1 + m[10]/p[118])^(-1*-1) 
                + (1 + m[11]/p[119])^1 
                - 1;


real Dr_reg_SOL3 = 0;
  real Tr_GND1 = p[23]*p[125]
                * (m[11]/p[121])^(-1*-1) * (p[5]/p[122])^(-1*-1) 
                - p[23]*p[125]/p[47]
                * (p[121])^-1 * (p[122])^-1 
                * (p[123])^1 * (p[124])^1 
                * (p[6]/p[123])^1 * (m[12]/p[124])^1 ;

real Dr_GND1 = (1 + m[11]/p[121])^(-1*-1) * (1 + p[5]/p[122])^(-1*-1) 
                + (1 + p[6]/p[123])^1 * (1 + m[12]/p[124])^1 
                - 1;


real Dr_reg_GND1 = 0;
  real Tr_RPE1 = p[24]*p[128]
                * (m[12]/p[126])^(-1*-1) 
                - p[24]*p[128]/p[48]
                * (p[126])^-1 
                * (p[127])^1 
                * (m[15]/p[127])^1 ;

real Dr_RPE1 = (1 + m[12]/p[126])^(-1*-1) 
                + (1 + m[15]/p[127])^1 
                - 1;


real Dr_reg_RPE1 = 0;
  real Tr_RKI1 = p[25]*p[131]
                * (m[12]/p[129])^(-1*-1) 
                - p[25]*p[131]/p[49]
                * (p[129])^-1 
                * (p[130])^1 
                * (m[13]/p[130])^1 ;

real Dr_RKI1 = (1 + m[12]/p[129])^(-1*-1) 
                + (1 + m[13]/p[130])^1 
                - 1;


real Dr_reg_RKI1 = 0;
  real Tr_TKL1 = p[26]*p[136]
                * (m[15]/p[132])^(-1*-1) * (m[13]/p[133])^(-1*-1) 
                - p[26]*p[136]/p[50]
                * (p[132])^-1 * (p[133])^-1 
                * (p[134])^1 * (p[135])^1 
                * (m[16]/p[134])^1 * (m[5]/p[135])^1 ;

real Dr_TKL1 = (1 + m[15]/p[132])^(-1*-1) * (1 + m[13]/p[133])^(-1*-1) 
                + (1 + m[16]/p[134])^1 * (1 + m[5]/p[135])^1 
                - 1;


real Dr_reg_TKL1 = 0;
  real Tr_TKL2 = p[27]*p[141]
                * (m[14]/p[137])^(-1*-1) * (m[15]/p[138])^(-1*-1) 
                - p[27]*p[141]/p[51]
                * (p[137])^-1 * (p[138])^-1 
                * (p[139])^1 * (p[140])^1 
                * (m[3]/p[139])^1 * (m[5]/p[140])^1 ;

real Dr_TKL2 = (1 + m[14]/p[137])^(-1*-1) * (1 + m[15]/p[138])^(-1*-1) 
                + (1 + m[3]/p[139])^1 * (1 + m[5]/p[140])^1 
                - 1;


real Dr_reg_TKL2 = 0;
  real Tr_TAL1 = p[28]*p[146]
                * (m[5]/p[142])^(-1*-1) * (m[16]/p[143])^(-1*-1) 
                - p[28]*p[146]/p[52]
                * (p[142])^-1 * (p[143])^-1 
                * (p[144])^1 * (p[145])^1 
                * (m[3]/p[144])^1 * (m[14]/p[145])^1 ;

real Dr_TAL1 = (1 + m[5]/p[142])^(-1*-1) * (1 + m[16]/p[143])^(-1*-1) 
                + (1 + m[3]/p[144])^1 * (1 + m[14]/p[145])^1 
                - 1;


real Dr_reg_TAL1 = 0;
  real Tr_PDC1 = p[29]*p[149]
                * (m[9]/p[147])^(-1*-1) 
                - p[29]*p[149]/p[53]
                * (p[147])^-1 
                * (p[148])^1 
                * (m[17]/p[148])^1 ;

real Dr_PDC1 = (1 + m[9]/p[147])^(-1*-1) 
                + (1 + m[17]/p[148])^1 
                - 1;


real Dr_reg_PDC1 = 0;
  real Tr_ADH1 = p[30]*p[154]
                * (m[17]/p[150])^(-1*-1) * (p[4]/p[151])^(-1*-1) 
                - p[30]*p[154]/p[54]
                * (p[150])^-1 * (p[151])^-1 
                * (p[152])^1 * (p[153])^1 
                * (p[3]/p[152])^1 * (p[9]/p[153])^1 ;

real Dr_ADH1 = (1 + m[17]/p[150])^(-1*-1) * (1 + p[4]/p[151])^(-1*-1) 
                + (1 + p[3]/p[152])^1 * (1 + p[9]/p[153])^1 
                - 1;


real Dr_reg_ADH1 = 0;
real free_enzyme_ratio_PFK = get_free_enzyme_ratio_modular_rate_law(Tr_PFK, Dr_PFK, Dr_reg_PFK);
real free_enzyme_ratio_PGK = get_free_enzyme_ratio_modular_rate_law(Tr_PGK, Dr_PGK, Dr_reg_PGK);
real free_enzyme_ratio_CDC19 = get_free_enzyme_ratio_modular_rate_law(Tr_CDC19, Dr_CDC19, Dr_reg_CDC19);
return [
  modular_rate_law(Tr_PGI, Dr_PGI, Dr_reg_PGI),
  modular_rate_law(Tr_PFK, Dr_PFK, Dr_reg_PFK)*get_regulatory_effect({m[3]},{p[1]},free_enzyme_ratio_PFK,{p[66]},{p[68]},p[67]),
  modular_rate_law(Tr_FBA, Dr_FBA, Dr_reg_FBA),
  modular_rate_law(Tr_TPI, Dr_TPI, Dr_reg_TPI),
  modular_rate_law(Tr_TDH1, Dr_TDH1, Dr_reg_TDH1)+modular_rate_law(Tr_TDH3, Dr_TDH3, Dr_reg_TDH3),
  modular_rate_law(Tr_PGK, Dr_PGK, Dr_reg_PGK)*get_regulatory_effect({p[2]},empty_array,free_enzyme_ratio_PGK,{p[95]},empty_array,p[96]),
  modular_rate_law(Tr_GPM1, Dr_GPM1, Dr_reg_GPM1),
  modular_rate_law(Tr_ENO1, Dr_ENO1, Dr_reg_ENO1)+modular_rate_law(Tr_ENO2, Dr_ENO2, Dr_reg_ENO2),
  modular_rate_law(Tr_CDC19, Dr_CDC19, Dr_reg_CDC19)*get_regulatory_effect({m[4]},empty_array,free_enzyme_ratio_CDC19,{p[111]},empty_array,p[112]),
  modular_rate_law(Tr_ZWF1, Dr_ZWF1, Dr_reg_ZWF1),
  modular_rate_law(Tr_SOL3, Dr_SOL3, Dr_reg_SOL3),
  modular_rate_law(Tr_GND1, Dr_GND1, Dr_reg_GND1),
  modular_rate_law(Tr_RPE1, Dr_RPE1, Dr_reg_RPE1),
  modular_rate_law(Tr_RKI1, Dr_RKI1, Dr_reg_RKI1),
  modular_rate_law(Tr_TKL1, Dr_TKL1, Dr_reg_TKL1),
  modular_rate_law(Tr_TKL2, Dr_TKL2, Dr_reg_TKL2),
  modular_rate_law(Tr_TAL1, Dr_TAL1, Dr_reg_TAL1),
  modular_rate_law(Tr_PDC1, Dr_PDC1, Dr_reg_PDC1),
  modular_rate_law(Tr_ADH1, Dr_ADH1, Dr_reg_ADH1),
]';
}
  real[] ode_func(real t, real[] m, real[] p, real[] xr, int[] xi){
  vector[22] fluxes = get_fluxes(m, p);
  return {
    1*fluxes[6]-1*fluxes[7],
    1*fluxes[7]-1*fluxes[8],
    1*fluxes[1]-1*fluxes[2]+1*fluxes[16]+1*fluxes[17],
    1*fluxes[2]-1*fluxes[3],
    1*fluxes[3]+1*fluxes[4]-1*fluxes[5]+1*fluxes[15]+1*fluxes[16]-1*fluxes[17],
    1*fluxes[5]-1*fluxes[6],
    1*fluxes[3]-1*fluxes[4],
    1*fluxes[8]-1*fluxes[9],
    1*fluxes[9]-1*fluxes[18],
    1*fluxes[10]-1*fluxes[11],
    1*fluxes[11]-1*fluxes[12],
    1*fluxes[12]-1*fluxes[13]-1*fluxes[14],
    1*fluxes[14]-1*fluxes[15],
    -1*fluxes[16]+1*fluxes[17],
    1*fluxes[13]-1*fluxes[15]-1*fluxes[16],
    1*fluxes[15]-1*fluxes[17],
    1*fluxes[18]-1*fluxes[19]
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
  int<lower=1> N_timepoint;
  // measurement indexes
  int<lower=1,upper=N_mic> unbalanced_mic_ix[N_unbalanced];
  int<lower=1,upper=N_mic> balanced_mic_ix[N_mic-N_unbalanced];
  // input parameters and initial values
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  vector<lower=0>[N_enzyme] enzyme_concentration[N_experiment];
  real<lower=0> conc_sim[N_experiment, N_mic];
  vector[N_enzyme] delta_g;
  // configuration
  real<lower=0> timepoints[N_timepoint];
}
transformed data {
  real xr[0];
  int xi[0];
  real minus_RT = - 0.008314 * 298.15;
}
parameters {
  real alpha;
}
transformed parameters {
  real<lower=0> conc[N_experiment, N_timepoint, N_mic];
  real flux[N_experiment, N_timepoint, N_reaction];
  real initial_time = 0;

  for (e in 1:N_experiment){
    vector[N_enzyme] keq = exp(delta_g / minus_RT);
    vector[N_unbalanced+N_enzyme+N_enzyme+N_kinetic_parameters] parameter = append_row(append_row(append_row(
      to_vector(conc_sim[e, unbalanced_mic_ix]), enzyme_concentration[e]), keq), kinetic_parameters);
    conc[e, ,unbalanced_mic_ix] = rep_array(conc_sim[e, unbalanced_mic_ix], N_timepoint);
    conc[e, ,balanced_mic_ix] = integrate_ode_bdf(
                                    ode_func,
                                    conc_sim[e, balanced_mic_ix],
                                    initial_time,
                                    timepoints,
                                    to_array_1d(parameter),
                                    xr,
                                    rep_array(0, 1),
                                    1e-9, 1e-12, 1e5
                                  );
    for (t in 1:N_timepoint){
      flux[e, t, ] = to_array_1d(get_fluxes(conc[e ,t ,balanced_mic_ix], to_array_1d(parameter)));
    }
  }
}
model {
  target += normal_lpdf(0 | alpha, 1);
}
generated quantities {

}