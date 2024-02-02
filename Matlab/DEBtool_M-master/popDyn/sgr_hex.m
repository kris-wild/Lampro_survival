%% sgr_hex
% Gets specific population growth rate for the hex model

%%
function [r, info] = sgr_hex (par, T_pop, f_pop)
  % created 2019/08/01 by Bas Kooijman
  
  %% Syntax
  % [r, info] = <../sgr_hex.m *sgr_hex*> (par, T_pop, f_pop)
  
  %% Description
  % Specific population growth rate for the hex model.
  % Hazard includes 
  %
  %  * thinning (optional, default: 1; otherwise specified in par.thinning), 
  %  * stage-specific background (optional, default: 0; otherwise specified in par.h_B0b, par.h_Bbj, par.h_Bje, par.h_Bei)
  %  * ageing (controlled by par.h_a and par.s_G)
  %
  % With thinning the hazard rate is such that the feeding rate of a cohort does not change during growth, in absence of other causes of death.
  % Survival of embryo due to ageing is taken for sure
  % Buffer handling rule: constant continuous reproduction is used in the hex model during mean life span as imago.
  % Food density and temperature are assumed to be constant; temperature is specified in par.T_typical.
  % The resulting specific growth rate r is solved from the characteristic equation 1 = \int_0^a_max S(a) R(a) exp(- r a) da
  %   with a_max such that S(a_max) = 1e-6
  %
  % Input
  %
  % * par: structure with parameters for individual (for hazard rates, see remarks)
  % * T_pop: optional temperature (in Kelvin, default C2K(20))
  % * f_pop: optional scalar with scaled functional response (overwrites value in par.f)
  %
  % Output
  %
  % * r: scalar with specific population growth rate
  % * info: scalar with indicator for failure (0) or success (1)
  %
  %% Remarks
  % See <ssd_hex.html *ssd_hex*> for mean age, length, squared length, cubed length and other statistics.
  % See <f_ris0_mod.html *f_ris0_mod*> for f at which r = 0.
  % par.thinning, par.h_B0b, par.h_Bbj, par.h_Bje and par.h_Bei are not standard in structure par; Add them before use if necessary.
  % par.reprodCode is not standard in structure par. Add it before use. If missing, "O" is assumed.

  % unpack par and compute statisitics
  cPar = parscomp_st(par); vars_pull(par);  vars_pull(cPar);  

  % defaults
  if exist('T_pop','var') && ~isempty(T_pop)
    T = T_pop;
  else
    T = C2K(20);
  end
  if exist('f_pop','var') && ~isempty(f_pop)
    f = f_pop;  % overwrites par.f
  end
  if ~exist('thinning','var')
    thinning = 1;
  end
  if ~exist('h_B0b', 'var')
    h_B0b = 0;
  end
  if ~exist('h_Bbj', 'var')
    h_Bbj = 0;
  end
  if ~exist('h_Bje', 'var')
    h_Bje = 0;
  end
  if ~exist('h_Bei', 'var')
    h_Bei = 0;
  end
  if (~exist('reprodCode', 'var') || strcmp(reprodCode, 'O')) && (~exist('genderCode', 'var') || strcmp(genderCode, 'D'))
    kap_R = kap_R/2; % take cost of male production into account
  end
  
  % temperature correction
  pars_T = T_A;
  if exist('T_L','var') && exist('T_AL','var')
    pars_T = [T_A; T_L; T_AL];
  end
  if exist('T_L','var') && exist('T_AL','var') && exist('T_H','var') && exist('T_AH','var')
    pars_T = [T_A; T_L; T_H; T_AL; T_AH]; 
  end
  TC = tempcorr(T, T_ref, pars_T);   % -, Temperature Correction factor
  kT_M = k_M * TC; vT = v * TC; hT_a = h_a * TC^2;   
  
  % supporting statistics
  pars_tj = [g k v_Hb v_He s_j kap kap_V];
  [tau_j, tau_e, tau_b, l_j, l_e, l_b, rho_j, v_Rj, u_Ee, info] = get_tj_hex(pars_tj, f);
  if exist('V_Hx','var')
    pars_UE0 = [V_Hx; g; k_J; k_M; v]; % compose parameter vector
    E_0 = p_Am * initial_scaled_reserve(f, pars_UE0); % J, initial energy
  else
    u_E0 = get_ue0([g k v_Hb], f); E_0 = u_E0 * E_G * L_m^3/ kap; % -, (scaled) cost for egg
  end
  aT_b = tau_b/ kT_M; tT_j = (tau_j - tau_b)/ kT_M; tT_e = (tau_e - tau_b)/ kT_M; % unscale
  L_b = L_m * l_b; L_j = L_m * l_j;  L_e = L_m * l_e;  % unscale
  S_b = exp( - aT_b * h_B0b); % - , survival prob at birth
  rT_j = kT_M * rho_j; % 1/d, growth rate
  
  % life span as imago
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  tau_m = get_tm_s(pars_tm, f, l_b);    % -, scaled mean life span at T_ref
  tT_im = tau_m / kT_M;                 % d, mean life span as imago
  
  % reproduction rate of imago
  E_R_ref = (1 - kap) * E_m * (g + l_b)/ (1 - l_b);  % J/cm^3, reference value for [E_R]
  E_Rj = E_R_ref * s_j * L_j^3;       % J, reproduction buffer
  N = kap_R * E_Rj/ E_0;              % #, number of eggs at emergence
  R = N/ tT_im;                       % #/d, reproduction rate
  tT_N0 = tT_e + tT_im;               % d, time since birth at which all eggs are produced

  pars_qhSC = {f, kT_M, vT, g, k, R, L_b, L_j, L_e, L_m, tT_j, tT_e, tT_N0, rT_j, s_G, hT_a, h_Bbj, h_Bje, h_Bei, thinning};
      
  % find r from char eq 1 = \int_0^infty S(t) R(t) exp(-r*t) dt
  if charEq(0, S_b, pars_qhSC{:}) > 0
    r = NaN; info = 0; % no positive r exists
  else
    [r, info] = nmfzero(@charEq, 0, [], S_b, pars_qhSC{:});
  end
end

% reproduction is continuous
function dqhSC = dget_qhSC(t, qhSC, sgr, f, k_M, v, g, k, R, L_b, L_j, L_e, L_m, t_j, t_e, t_N0, r_j, s_G, h_a, h_Bbj, h_Bje, h_Bei, thinning)
  % t: time since birth
  q   = max(0, qhSC(1)); % 1/d^2, aging acceleration
  h_A = max(0, qhSC(2)); % 1/d^2, hazard rate due to aging
  S   = max(0, qhSC(3)); % -, survival prob
  
  if t < t_j % larva
    h_B = h_Bbj;
    h_X = thinning * r_j; % 1/d, hazard due to thinning
    dq = 0;
    dh_A = 0;
  elseif t < t_e % pupa
    h_B = h_Bje;
    h_X = 0; % 1/d, hazard due to thinning
    dq = 0;
    dh_A = 0;
  else % imago
    h_B = h_Bei;
    s_M = L_j/ L_b;
    L = L_e;
    r = 0; % 1/d, spec growth rate of structure
    h_X = 0; % 1/d, hazard due to thinning
    dq = (q * s_G + h_a) * f * v * s_M/ L;
    dh_A = q; % 1/d^2, change in hazard due to aging
  end

  h = h_A + h_B + h_X; 
  dS = - h * S; % 1/d, change in survival prob
      
  dCharEq = (t > t_e) * (t < t_N0) * S * R  * exp(- sgr * t);
  
  dqhSC = [dq; dh_A; dS; dCharEq]; 
end

function value = charEq (sgr, S_b, f, k_M, v, g, k, R, L_b, L_j, L_e, L_m, t_j, t_e, t_N0, r_j, s_G, h_a, h_Bbj, h_Bje, h_Bei, thinning)
  options = odeset('AbsTol',1e-8, 'RelTol',1e-8);  
  [t, qhSC] = ode45(@dget_qhSC, [0 t_N0], [0 0 S_b 0], options, sgr, f, k_M, v, g, k, R, L_b, L_j, L_e, L_m, t_j, t_e, t_N0, r_j, s_G, h_a, h_Bbj, h_Bje, h_Bei, thinning);
  value = 1 - qhSC(end,4);
end
