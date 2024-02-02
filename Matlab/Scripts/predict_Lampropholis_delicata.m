function [prdData, info] = predict_Lampropholis_delicata(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

%% compute temperature correction factors
pars_T = [T_A; T_L; T_H; T_AL; T_AH];
TC_ab = tempcorr(temp.ab, T_ref, pars_T);
TC_tp = tempcorr(temp.tp, T_ref, pars_T);
TC_am = tempcorr(temp.am, T_ref, pars_T);
TC_Ri = tempcorr(temp.Ri, T_ref, pars_T);
TC_tL = tempcorr(temp.tL, T_ref, pars_T);
% TC_WO = tempcorr(temp.WJO, T_ref, pars_T); 
TC_tW = tempcorr(temp.tW, T_ref, pars_T); % 
%TC_Tab = tempcorr(C2K(Tab(:,1)), T_ref, pars_T);
  
%% zero-variate data
% life cycle
pars_tp = [g; k; l_T; v_Hb; v_Hp];               % compose parameter vector
[tau_p, tau_b, l_p, l_b, info] = get_tp(pars_tp, f); % -, scaled times & lengths at f

% birth
L_b = L_m * l_b;                  % cm, structural length at birth at f
Lw_b = L_b/ del_M;                % cm, physical length at birth at f
Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f
aT_b = tau_b/ k_M/ TC_ab;           % d, age at birth at f and T

% puberty 
L_p = L_m * l_p;                  % cm, structural length at puberty at f
Lw_p = L_p/ del_M;                % cm, physical length at puberty at f
Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty
tT_p = (tau_p - tau_b)/ k_M/ TC_tp;   % d, time since birth at puberty at f and T

% ultimate
l_i = f - l_T;                    % -, scaled ultimate length
L_i = L_m * l_i;                  % cm, ultimate structural length at f
Lw_i = L_i/ del_M;                % cm, ultimate physical length at f
Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight; egg calculation 0.2466 = clutch 30*3/365

% reproduction
pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
RT_i = TC_Ri * reprod_rate(L_i, f, pars_R);             % #/d, ultimate reproduction rate at T

% life span
pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
%% pack to output
prdData.ab = aT_b;
prdData.tp = tT_p;
prdData.am = aT_m;
prdData.Lb = Lw_b;
prdData.Lp = Lw_p;
prdData.Li = Lw_i;
prdData.Wwb = Ww_b;
prdData.Wwp = Ww_p;
prdData.Wwi = Ww_i;
prdData.Ri = RT_i;
  
%% uni-variate data
% temperature-development time
%aT_b2 = (tau_b/ k_M)./ TC_Tab;           % d, age at birth at f and T

% length-weight
EWw = (LW(:,1) * del_M).^3 * (1 + f * w); % g, expected wet weight at length

% time-length
rT_B = TC_tL * k_M/ 3/ (1 + f/ g); % 1/d, von Bert growth rate
L = L_i - (L_i - L_b) * exp(- rT_B * tL(:,1));
Et_L = L/del_M; % cm, length

% time-weight
L = L_i - (L_i - L_b) * exp(- rT_B * tW(:,1));
Et_W = L.^3 * (1 + E_m /mu_E * w_E/d_E);

%  % O2 consumption rate in ml/min 
% n_M(:,4)=[1, 4/5, 3/5, 4/5]; % uric acid
% O2M = (- n_M\n_O)'; % -, matrix that converts organic to mineral fluxes  O2M is prepared for post-multiplication eq. 4.35
% p_ref = p_Am * L_m^2;               % J/d, max assimilation power at max size
% pars_power = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp];
% L = (WJO(:,1) / (1 + f * w)).^(1/3); % cm, structural length
% pACSJGRD = p_ref * scaled_power(L, f, pars_power, l_b, l_p);
% pADG = pACSJGRD(:, [1 7 5]);  pADG(:,1) = 0; % exclude assim contribution
% JM = pADG * eta_O' * O2M; 
% X_gas = 24.06;  % gas correction factor, L/mol
% EO = (- 1 * (JM(:,3) * X_gas) .* TC_WO)/24*1000; %/ 60; % L/day to ml/min


%% pack to output
%prdData.Tab = aT_b2;
prdData.LW = EWw;
prdData.tL = Et_L;
prdData.tW = Et_W;
% prdData.WJO = EO;
