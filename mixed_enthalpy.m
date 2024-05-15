function H_diff = mixed_enthalpy(T, w_fracs, H_mix_target)

[H_H2, mu_H2, k_H2, C_H2] = hydrogen(T);
[H_H2O, mu_H2O, k_H2O, C_H2O] = water(T);
[H_N2, mu_N2, k_N2, C_N2] = nitrogen(T);

% Create vectors of properties for prop_mix function
Hs = [H_H2, H_H2O, H_N2];
mus = [mu_H2, mu_H2O, mu_N2];
ks = [k_H2, k_H2O, k_N2];
Cps = [C_H2, C_H2O, C_N2];

% get vicosicty and specific heat of gas mixture in
[H, ~, ~,  ~] = prop_mix(Hs, mus, ks, Cps, w_fracs);

H_diff = H - H_mix_target;