my = polyval(coeff_dMdtheta2, theta_v);
DTU10MW_ref = -2.8*(1 + theta_v/164.13 + theta_v.^2/702.09);

(my - DTU10MW_ref)./theta_v

figure()
% plot(theta_v, my - DTU10MW_ref)
% hold on
plot(theta_v, (my - DTU10MW_ref)./theta_v)