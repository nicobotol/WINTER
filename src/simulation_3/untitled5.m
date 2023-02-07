clear s
syms B I p Lambda R L s

yq1 = (B + s*I)/(1.5*(p*Lambda)^2*(1 + 2/3*(R + s*L)*(B + s*I)/(p*Lambda)^2));

yq2 = (B+s*I)/(L*I*s^2+(R*I+L*B)*s+R*B+1.5*(p*Lambda)^2);

simplify(yq1 - yq2)
