function [I_blade, I_rotor, m_blade, m_rotor] = rotor_inertia(r, rho, name)

I_blade = 0;  % initialize blade inertia
I_rotor = 0;  % initialize rotor inertia
dI = 0;       % initialize infinitesimal inertia [kgm^2]
dm = 0;       % initialize infinitesimal mass [kg]
m_blade = 0;  % initialize blade mass

for i = 1:length(r) - 1
  dm = (rho(i) + rho(i + 1))*(r(i + 1) - r(i))/2; % infinitesimal mass [kg]
  m_blade = m_blade + dm;
  dI = dm*((r(i) + r(i + 1))/2)^2; % inertia of the section [kgm^2]
  I_blade = I_blade + dI;
end

I_rotor = 3*I_blade;
m_rotor = 3*m_blade;

fprintf('%s: The blade mass is m = %6.4e [kg]\n',name, m_blade);
fprintf('%s: The blade inertia is I = %6.4e [kgm^2]\n', name, I_blade);
fprintf('%s: The rotor mass is m = %6.4e [kg]\n',name, m_rotor);
fprintf('%s: The rotor inertia is I = %6.4e [kgm^2]\n', name, I_rotor);
end