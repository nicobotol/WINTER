function [cl, cd] = interpolate_cl_cd(aoa, cl_mat, cd_mat, thick_prof, alpha, thick)
% This function interpolates the data of cd and cl from the .txt file.
% cl -> cl coefficient
% cd -> cd coefficient
% aoa -> angle of attack (loaded from the .txt file)
% cl_mat -> cl matrix (loaded from the .txt file)
% cd_mat -> cd matrix (loaded from the .txt file)
% thick_prof -> vector of thickness
% alpha -> angle of attack (rad)
% thick -> t/c ratio

% initialize the vectors
clthick = zeros(1,6);
cdthick = zeros(1,6);

%interpolate the values to the different thicknesses
for k=1:6 % k indicate the airfoil
  clthick(k) = interp1(aoa(:,k),cl_mat(:,k), alpha);
  cdthick(k) = interp1(aoa(:,k),cd_mat(:,k), alpha);
end

% then interpolate to the actual thickness
% thick_prof =(100,60,48,36,30.1,24.1), i indicates the element nr .

cl = interp1(thick_prof(:), clthick(:), thick);
cd = interp1(thick_prof(:), cdthick(:), thick);

end