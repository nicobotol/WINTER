function [aoa, cl_mat, cd_mat] = load_airfoil_data(filenames)
% This function loads the data of the different airfoils and put them into
% matrices. Each column of these matrices corresponds to one specific
% airfoil t/c
% aoa ->  matrix (105x6): by columns angle used to define the airfoil
%         parameters
%
%        ATTENTION: ANGLES ARE IMMEDIATELY CONVERTED IN RADIANS
%
% cl_mat -> matrix (105x6): by columns cl coefficient for differen angle of
%           attack
% cd_mat -> matrix (105x6): by columns cl coefficient for differen angle of
%           attack
% filenames -> vector of strings containing the file to load

% initialization
aoa = zeros(105,6);
cl_mat = zeros(105,6);
cd_mat = zeros(105,6);

for i = 1:6
  mat = readmatrix(filenames(i));
  aoa(:,i) = deg2rad(mat(:,1)); % conversion from deg to rad
  cl_mat(:,i) = mat(:,2);
  cd_mat(:,i) = mat(:,3);
end

end