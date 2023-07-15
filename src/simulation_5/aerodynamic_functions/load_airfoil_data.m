function [aoa, cl_mat, cd_mat] = load_airfoil_data_NREL5MW(filenames)
% This function loads the data of the different airfoils and put them into
% matrices. Each column of these matrices corresponds to one specific
% airfoil t/c
% aoa ->  matrix : by columns angle used to define the airfoil
%         parameters
%
%        ATTENTION: ANGLES ARE IMMEDIATELY CONVERTED IN RADIANS
%
% cl_mat -> matrix : by columns cl coefficient for differen angle of
%           attack
% cd_mat -> matrix : by columns cl coefficient for differen angle of
%           attack
% filenames -> vector of strings containing the file to load

% initialization
% aoa = zeros(105,6);
% cl_mat = zeros(105,6);
% cd_mat = zeros(105,6);

mat = cell(1, length(filenames));
aoa_str = cell(1, length(filenames));
cl_str = cell(1, length(filenames));
cd_str = cell(1, length(filenames));

for i = 1:length(filenames)
  mat{i} = readmatrix(filenames(i));
  sizes(i) = [length(mat{i})];
  aoa_str{i}(:) = deg2rad(mat{i}(:, 1)); % conversion from deg to rad
  cl_str{i}(:) = mat{i}(:,2);
  cd_str{i}(:) = mat{i}(:,3);
end

[max_size, pos] = max(sizes);
aoa_tmp = ones(max_size, length(filenames));
aoa = aoa_tmp.*aoa_str{pos}(:);

% interpolate
cl_mat = zeros(max_size, length(filenames));
cd_mat = zeros(max_size, length(filenames));
for i = 1:length(filenames)
  cl_mat(:, i) = interp1(aoa_str{i}(:), cl_str{i}(:)', aoa(:, i));
  cd_mat(:, i) = interp1(aoa_str{i}(:), cd_str{i}(:)', aoa(:, i));
end

end