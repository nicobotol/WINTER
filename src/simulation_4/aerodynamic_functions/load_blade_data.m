function [r, c, beta, thick] = load_blade_data(blade_filename)
% Function to upload data of the airfoilfrom the .txt file
% r -> radius of the sections (m)
% c -> chord of the section (m)
% beta -> twist angle of the section (deg)
% thick -> t/c ratio
% blade_filaname -> string of the file with the blade data

mat = readmatrix(blade_filename);
r(:) = mat(:, 1);
beta(:) = deg2rad(mat(:, 2));
c(:) = mat(:, 3);
thick(:) = mat(:, 4);

end