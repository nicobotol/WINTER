function [r, c, beta, thick] = load_blade_data(blade_filename)
% Function to upload data of the airfoilfrom the .txt file
% r -> radius of the sections (m)
% c -> chord of the section (m)
% beta -> twist angle of the section (deg)
% thick -> t/c ratio
% blade_filaname -> string of the file with the blade data
% opts.CommentStyle = '%';
mat = readmatrix(blade_filename, "CommentStyle", "%");
r(:) = mat(:, 1);
beta(:) = deg2rad(mat(:, 3));
c(:) = mat(:, 2);
thick_ID(:) = mat(:, 4); % identifier of the airfoil profile

% t/c for the different airfoil along the bladespan
% cylinder1 cylinder2 DU40_A17 DU35_A17 DU30_A17 DU25_A17 DU20_A17 NACA64_A17 
thick_table = [100 100 40.5 35.09 30 25 21 18];

thick = thick_table(thick_ID); % t/c

end