read = readmatrix("airfoil_data_NREL5MW\DU21_A17");

len = size(read, 1);
cylinder1 = zeros(len, 4);
cylinder1(:,1) = read(:,1);
cylinder1(:, 3) = 0.5;

cylinder2 = zeros(len, 4);
cylinder1(:,1) = read(:,1);
cylinder2(:, 3) = 0.35;

