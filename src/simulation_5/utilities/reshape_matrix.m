function vec = reshape_matrix(mat, x, y)
  %% This function reshapes the value stored in a matrix as function of two parameters in a vectorial form suitable for plots

  len_x = length(x);
  len_y = length(y);

  vec(1, :) = kron(ones(1, len_y), x);
  vec(2, :) = kron(y, ones(1, len_x));
  vec(3, :) = reshape(mat', 1, len_x*len_y);