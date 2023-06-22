function RMS = RMS_error(x,y)
  % This function calculates the RMS error between two vectors

  RMS = sqrt(mean((x-y).^2));

end