figure()
hold on
for i = 1:length(pitch_vector)
plot3(pitch_vector(i)*ones(1,length(lambda_vector)), lambda_vector, lookup_cP(i,:), 'o')
end

figure()
plot(lambda_vector, lookup_cP(theta_pos, :))