figure();hold on;grid on;
for i=1:IMM.n_models
  data = out_store{1}.model_x_est.Data(1,i,:);
  data = reshape(data, 1, size(data,3));
  plot(data, 'LineWidth', 2, 'DisplayName', ['Model ', num2str(i)]);
end
plot(out_store{1}.omega_R.Data, 'LineWidth', 2, 'DisplayName', 'omega_R');
plot(out_store{1}.omega_tilde.Data, 'LineWidth', 2, 'DisplayName', 'omega_tilde');
legend()
title('$\omega$')

figure();hold on;grid on;
for i=1:IMM.n_models
  data = out_store{1}.model_x_est.Data(2,i,:);
  data = reshape(data, 1, size(data,3));
  plot(data, 'LineWidth', 2, 'DisplayName', ['Model ', num2str(i)]);
end
plot(out_store{1}.wind.Data, 'LineWidth', 2, 'DisplayName', 'wind');
plot(out_store{1}.V0_tilde.Data, 'LineWidth', 2, 'DisplayName', 'V0_tilde');
legend()
title('$V_0$')

figure();grid on;
plot(out_store{1}.K_opt.Data, 'LineWidth',2)

figure();hold on;grid on;
for i=1:IMM.n_models
  data = out_store{1}.mu.Data(i,1,:);
  data = reshape(data, 1, size(data,3));
  plot(data, 'LineWidth', 2, 'DisplayName', ['Model ', num2str(i)]);
end
legend()
title('$Probability$')