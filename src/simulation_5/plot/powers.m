figure();
  
hold on
P = out_store{1}.P_GE.Data + out_store{1}.P_GJoule.Data + out_store{1}.P_GInductance.Data;
% plot(out_store{1}.P_GE.Time, out_store{1}.P_GSum.Data, '--', 'DisplayName','SUM')
plot(out_store{1}.P_GE.Time, out_store{1}.P_GE.Data, '--', 'DisplayName','GE', 'LineWidth',line_width)
plot(out_store{1}.P_GE.Time, out_store{1}.P_GJoule.Data, '--', 'DisplayName','Joule', 'LineWidth',line_width)
plot(out_store{1}.P_GE.Time, P, '--', 'DisplayName','GE + Inductance + Joule', 'LineWidth',line_width, 'Marker','o')
plot(out_store{1}.P_GE.Time, out_store{1}.P_G.Data, '-.', 'DisplayName','Input', 'LineWidth',line_width)
% plot(out_store{1}.P_GE.Time, out_store{1}.P_GInductance.Data, '-.', 'DisplayName','Inductance')
% plot(out_store{1}.P_GE.Time, out_store{1}.P_GMech.Data, '-.', 'DisplayName','Mech-to-electric', 'LineWidth',line_width)
legend('location', 'northwest' )
grid on
xlabel('Time [s]')
ylabel('Power [W]')