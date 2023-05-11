figure();
  
hold on
P = out_store{1}.P_GMech.Data + out_store{1}.P_GJoule.Data + out_store{1}.P_GInductance.Data;
% plot(out_store{1}.P_GE.Time, out_store{1}.P_GSum.Data, '--', 'DisplayName','SUM')
plot(out_store{1}.P_GE.Time, out_store{1}.P_GE.Data, '--', 'DisplayName','GE')
plot(out_store{1}.P_GE.Time, P, '--', 'DisplayName','P')
plot(out_store{1}.P_GE.Time, out_store{1}.P_G.Data, '-.', 'DisplayName','Input')
plot(out_store{1}.P_GE.Time, out_store{1}.P_GInductance.Data, '-.', 'DisplayName','Inductance')
plot(out_store{1}.P_GE.Time, out_store{1}.P_GMech.Data, '-.', 'DisplayName','Mech-to-electric')
legend('location', 'northwest' )
