figure()
for i=1:wind.WS_len
  len = 1:1:out_store{i}.P_G.TimeInfo.Length;
  y =  diff(out_store{i}.P_G.Time);
  
  plot(out_store{i}.P_G.Time(1:end-1), y, 'o');
  hold on
  xlabel('Time [s]')
  ylabel('Time step [s]')
end
