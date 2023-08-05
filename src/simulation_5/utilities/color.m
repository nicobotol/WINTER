%% Define new colors
function col = color(name)
  switch name
      case {1}
          col = [0 0.4470 0.7410];
      case {2}
          col = [0.8500 0.3250 0.0980];
      case {2}
          col = [0.9290 0.6940 0.1250];
      case {4}
          col = [0.4940 0.1840 0.5560];
      case {5}
          col = [0.4660 0.6740 0.1880];
      case {6}
          col = [0.3010 0.7450 0.9330];
      case {7}
          col = [0.6350 0.0780 0.1840];
      case {'green' , 8}
          col = [25   125 44] ./ 255;
      case {'orange', 9}
          col = [239  107 68] ./ 255;
      case {'blue' , 10}
          col = [65   46  248] ./ 255;
      case {'red', 11}
          col = [214  0   18] ./ 255;
      case {'grey', 12}
          col = [0.7  0.7 0.7];
      case {'indian_red', 13}
          col = [176   23  31]./ 255;
      case {'violet_red', 14}
          col = [208	32	144]./ 255;
      case {'orchid', 15}
          col = [139	71	137]./ 255;
      case {'purple', 16}
          col = [128	0	128]./ 255;
      case {'indigo', 17}
          col = [75	0	130]./ 255;
      case {'midnight_blue', 18}
          col = [25	25	112]./ 255;
      case {'dodger_blue', 19}
          col = [30	144	255]./ 255;
      case {'deepsky_blue', 20}
          col = [0	191	255]./ 255;
      case {'turquise', 21}
          col = [0	229	238]./ 255;
      case {'darkturquise', 22}
          col = [0	206	209]./ 255;
      case {'darkslate_grey', 23}
          col = [47	79	79]./ 255;
      case {'aquamarine', 24}
          col = [69	139	116]./ 255;
      case {'cobalt_green', 25}
          col = [61	145	64]./ 255;
      case {'dark_green', 26}
          col = [0	100	0]./ 255;
      case {'yellow', 27}
          col = [205	205	0]./ 255;
      case {'gold', 28}
          col = [238	201	0]./ 255;
      case {'comsilk', 29}
          col = [139	136	120]./ 255;
      otherwise
          warning('Color name not recognized. Color set to black.')
          col = [0	0	0]./ 255;
  end
  end