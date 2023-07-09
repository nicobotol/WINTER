%% Define new colors
function col = color(name)
switch name
    case {'green' , 1}
        col = [25   125 44] ./ 255;
    case {'orange', 2}
        col = [239  107 68] ./ 255;
    case {'blue' , 3}
        col = [65   46  248] ./ 255;
    case {'red', 4}
        col = [214  0   18] ./ 255;
    case {'grey', 5}
        col = [0.7  0.7 0.7];
    case {'indian_red', 6}
        col = [176   23  31]./ 255;
    case {'violet_red', 7}
        col = [208	32	144]./ 255;
    case {'orchid', 8}
        col = [139	71	137]./ 255;
    case {'purple', 9}
        col = [128	0	128]./ 255;
    case {'indigo', 10}
        col = [75	0	130]./ 255;
    case {'midnight_blue', 11}
        col = [25	25	112]./ 255;
    case {'dodger_blue', 12}
        col = [30	144	255]./ 255;
    case {'deepsky_blue', 13}
        col = [0	191	255]./ 255;
    case {'turquise', 14}
        col = [0	229	238]./ 255;
    case {'darkturquise', 15}
        col = [0	206	209]./ 255;
    case {'darkslate_grey', 16}
        col = [47	79	79]./ 255;
    case {'aquamarine', 17}
        col = [69	139	116]./ 255;
    case {'cobalt_green', 18}
        col = [61	145	64]./ 255;
    case {'dark_green', 19}
        col = [0	100	0]./ 255;
    case {'yellow', 20}
        col = [205	205	0]./ 255;
    case {'gold', 21}
        col = [238	201	0]./ 255;
    case {'comsilk', 22}
        col = [139	136	120]./ 255;
    otherwise
        warning('Color name not recognized. Color set to black.')
        col = [0	0	0]./ 255;
end
end