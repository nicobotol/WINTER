
function export_figure(fig, file_name, folder_name, v_size) 

  figure_folder_name = folder_name;
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  
  if nargin == 2
      v_size = 4;
  end
  
  status = mkdir(figure_folder_name);
  status = mkdir(fullfile(figure_folder_name, "vectorial"));
  status = mkdir(fullfile(figure_folder_name, "presentation"));
  status = mkdir(fullfile(figure_folder_name, "matlab"));
  
  export_fig(fig, fullfile(figure_folder_name, "vectorial" ,file_name));

%   file_name = char(file_name);
%   file_name = string(file_name(1:end-4)) + ".emf";
%   export_fig(fig, fullfile(figure_folder_name, "presentation" ,file_name));
%   file_name = char(file_name);
%   file_name = string(file_name(1:end-4)) + ".fig";
%   savefig(fig, fullfile(figure_folder_name, "matlab", file_name));

end
