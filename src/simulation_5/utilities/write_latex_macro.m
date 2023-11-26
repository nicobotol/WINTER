function write_latex_macro(file_name, name, data, mode)

  fileID = fopen(file_name,mode);
  formatSpec = '\\def\\%s{%5.1e}\n';
  for i=1:length(data)
    fprintf(fileID,formatSpec, name(i), data(i));
  end
  fclose(fileID);
  
  end