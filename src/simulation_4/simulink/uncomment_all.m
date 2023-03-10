function uncomment_all(simulation_mdl)
% Uncomment all the blocks in the model

% Close the Simulink model
save_system(simulation_mdl)
close_system(simulation_mdl);

% Load the Simulink model
load_system(simulation_mdl);

% Get the contents of the model
model_contents = fileread(['simulink\' simulation_mdl '.slx']);

% Define the regular expression pattern to match the commented blocks
pattern = '%\s*(?<block>\w+)\s*=\s*(?<type>\w+)\(';

% Find all the matches of the pattern in the model contents
matches = regexp(model_contents, pattern, 'names');

% Loop through the matches and uncomment the corresponding lines
for i = 1:numel(matches)
    blockName = matches(i).block;
    model_contents = regexprep(model_contents, sprintf('%%\\s*%s', ...
      blockName), blockName);
end

% Save the modified model contents to the model file
fid = fopen([simulation_mdl '.slx'], 'w');
fwrite(fid, model_contents, 'char');
fclose(fid);

% Close the Simulink model
% save_system(simulation_mdl)
close_system(simulation_mdl);

end