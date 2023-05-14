function [] = set_simulink_parameters(simulation_mdl, simulation_type)
% Set the parameters for the simulink simulation

% Decrease the simulation time for the electrical generator
set_param(strcat(simulation_mdl,'/PMSM/Gain'),'SampleTime', ...
  'simulation.time_step_L');
set_param(strcat(simulation_mdl,'/PMSM/Gain1'),'SampleTime', ...
  'simulation.time_step_L');
set_param(strcat(simulation_mdl,'/PMSM/Sum1'),'SampleTime', ...
  'simulation.time_step_L');

% In case of generator step response set also other parameters
if simulation_type == 4
  % Get all blocks in the current model
  % blocks = find_system(simulation.mdl, 'Type', 'Block');
  blocks = find_system(simulation_mdl,  'MatchFilter', @nonInOutBlocks);

  % Find the index of the element that you want to remove
  indexToRemove = find(strcmp(blocks, simulation_mdl));

  % Check whether the element is present
  if ~isempty(indexToRemove)
      % Remove the element using setdiff
      blocks = setdiff(blocks, simulation_mdl);
  end
  
  % Comment out each block
  for i = 1:length(blocks)
      set_param(blocks{i}, 'Commented', 'on');
  end
  
  % Uncomment blocks necessary for the step response
  set_param(strcat(simulation_mdl, '/PMSM'), 'Commented', 'off');
  block_step = strings(1, 13);
  block_step = ["Step", "Sum", "high_to_low_time", "Gain1", "Sum1", ...
    "Riq","Gc", "Yq", "Gain", "low_to_high_time", "To Workspace", ...
    "To Workspace1", "Scope"];
  
  for i=1:length(block_step)
    set_param(strcat(simulation_mdl, '/PMSM/', block_step(1, i)), ...
      'Commented', 'off');
  end

end
  

end

function match = nonInOutBlocks(handle)
  match = true;
  if strcmp(get_param(handle, 'Type'), 'block')
    blockType = get_param(handle, 'BlockType');
    if strcmp(blockType, 'Inport') || ...
          strcmp(blockType,  'Outport') 
        match = false;
    end
  end
end
