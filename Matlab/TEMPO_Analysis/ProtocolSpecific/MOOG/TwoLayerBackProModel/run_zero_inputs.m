
range_az = [-pi pi];  % does nothing
range_el = [-pi/2 pi/2];  % does nothing
range_gaze_az = 2*[-pi/9 pi/9];
range_gaze_el = 2*[-pi/9 pi/9];

sprintf('making zero training set...')
[P num_neurons P_eye num_eye_neurons] = init_params(range_az,range_el,range_gaze_az,range_gaze_el);
[training_zero_inputs] = make_training_set_zero(@Curvefit_cos_tuning_5p,@Curvefit_linear_tuning_2p,P,P_eye,...
    range_az,range_el,range_gaze_az,range_gaze_el);
sprintf('done making zero training set...')

sprintf('simulating network with only gaze angle inputs...')
hidden_outputs_zero = zeros(num_gaze_angles,num_hidden_units);
for j=1:num_gaze_angles            
    ga_az = unique_gaze_angles_xr(j);
    ga_el = unique_gaze_angles_yr(j);
        
    zero_inputs = training_zero_inputs(:,j);
    
    % normalize to [-1 1]
    pnewn = tramnmxr(norm_rng,zero_inputs,minp,maxp);

    % use neural network toolbox to simulate inputs on network
    [all_outputs] = msim(net,pnewn);
    
    anewn = all_outputs(num_hidden_units+1:end);
    
    % map back to original range
    zero_outputs = postmnmxr(norm_rng,anewn,mint,maxt);            
    
    hidden_outputs_zero(j,1:num_hidden_units) = all_outputs(1:num_hidden_units,1)';
end % for each sampled gaze angle            
sprintf('done simulating network with only gaze angle inputs...')

