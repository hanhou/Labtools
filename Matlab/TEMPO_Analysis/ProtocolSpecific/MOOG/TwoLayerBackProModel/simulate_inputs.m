
for i=1:num_unique_points
    
    % this is idiotic, not sure why i did this...
    %fix_az = pi/2 + ga_az;   % because horizontal gaze angles vary from straight ahead
    fix_az = ga_az;  % doing?
    fix_el = ga_el;
    % gaze angles vary hc coords
    [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(unique_point_azimuth_r(i),...
        unique_point_elevation_r(i),-fix_az,-fix_el);
    az_hc = azimuth_gaze; el_hc = elevation_gaze;
    az_ec = unique_point_azimuth_r(i); el_ec = unique_point_elevation_r(i);
%     % gaze angles vary ec coords
%     [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(unique_point_azimuth_r(i),...
%         unique_point_elevation_r(i),fix_az,fix_el);
%     az_hc = unique_point_azimuth_r(i); el_hc = unique_point_elevation_r(i);
%     az_ec = azimuth_gaze; el_ec = elevation_gaze;
    
    ec_inputs = training_inputs(:,beg_i + i);
    
    % normalize to [-1 1]
    pnewn = tramnmxr(norm_rng,ec_inputs,minp,maxp);
    
    %[all_outputs] = msim(net,ec_inputs);
    [all_outputs] = msim(net,pnewn);
    
    %net_outputs = all_outputs(num_hidden_units+1:end);
    anewn = all_outputs(num_hidden_units+1:end);
    
    % map back to original range
    net_outputs = postmnmxr(norm_rng,anewn,mint,maxt);            
    
    hidden_outputs(m,j,i,1:num_hidden_units) = ...
        all_outputs(1:num_hidden_units);

    if plot_training_results %& m > 1 & j > 5
        % for testing training
        hc_outputs = training_outputs(:,beg_i + i);
        basis_model_contour_plots3(unique_azimuth, unique_elevation, ...
            unique_point_azimuth, unique_point_elevation, ec_inputs(1:num_neurons), ...
            ec_inputs(num_neurons+2*num_eye_neurons+1:end), net_outputs, hc_outputs);
        sprintf('ec_az %.2f ec_el %.2f\n gaze_az %.2f gaze_el %.2f\n hc_az %.2f hc_el %.2f\n %s', ...
            az_ec/pi*180, el_ec/pi*180, ga_az/pi*180, ga_el/pi*180, az_hc/pi*180, el_hc/pi*180, curvefit_stim_string{m})
        pause;
    end
end
