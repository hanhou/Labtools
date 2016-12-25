
Curvefit_defines;

remake_training_set = 1;
retrain_network = 1;
simulate_network = 1;
plot_training_results = 0;
plot_hidden_layer = 0;

% this code could stand a rewrite (have fun)

if remake_training_set
    range_az = [-pi pi];  % does nothing
    range_el = [-pi/2 pi/2];  % does nothing
    %range_gaze_az = [-pi/4 pi/4];
    %range_gaze_el = [-pi/4 pi/4];
    range_gaze_az = 2*[-pi/9 pi/9];
    range_gaze_el = 2*[-pi/9 pi/9];
    
    global num_unique_elevation;  % this is a hack to set this param from a top level looping script
    global num_unique_elevation_load;
    if ~isempty(num_unique_elevation) & num_unique_elevation == 0
        num_unique_elevation = num_unique_elevation_load;
    else
        % num_unique_elevation = 9;  % number of latitudes to generate (including poles)
        num_unique_elevation = 5;  % number of latitudes to generate (including poles)
    end
    global num_unique_azimuth;  % this is a hack to set this param from a top level looping script
    global num_unique_azimuth_load;
    if ~isempty(num_unique_azimuth) & num_unique_azimuth == 0
        num_unique_azimuth = num_unique_azimuth_load;
    else
        % num_unique_azimuth = 16;    % number of longitutdes to generate
        num_unique_azimuth = 8;    % number of longitutdes to generate
    end
    
    sprintf('making training set...')
    [Pvis Pves num_neurons P_eye num_eye_neurons] = init_params(range_az,range_el,range_gaze_az,range_gaze_el);
%     [training_inputs training_outputs] = make_training_set(@Curvefit_cos_tuning_5p,@Curvefit_gaus1_tuning_3p,P,P_eye,...
%         range_az,range_el,range_gaze_az,range_gaze_el);
    [training_inputs training_outputs] = make_training_set(@Curvefit_cos_tuning_5p,@Curvefit_linear_tuning_2p,Pvis,Pves,P_eye,...
        range_az,range_el,range_gaze_az,range_gaze_el);
    sprintf('done making training set...')
end

%return

if retrain_network
    norm_rng = [-1 1];  % range that the inputs will be normalized to
    %input_rng = [repmat([0 35],num_neurons,1); repmat([0 15],num_eye_neurons,1); repmat([0 35],num_neurons,1)];  %gaussian
    %input_rng = [repmat([0 35],num_neurons,1); repmat([5 25],num_eye_neurons,1); repmat([0 35],num_neurons,1)];  %linear
    input_rng = repmat(norm_rng,2*num_neurons+2*num_eye_neurons,1);
    global num_hidden_units;  % this is a hack to set this param from a top level looping script
    global num_hidden_units_load;
    if ~isempty(num_hidden_units) & num_hidden_units == 0
        num_hidden_units = num_hidden_units_load;
    else
        %num_hidden_units = num_neurons*num_eye_neurons;
        num_hidden_units = 150;
    end
%     training_alg = 'trainbfg';  % out of mem
%     training_alg = 'trainbr';   % out of mem
%%    training_alg = 'traincgb';  % 2.75e-3
%     training_alg = 'traincgf';  % 2.76e-3
%     training_alg = 'traincgp';  % 2.78e-3
%     training_alg = 'traingd';   % 8.57e-2
%     training_alg = 'traingda';  % 9.32e-2
%     training_alg = 'traingdm';  % 8.88e-2
%     training_alg = 'traingdx';  % 7.34e-2
%     training_alg = 'trainlm';   % out of mem
%     training_alg = 'trainoss';  % 2.81e-3
%     training_alg = 'trainrp';   % 4.78e-2
    training_alg = 'trainscg';  % 2.83e-3
    net = newff(input_rng,[num_hidden_units,num_neurons],...
        {'tansig','purelin'},training_alg);
    net = init(net);
    
    net.trainParam.show = 50;
    net.trainParam.lr = 0.1;
    net.trainParam.lr_inc = 1.05;
    net.trainParam.mc = 0.8;
    net.trainParam.epochs = 1000;
    net.trainParam.goal = 1e-5;
    net.trainParam.srchFcn   = 'srchcha';
    net.trainParam.searchFcn = 'srchcha';
%     net.trainParam.srchFcn   = 'srchgol';
%     net.trainParam.searchFcn = 'srchgol';
%     net.trainParam.srchFcn   = 'srchbre';
%     net.trainParam.searchFcn = 'srchbre';
%     net.trainParam.srchFcn   = 'srchhyb';
%     net.trainParam.searchFcn = 'srchhyb';
%     net.trainParam.srchFcn   = 'srchbac';
%     net.trainParam.searchFcn = 'srchbac';

    % use weights as a portion of the performance to minimize weightages
    net.performFcn = 'msereg';
    net.performParam.ratio = 0.5;

    % normalize to some range
    [pn,minp,maxp,tn,mint,maxt] = premnmxr(norm_rng,training_inputs,training_outputs);

    %[net,tr]=train(net,training_inputs,training_outputs);
    [net,tr]=train(net,pn,tn);
    
%     % radial basis network
%     goal = 0.1; spread = 5;
%     [net,tr] = newrb(pn,tn,goal,spread);
%     num_hidden_units = net.inputWeights{1,1}.size(1);
end

create_sample_points;

num_stim_types = 3;

if simulate_network
    sprintf('simulating network...')
    hidden_outputs = zeros(num_stim_types,num_gaze_angles,num_unique_points,num_hidden_units);
    for m=1:num_stim_types
        for j=1:num_gaze_angles            
            % start and stop indices for this trial set
            beg_i = (m-1)*num_gaze_angles*num_unique_points + (j-1)*num_unique_points;
            end_i = (m-1)*num_gaze_angles*num_unique_points + j*num_unique_points;
            
            ga_az = unique_gaze_angles_xr(j);
            ga_el = unique_gaze_angles_yr(j);
            
            simulate_inputs;        
        end % for each sampled gaze angle            
    end % for each stim type
    sprintf('done simulating network...')
end

% calculate errors
pnewn = tramnmxr(norm_rng,training_inputs,minp,maxp);  % normalize
[all_outputs] = msim(net,pnewn);
anewn = all_outputs(num_hidden_units+1:end,:);
net_outputs = postmnmxr(norm_rng,anewn,mint,maxt);  % back to original range
onewn = tramnmxr(norm_rng,training_outputs,mint,maxt);  % normalize
network_norm_mse = mse(anewn-onewn);
network_mse = mse(net_outputs-training_outputs);
network_norm_msereg = msereg(anewn-onewn,net);
network_msereg = msereg(net_outputs-training_outputs,net);
% network_norm_msereg = 0;
% network_msereg = 0;

% compare preferred and max directions for vestibular only and visual only conditions
pref_dir_visual_az = zeros(1,num_hidden_units); pref_dir_visual_el = zeros(1,num_hidden_units);
pref_dir_vestib_az = zeros(1,num_hidden_units); pref_dir_vestib_el = zeros(1,num_hidden_units);
pref_dir_combin_az = zeros(1,num_hidden_units); pref_dir_combin_el = zeros(1,num_hidden_units);
pref_n45_dir_visual_az = zeros(1,num_hidden_units); pref_n45_dir_visual_el = zeros(1,num_hidden_units);
pref_n45_dir_vestib_az = zeros(1,num_hidden_units); pref_n45_dir_vestib_el = zeros(1,num_hidden_units);
pref_n45_dir_combin_az = zeros(1,num_hidden_units); pref_n45_dir_combin_el = zeros(1,num_hidden_units);
pref_visual_mag = zeros(1,num_hidden_units);
pref_vestib_mag = zeros(1,num_hidden_units);
pref_combin_mag = zeros(1,num_hidden_units);
max_dir_visual_az = zeros(1,num_hidden_units); max_dir_visual_el = zeros(1,num_hidden_units);
max_dir_vestib_az = zeros(1,num_hidden_units); max_dir_vestib_el = zeros(1,num_hidden_units);
max_dir_combin_az = zeros(1,num_hidden_units); max_dir_combin_el = zeros(1,num_hidden_units);
max_visual_mag = zeros(1,num_hidden_units);
max_vestib_mag = zeros(1,num_hidden_units);
max_combin_mag = zeros(1,num_hidden_units);

% normalized dot products in response space.
ndot_visual_vestib = zeros(1,num_hidden_units);
ndot_visual_combin = zeros(1,num_hidden_units);
ndot_vestib_combin = zeros(1,num_hidden_units);

% compute heading tuning index for each hidden unit.
% HTI = norm of vector sums / sum of vector norms
HTI_vestib = zeros(1,num_hidden_units); 
HTI_visual = zeros(1,num_hidden_units);
HTI_combin = zeros(1,num_hidden_units);

regr_b = zeros(2,num_hidden_units);

ga1_out = zeros(1,num_unique_points);
ga2_out = zeros(1,num_unique_points);
ga3_out = zeros(1,num_unique_points);
ga4_out = zeros(1,num_unique_points);
ga5_out = zeros(1,num_unique_points);
ga6_out = zeros(1,num_unique_points);
ga7_out = zeros(1,num_unique_points);
ga8_out = zeros(1,num_unique_points);
ga9_out = zeros(1,num_unique_points);
ga10_out = zeros(1,num_unique_points);
str = {'vary horz' 'vary vert'};
for i=1:num_hidden_units

    % these should be renamed
    ga1_out(1:num_unique_points) = hidden_outputs(CURVEFIT_VESTIBULAR_STIM,3,:,i);
    ga2_out(1:num_unique_points) = hidden_outputs(CURVEFIT_OCCULAR_STIM,3,:,i);
    ga3_out(1:num_unique_points) = hidden_outputs(CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM,3,:,i);
    ga7_out(1:num_unique_points) = hidden_outputs(CURVEFIT_VESTIBULAR_STIM,1,:,i);
    ga8_out(1:num_unique_points) = hidden_outputs(CURVEFIT_OCCULAR_STIM,1,:,i);
    ga9_out(1:num_unique_points) = hidden_outputs(CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM,1,:,i);
    
    % horz and vert at gaze angle 0 should be the same.
    ga4_out(1:num_unique_points) = hidden_outputs(CURVEFIT_VESTIBULAR_STIM,8,:,i);
    ga5_out(1:num_unique_points) = hidden_outputs(CURVEFIT_OCCULAR_STIM,8,:,i);
    ga6_out(1:num_unique_points) = hidden_outputs(CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM,8,:,i);
    %if any(ga1_out ~= ga4_out) | any(ga2_out ~= ga5_out) | any(ga3_out ~= ga6_out)
    %    error('the thing that should not be.');
    %end
    
    % this transformation effects max response points and HTI calculation,
    % but does not change the preferred direction vectors.  this considers
    % the min response of the cell to be the spontaneous firing rate.
    ga1_out = ga1_out - min(ga1_out);
    ga2_out = ga2_out - min(ga2_out);
    ga3_out = ga3_out - min(ga3_out);
%     ga1_out = ga1_out + 1;
%     ga2_out = ga2_out + 1;
%     ga3_out = ga3_out + 1;
    ga7_out = ga7_out - min(ga7_out);
    ga8_out = ga8_out - min(ga8_out);
    ga9_out = ga9_out - min(ga9_out);

    ndot_visual_vestib(i) = sum(ga1_out.*ga2_out)/norm(ga1_out)/norm(ga2_out);
    ndot_visual_combin(i) = sum(ga1_out.*ga3_out)/norm(ga1_out)/norm(ga3_out);
    ndot_vestib_combin(i) = sum(ga2_out.*ga3_out)/norm(ga2_out)/norm(ga3_out);

    [x y z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, ga1_out);
    x = sum(x); y = sum(y); z = sum(z);
    [pref_dir_vestib_az(i) pref_dir_vestib_el(i) pref_vestib_mag(i)] = cart2sph(x,y,z);
    HTI_vestib(i) = pref_vestib_mag(i) / sum(abs(ga1_out));  % norm of vector sums divided by sum of vector norms
    
    [x y z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, ga2_out);
    x = sum(x); y = sum(y); z = sum(z);
    [pref_dir_visual_az(i) pref_dir_visual_el(i) pref_visual_mag(i)] = cart2sph(x,y,z);
    HTI_visual(i) = pref_visual_mag(i) / sum(abs(ga2_out));  % norm of vector sums divided by sum of vector norms

    [x y z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, ga3_out);
    x = sum(x); y = sum(y); z = sum(z);
    [pref_dir_combin_az(i) pref_dir_combin_el(i) pref_combin_mag(i)] = cart2sph(x,y,z);
    HTI_combin(i) = pref_combin_mag(i) / sum(abs(ga3_out));  % norm of vector sums divided by sum of vector norms

    [x y z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, ga7_out);
    x = sum(x); y = sum(y); z = sum(z);
    [pref_n45_dir_vestib_az(i) pref_n45_dir_vestib_el(i) r] = cart2sph(x,y,z);
    
    [x y z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, ga8_out);
    x = sum(x); y = sum(y); z = sum(z);
    [pref_n45_dir_visual_az(i) pref_n45_dir_visual_el(i) r] = cart2sph(x,y,z);

    [x y z] = sph2cart(unique_point_azimuth_r, unique_point_elevation_r, ga9_out);
    x = sum(x); y = sum(y); z = sum(z);
    [pref_n45_dir_combin_az(i) pref_n45_dir_combin_el(i) r] = cart2sph(x,y,z);
        
    sel = (ga1_out == max(ga1_out));
    if sum(sel) == 1
        max_dir_vestib_az(i) = unique_point_azimuth_r(sel);
        max_dir_vestib_el(i) = unique_point_elevation_r(sel);
        max_vestib_mag(i) = max(ga1_out);
    end
    sel = (ga2_out == max(ga2_out));
    if sum(sel) == 1
        max_dir_visual_az(i) = unique_point_azimuth_r(sel);
        max_dir_visual_el(i) = unique_point_elevation_r(sel);
        max_visual_mag(i) = max(ga2_out);
    end
    sel = (ga3_out == max(ga3_out));
    if sum(sel) == 1
        max_dir_combin_az(i) = unique_point_azimuth_r(sel);
        max_dir_combin_el(i) = unique_point_elevation_r(sel);
        max_combin_mag(i) = max(ga3_out);
    end

    fit_resp = ga3_out - ga2_out;
    % subtracting the min response does not effect the slope term 
    % from the regression (just the offset).
    %fit_resp = fit_resp - min(fit_resp);
    regr_combined_minus_visual = [fit_resp' ones(num_unique_points,1)];
    regr_vestib = ga1_out';
    [regr_b(1:2,i) b_int] = regress(regr_vestib, regr_combined_minus_visual, 0.05);

%     tmp = [ga1_out' ga2_out' ones(num_unique_points,1)];
%     [regr_b(1:3,i) b_int(1:3,1:2,i)] = regress(ga3_out', tmp, 0.05);
    
    if plot_hidden_layer %& i == 57
        for m=1:num_stim_types
            if m == CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM
                pref_az = pref_dir_combin_az(i)/pi*180;
                pref_el = pref_dir_combin_el(i)/pi*180;
                max_az = max_dir_combin_az(i)/pi*180;
                max_el = max_dir_combin_el(i)/pi*180;
                HTI = HTI_combin(i);
                sel = pref_az < -90; pref_az(sel) = pref_az(sel) + 360;
                sel = max_az < -90; max_az(sel) = max_az(sel) + 360;
                for j=1:2
                    ga1_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+1,:,i);
                    ga2_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+2,:,i);
                    ga3_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+3,:,i);
                    ga4_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+4,:,i);
                    ga5_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+5,:,i);
                    ga6_out(1:num_unique_points) = hidden_outputs(1,(j-1)*5+1,:,i) + hidden_outputs(2,(j-1)*5+1,:,i);
                    ga7_out(1:num_unique_points) = hidden_outputs(1,(j-1)*5+2,:,i) + hidden_outputs(2,(j-1)*5+2,:,i);
                    ga8_out(1:num_unique_points) = hidden_outputs(1,(j-1)*5+3,:,i) + hidden_outputs(2,(j-1)*5+3,:,i);
                    ga9_out(1:num_unique_points) = hidden_outputs(1,(j-1)*5+4,:,i) + hidden_outputs(2,(j-1)*5+4,:,i);
                    ga10_out(1:num_unique_points) = hidden_outputs(1,(j-1)*5+5,:,i) + hidden_outputs(2,(j-1)*5+5,:,i);
                    basis_model_contour_plots4(unique_azimuth, unique_elevation, ...
                        unique_point_azimuth, unique_point_elevation, ...
                        ga1_out, ga2_out, ga3_out, ga4_out, ga5_out, ...
                        ga6_out, ga7_out, ga8_out, ga9_out, ga10_out, (m-1)*num_stim_types + j, ...
                        sprintf('cell %d - %s - %s - pref %.2f %.2f - max %.2f %.2f - HTI %.3f', ...
                        i, str{j}, curvefit_stim_string{m}, pref_az, pref_el, max_az, max_el, HTI) );
                end            
            else
                if m == CURVEFIT_OCCULAR_STIM
                    pref_az = pref_dir_visual_az(i)/pi*180;
                    pref_el = pref_dir_visual_el(i)/pi*180;
                    max_az = max_dir_visual_az(i)/pi*180;
                    max_el = max_dir_visual_el(i)/pi*180;
                    HTI = HTI_visual(i);
                else
                    pref_az = pref_dir_vestib_az(i)/pi*180;
                    pref_el = pref_dir_vestib_el(i)/pi*180;
                    max_az = max_dir_vestib_az(i)/pi*180;
                    max_el = max_dir_vestib_el(i)/pi*180;
                    HTI = HTI_vestib(i);
                end
                sel = pref_az < -90; pref_az(sel) = pref_az(sel) + 360;
                sel = max_az < -90; max_az(sel) = max_az(sel) + 360;
                for j=1:2
                    ga1_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+1,:,i);
                    ga2_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+2,:,i);
                    ga3_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+3,:,i);
                    ga4_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+4,:,i);
                    ga5_out(1:num_unique_points) = hidden_outputs(m,(j-1)*5+5,:,i);
                    basis_model_contour_plots2(unique_azimuth, unique_elevation, ...
                        unique_point_azimuth, unique_point_elevation, ...
                        ga1_out, ga2_out, ga3_out, ga4_out, ga5_out, (m-1)*num_stim_types + j, ...
                        sprintf('cell %d - %s - %s - pref %.2f %.2f - max %.2f %.2f - HTI %.3f', ...
                        i, str{j}, curvefit_stim_string{m}, pref_az, pref_el, max_az, max_el, HTI) );
                end
            end
        end
        pause;
    end % plot_hidden_layer
end

sel = pref_dir_visual_az < 0;
pref_dir_visual_az(sel) = pref_dir_visual_az(sel) + 2*pi;
sel = pref_dir_vestib_az < 0;
pref_dir_vestib_az(sel) = pref_dir_vestib_az(sel) + 2*pi;
sel = pref_dir_combin_az < 0;
pref_dir_combin_az(sel) = pref_dir_combin_az(sel) + 2*pi;

figure(500);
xtick = [0:45:360]; ytick = [-90:45:90];
subplot(2,2,1);
scatter(pref_dir_vestib_az/pi*180, pref_dir_vestib_el/pi*180);
xlabel('pref dir az'); ylabel('pref dir el'); 
title(sprintf('vestib - mse %.2e msenorm - %.2e',network_mse,network_norm_mse));
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);
subplot(2,2,2);
scatter(pref_dir_visual_az/pi*180, pref_dir_visual_el/pi*180);
xlabel('pref dir az'); ylabel('pref dir el');
title(sprintf('visual - msereg %.2e mseregnorm - %.2e',network_msereg,network_norm_msereg));
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);
subplot(2,2,3);
scatter(max_dir_vestib_az/pi*180, max_dir_vestib_el/pi*180);
xlabel('max dir az'); ylabel('max dir el'); title('vestib');
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);
subplot(2,2,4);
scatter(max_dir_visual_az/pi*180, max_dir_visual_el/pi*180);
xlabel('max dir az'); ylabel('max dir el'); title('visual');
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);

figure(501);
e = 5e-16;  % some dot products come out minutely larger than 1
ebins = [0:15:180]; cbins = [7.5:15:172.5]; nbins = 12;

[x y z] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
dot_pref = x.*xx + y.*yy + z.*zz;
%dot_pref = ndot_visual_vestib;
sel = dot_pref > 1; dot_pref(sel) = dot_pref(sel)-e; sel = dot_pref < -1; dot_pref(sel) = dot_pref(sel)+e;
ang_pref = acos(dot_pref)/pi*180;
subplot(3,3,1);
n = histc(ang_pref, ebins); bar(cbins, n(1:nbins));
%hist(dot_pref);
ylabel('# neurons'); xlabel('angle (deg)'); title('pref vestib vs pref visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180]);

[x y z] = sph2cart(max_dir_vestib_az, max_dir_vestib_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(max_dir_visual_az, max_dir_visual_el, ones(1,num_hidden_units));
dot_max = x.*xx + y.*yy + z.*zz;
sel = dot_max > 1; dot_max(sel) = dot_max(sel)-e; sel = dot_max < -1; dot_max(sel) = dot_max(sel)+e;
ang_max = acos(dot_max)/pi*180;
subplot(3,3,2);
n = histc(ang_max, ebins); bar(cbins, n(1:nbins));
%hist(dot_max);
ylabel('# neurons'); xlabel('angle (deg)'); title('max vestib vs max visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180]);

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
dot_prefc = x.*xx + y.*yy + z.*zz;
%dot_prefc = ndot_visual_combin;
sel = dot_prefc > 1; dot_prefc(sel) = dot_prefc(sel)-e; sel = dot_prefc < -1; dot_prefc(sel) = dot_prefc(sel)+e;
ang_prefc = acos(dot_prefc)/pi*180;
subplot(3,3,4);
n = histc(ang_prefc, ebins); bar(cbins, n(1:nbins));
%hist(dot_prefc);
ylabel('# neurons'); xlabel('angle (deg)'); title('pref combin vs pref visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180]);

[x y z] = sph2cart(max_dir_combin_az, max_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(max_dir_visual_az, max_dir_visual_el, ones(1,num_hidden_units));
dot_maxc = x.*xx + y.*yy + z.*zz;
sel = dot_maxc > 1; dot_maxc(sel) = dot_maxc(sel)-e; sel = dot_maxc < -1; dot_maxc(sel) = dot_maxc(sel)+e;
ang_maxc = acos(dot_maxc)/pi*180;
subplot(3,3,5);
n = histc(ang_maxc, ebins); bar(cbins, n(1:nbins));
%hist(dot_maxc);
ylabel('# neurons'); xlabel('angle (deg)'); title('max combin vs max visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180]);

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
dot_prefcv = x.*xx + y.*yy + z.*zz;
%dot_prefcv = ndot_vestib_combin;
sel = dot_prefcv > 1; dot_prefcv(sel) = dot_prefcv(sel)-e; sel = dot_prefcv < -1; dot_prefcv(sel) = dot_prefcv(sel)+e;
ang_prefcv = acos(dot_prefcv)/pi*180;
subplot(3,3,7);
n = histc(ang_prefcv, ebins); bar(cbins, n(1:nbins));
%hist(dot_prefc);
ylabel('# neurons'); xlabel('angle (deg)'); title('pref combin vs pref vestib');
set(gca,'xtick',[0:30:180],'xlim',[0 180]);

[x y z] = sph2cart(max_dir_combin_az, max_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(max_dir_vestib_az, max_dir_vestib_el, ones(1,num_hidden_units));
dot_maxcv = x.*xx + y.*yy + z.*zz;
sel = dot_maxcv > 1; dot_maxcv(sel) = dot_maxcv(sel)-e; sel = dot_maxcv < -1; dot_maxcv(sel) = dot_maxcv(sel)+e;
ang_maxcv = acos(dot_maxcv)/pi*180;
subplot(3,3,8);
n = histc(ang_maxcv, ebins); bar(cbins, n(1:nbins));
%hist(dot_maxc);
ylabel('# neurons'); xlabel('angle (deg)'); title('max combin vs max vestib');
set(gca,'xtick',[0:30:180],'xlim',[0 180]);

subplot(3,3,6);
scatter(ang_prefc,ang_prefcv);
cnt = sum(ang_prefcv > ang_prefc); ang_prefcv_gt_ang_prefc = cnt;
cnt0 = sum(ang_prefcv == ang_prefc);
x = [min([min(ang_prefcv) min(ang_prefc)]):0.1:max([max(ang_prefcv) max(ang_prefc)])];
hold on; plot(x,x,'r'); hold off;
ylabel('vestib-combin angle'); xlabel('visual-combin angle'); title(sprintf('pref %d(%d) vestib > visual',cnt,cnt0));
set(gca,'xtick',[0:30:180],'ytick',[0:30:180],'xlim',[0 180],'ylim',[0 180]);
subplot(3,3,9);
scatter(ang_maxc,ang_maxcv);
cnt = sum(ang_maxcv > ang_maxc);
cnt0 = sum(ang_maxcv == ang_maxc);
x = [min([min(ang_maxcv) min(ang_maxc)]):0.1:max([max(ang_maxcv) max(ang_maxc)])];
hold on; plot(x,x,'r'); hold off;
ylabel('vestib-combin angle'); xlabel('visual-combin angle'); title(sprintf('max %d(%d) vestib > visual',cnt,cnt0));
set(gca,'xtick',[0:30:180],'ytick',[0:30:180],'xlim',[0 180],'ylim',[0 180]);

subplot(3,3,3);
% scatter3(ang_prefc,ang_prefcv,ang_pref,'r');
% [b bint] = regress(ang_pref', [ang_prefc' ang_prefcv']);
% hold on; surf([0 180; 0 0], [0 0;180 180], [0 b(1)*180; b(2)*180 b(1)*180]); hold off;
% xlabel('visual-combin'); ylabel('vestib-combin'); zlabel('visual-vestib');
% set(gca,'xtick',[0:30:180],'ytick',[0:30:180],'ztick',[0:30:180],...
%     'xlim',[0 180],'ylim',[0 180],'zlim',[0 180]);
scatter(ang_prefc./ang_pref, ang_prefcv./ang_pref);
xlabel('vis-com / vis-ves'); ylabel('ves-com / vis-ves');

figure(502);
nbins = floor(num_hidden_units/10);
subplot(2,3,1);
hist(HTI_vestib,nbins);
ylabel('# neurons'); xlabel('HTI'); title('vestib');
subplot(2,3,2);
hist(HTI_visual,nbins);
ylabel('# neurons'); xlabel('HTI'); title('visual');
subplot(2,3,3);
hist(HTI_combin,nbins);
ylabel('# neurons'); xlabel('HTI'); title('combined');

x = [0:1:1];
subplot(2,3,4);
scatter(HTI_vestib,HTI_combin);
hold on; plot(x,x,'r'); hold off;
cnt = sum(HTI_combin > HTI_vestib); cnt0 = sum(HTI_combin == HTI_vestib);
ylabel('HTI combined'); xlabel('HTI vestib'); title(sprintf('%d(%d) combin > vestib',cnt,cnt0));
subplot(2,3,5);
scatter(HTI_visual,HTI_combin);
hold on; plot(x,x,'r'); hold off;
cnt = sum(HTI_combin > HTI_visual); cnt0 = sum(HTI_combin == HTI_visual);
ylabel('HTI combined'); xlabel('HTI visual'); title(sprintf('%d(%d) combin > visual',cnt,cnt0));
subplot(2,3,6);
scatter(HTI_vestib,HTI_visual);
hold on; plot(x,x,'r'); hold off;
cnt = sum(HTI_visual > HTI_vestib); cnt0 = sum(HTI_visual == HTI_vestib);
ylabel('HTI visual'); xlabel('HTI vestib'); title(sprintf('%d(%d) visual > vestib',cnt,cnt0));

figure(503);
a = regr_b(1,1:num_hidden_units); sel = a > 1000; a(sel) = 0;
subplot(1,3,1);
hist(a,nbins);
ylabel('# neurons'); xlabel('vestib gain');
subplot(1,3,2);
scatter(HTI_vestib./HTI_visual,a);
ylabel('vestib gain'); xlabel('HTI vestib / HTI visual');
subplot(1,3,3);
scatter(regr_b(1,1:num_hidden_units),regr_b(2,1:num_hidden_units));

figure(504);
subplot(3,3,1);
hist(pref_vestib_mag,nbins);
ylabel('# neurons'); xlabel('mag pref'); title('vestib');
subplot(3,3,2);
hist(pref_visual_mag,nbins);
ylabel('# neurons'); xlabel('mag pref'); title('visual');
subplot(3,3,3);
hist(pref_combin_mag,nbins);
ylabel('# neurons'); xlabel('mag pref'); title('combin');
subplot(3,3,4);
hist(max_vestib_mag,nbins);
ylabel('# neurons'); xlabel('mag max'); title('vestib');
subplot(3,3,5);
hist(max_visual_mag,nbins);
ylabel('# neurons'); xlabel('mag max'); title('visual');
subplot(3,3,6);
hist(max_combin_mag,nbins);
ylabel('# neurons'); xlabel('mag max'); title('combin');
subplot(3,2,5);
scatter(pref_vestib_mag, pref_visual_mag);
x = [min([min(pref_vestib_mag) min(pref_visual_mag)]):0.1:max([max(pref_vestib_mag) max(pref_visual_mag)])];
cnt = sum(pref_visual_mag > pref_vestib_mag); cnt0 = sum(pref_visual_mag == pref_vestib_mag);
pref_visual_mag_gt_pref_vestib_mag = cnt;
hold on; plot(x,x,'r'); hold off;
ylabel('visual mag'); xlabel('vestib mag'); title(sprintf('pref - %d(%d) visual > vestib', cnt,cnt0));
subplot(3,2,6);
scatter(max_vestib_mag, max_visual_mag);
x = [min([min(max_vestib_mag) min(max_visual_mag)]):0.1:max([max(max_vestib_mag) max(max_visual_mag)])];
cnt = sum(max_visual_mag > max_vestib_mag); cnt0 = sum(max_visual_mag == max_vestib_mag);
hold on; plot(x,x,'r'); hold off;
ylabel('visual mag'); xlabel('vestib mag'); title(sprintf('max - %d(%d) visual > vestib', cnt,cnt0));

figure(505);
e = 5e-16;  % some dot products come out minutely larger than 1
nbins = 12; gaze_shift = 45;

[x y z] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_n45_dir_visual_az, pref_n45_dir_visual_el, ones(1,num_hidden_units));
tmp = x.*xx + y.*yy + z.*zz;
sel = tmp > 1; tmp(sel) = tmp(sel)-e; sel = tmp < -1; tmp(sel) = tmp(sel)+e;
ang_gaze_visual = acos(tmp)/pi*180;
rf_shift_ratio_horz_visual = ang_gaze_visual / gaze_shift;

[x y z] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_n45_dir_vestib_az, pref_n45_dir_vestib_el, ones(1,num_hidden_units));
tmp = x.*xx + y.*yy + z.*zz;
sel = tmp > 1; tmp(sel) = tmp(sel)-e; sel = tmp < -1; tmp(sel) = tmp(sel)+e;
ang_gaze_vestib = acos(tmp)/pi*180;
rf_shift_ratio_horz_vestib = ang_gaze_vestib / gaze_shift;

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_n45_dir_combin_az, pref_n45_dir_combin_el, ones(1,num_hidden_units));
tmp = x.*xx + y.*yy + z.*zz;
sel = tmp > 1; tmp(sel) = tmp(sel)-e; sel = tmp < -1; tmp(sel) = tmp(sel)+e;
ang_gaze_combin = acos(tmp)/pi*180;
rf_shift_ratio_horz_combin = ang_gaze_combin / gaze_shift;

subplot(1,3,1);
hist(rf_shift_ratio_horz_visual,nbins);
subplot(1,3,2);
hist(rf_shift_ratio_horz_vestib,nbins);
subplot(1,3,3);
hist(rf_shift_ratio_horz_combin,nbins);


%%%%% start paper plots

figure(601);
e = 5e-16;  % some dot products come out minutely larger than 1
ebins = [0:15:180]; cbins = [7.5:15:172.5]; nbins = 12;

[x y z] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
dot_pref = x.*xx + y.*yy + z.*zz;
%dot_pref = ndot_visual_vestib;
sel = dot_pref > 1; dot_pref(sel) = dot_pref(sel)-e; sel = dot_pref < -1; dot_pref(sel) = dot_pref(sel)+e;
ang_pref = acos(dot_pref)/pi*180;
subplot(2,3,1);
n = histc(ang_pref, ebins); bar(cbins, n(1:nbins));
%hist(dot_pref);
ylabel('# neurons'); xlabel('angle (deg)'); title('angle pref vestib to pref visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180],'ylim',[0 65]);

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
dot_prefc = x.*xx + y.*yy + z.*zz;
%dot_prefc = ndot_visual_combin;
sel = dot_prefc > 1; dot_prefc(sel) = dot_prefc(sel)-e; sel = dot_prefc < -1; dot_prefc(sel) = dot_prefc(sel)+e;
ang_prefc = acos(dot_prefc)/pi*180;
subplot(2,3,2);
n = histc(ang_prefc, ebins); bar(cbins, n(1:nbins));
%hist(dot_prefc);
ylabel('# neurons'); xlabel('angle (deg)'); title('angle pref combin to pref visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180],'ylim',[0 65]);

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
dot_prefcv = x.*xx + y.*yy + z.*zz;
%dot_prefcv = ndot_vestib_combin;
sel = dot_prefcv > 1; dot_prefcv(sel) = dot_prefcv(sel)-e; sel = dot_prefcv < -1; dot_prefcv(sel) = dot_prefcv(sel)+e;
ang_prefcv = acos(dot_prefcv)/pi*180;
subplot(2,3,3);
n = histc(ang_prefcv, ebins); bar(cbins, n(1:nbins));
%hist(dot_prefc);
ylabel('# neurons'); xlabel('angle (deg)'); title('angle pref combin to pref vestib');
set(gca,'xtick',[0:30:180],'xlim',[0 180],'ylim',[0 65]);

subplot(2,2,3);
scatter(ang_prefc,ang_prefcv);
cnt = sum(ang_prefcv > ang_prefc); ang_prefcv_gt_ang_prefc = cnt;
cnt0 = sum(ang_prefcv == ang_prefc);
x = [min([min(ang_prefcv) min(ang_prefc)]):0.1:max([max(ang_prefcv) max(ang_prefc)])];
hold on; plot(x,x,'r'); hold off;
ylabel('vestib-combin angle'); xlabel('visual-combin angle'); title(sprintf('%d vestib-combin > visual-combin',cnt));
set(gca,'xtick',[0:30:180],'ytick',[0:30:180],'xlim',[0 180],'ylim',[0 180]);

subplot(2,2,4);
% scatter3(ang_prefc,ang_prefcv,ang_pref,'r');
% [b bint] = regress(ang_pref', [ang_prefc' ang_prefcv']);
% hold on; surf([0 180; 0 0], [0 0;180 180], [0 b(1)*180; b(2)*180 b(1)*180]); hold off;
% xlabel('visual-combin'); ylabel('vestib-combin'); zlabel('visual-vestib');
% set(gca,'xtick',[0:30:180],'ytick',[0:30:180],'ztick',[0:30:180],...
%     'xlim',[0 180],'ylim',[0 180],'zlim',[0 180]);
scatter(ang_prefc./ang_pref, ang_prefcv./ang_pref);
xlabel('visual-combin / visual-vestib'); ylabel('vestib-combin / visual-vestid');
title('visual-combin + vestib-combin = visual-vestib');

figure(6011);
e = 5e-16;  % some dot products come out minutely larger than 1
ebins = [0:15:180]; cbins = [7.5:15:172.5]; nbins = 12;

[x y z] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
dot_pref = x.*xx + y.*yy + z.*zz;
%dot_pref = ndot_visual_vestib;
sel = dot_pref > 1; dot_pref(sel) = dot_pref(sel)-e; sel = dot_pref < -1; dot_pref(sel) = dot_pref(sel)+e;
ang_pref = acos(dot_pref)/pi*180;
subplot(2,1,1);
n = histc(ang_pref, ebins); bar(cbins, n(1:nbins));
%hist(dot_pref);
ylabel('# neurons'); xlabel('angle (deg)'); title('angle pref vestib to pref visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180],'ylim',[0 25]);

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_visual_az, pref_dir_visual_el, ones(1,num_hidden_units));
dot_prefc = x.*xx + y.*yy + z.*zz;
%dot_prefc = ndot_visual_combin;
sel = dot_prefc > 1; dot_prefc(sel) = dot_prefc(sel)-e; sel = dot_prefc < -1; dot_prefc(sel) = dot_prefc(sel)+e;
ang_prefc = acos(dot_prefc)/pi*180;
subplot(2,2,3);
n = histc(ang_prefc, ebins); bar(cbins, n(1:nbins));
%hist(dot_prefc);
ylabel('# neurons'); xlabel('angle (deg)'); title('angle pref combin to pref visual');
set(gca,'xtick',[0:30:180],'xlim',[0 180],'ylim',[0 65]);

[x y z] = sph2cart(pref_dir_combin_az, pref_dir_combin_el, ones(1,num_hidden_units));
[xx yy zz] = sph2cart(pref_dir_vestib_az, pref_dir_vestib_el, ones(1,num_hidden_units));
dot_prefcv = x.*xx + y.*yy + z.*zz;
%dot_prefcv = ndot_vestib_combin;
sel = dot_prefcv > 1; dot_prefcv(sel) = dot_prefcv(sel)-e; sel = dot_prefcv < -1; dot_prefcv(sel) = dot_prefcv(sel)+e;
ang_prefcv = acos(dot_prefcv)/pi*180;
subplot(2,2,4);
n = histc(ang_prefcv, ebins); bar(cbins, n(1:nbins));
%hist(dot_prefc);
ylabel('# neurons'); xlabel('angle (deg)'); title('angle pref combin to pref vestib');
set(gca,'xtick',[0:30:180],'xlim',[0 180],'ylim',[0 65]);

figure(604);
subplot(2,2,1);
hist(pref_vestib_mag,nbins);
ylabel('# neurons'); xlabel('||pref||'); title('vestib');
subplot(2,2,2);
hist(pref_visual_mag,nbins);
ylabel('# neurons'); xlabel('||pref||'); title('visual');
subplot(2,2,4);
hist(pref_combin_mag,nbins);
ylabel('# neurons'); xlabel('||pref||'); title('combin');
subplot(2,2,3);
scatter(pref_vestib_mag, pref_visual_mag);
x = [min([min(pref_vestib_mag) min(pref_visual_mag)]):0.1:max([max(pref_vestib_mag) max(pref_visual_mag)])];
cnt = sum(pref_visual_mag > pref_vestib_mag); cnt0 = sum(pref_visual_mag == pref_vestib_mag);
pref_visual_mag_gt_pref_vestib_mag = cnt;
hold on; plot(x,x,'r'); hold off;
ylabel('||visual pref||'); xlabel('||vestib pref||'); title(sprintf('%d ||visual pref|| > ||vestib pref||', cnt));

figure(600);
xtick = [0:45:360]; ytick = [-90:45:90];
subplot(2,2,1);
scatter(pref_dir_vestib_az/pi*180, pref_dir_vestib_el/pi*180);
xlabel('pref dir az'); ylabel('pref dir el'); 
%title(sprintf('vestib - mse %.2e msenorm - %.2e',network_mse,network_norm_mse));
title('vestib');
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);
subplot(2,2,2);
scatter(pref_dir_visual_az/pi*180, pref_dir_visual_el/pi*180);
xlabel('pref dir az'); ylabel('pref dir el');
%title(sprintf('visual - msereg %.2e mseregnorm - %.2e',network_msereg,network_norm_msereg));
title('visual');
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);
subplot(2,2,3);
scatter(max_dir_vestib_az/pi*180, max_dir_vestib_el/pi*180);
xlabel('max dir az'); ylabel('max dir el'); title('vestib');
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);
subplot(2,2,4);
scatter(max_dir_visual_az/pi*180, max_dir_visual_el/pi*180);
xlabel('max dir az'); ylabel('max dir el'); title('visual');
set(gca,'xtick',xtick,'ytick',ytick); axis([0 360 -90 90]);

figure(602);
nbins = floor(num_hidden_units/10);
subplot(2,3,1);
hist(HTI_vestib,nbins);
ylabel('# neurons'); xlabel('HTI'); title('vestib');
subplot(2,3,2);
hist(HTI_visual,nbins);
ylabel('# neurons'); xlabel('HTI'); title('visual');
subplot(2,3,3);
hist(HTI_combin,nbins);
ylabel('# neurons'); xlabel('HTI'); title('combined');

x = [0:1:1];
subplot(2,3,4);
scatter(HTI_vestib,HTI_combin);
hold on; plot(x,x,'r'); hold off;
cnt = sum(HTI_combin > HTI_vestib); cnt0 = sum(HTI_combin == HTI_vestib);
ylabel('HTI combined'); xlabel('HTI vestib'); title(sprintf('%d combin > vestib',cnt));
subplot(2,3,5);
scatter(HTI_visual,HTI_combin);
hold on; plot(x,x,'r'); hold off;
cnt = sum(HTI_combin > HTI_visual); cnt0 = sum(HTI_combin == HTI_visual);
ylabel('HTI combined'); xlabel('HTI visual'); title(sprintf('%d combin > visual',cnt));
subplot(2,3,6);
scatter(HTI_vestib,HTI_visual);
hold on; plot(x,x,'r'); hold off;
cnt = sum(HTI_visual > HTI_vestib); cnt0 = sum(HTI_visual == HTI_vestib);
ylabel('HTI visual'); xlabel('HTI vestib'); title(sprintf('%d visual > vestib',cnt));
