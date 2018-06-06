% to test if the cell have temporal tuning
% if it shows a significant monophasic or biphasic response to at least 2
% nearby directions, we say this neuron is responsive.
% LBY 20171123

function resp_T = respon_True(sigT)

sigT_trans = nan*ones(5,8);
sigT_trans2 = sigT_trans;
sigT_trans3 = sigT_trans;
sigT_trans(2:4,:) = reshape(sigT(2:25),[3 8]);
sigT_trans(1,:) = sigT(1);
sigT_trans(5,:) = sigT(26);

sigT_trans2 = sigT_trans(:,2:end);
sigT_trans3 = sigT_trans(2:end,:);

resp1 = sigT_trans(2:4,1:7)+sigT_trans2(2:4,:);
resp2 = sigT_trans(1:4,:)+sigT_trans3(:,:);

if sum(find(resp1 == 2))>0 || sum(find(resp2 == 2))>0
    resp_T = 1;
else
    resp_T = 0;
end

end