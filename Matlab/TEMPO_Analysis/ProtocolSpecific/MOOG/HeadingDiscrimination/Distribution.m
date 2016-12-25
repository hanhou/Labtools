% run ROC and test significance from chance (0.5)
congruent = dlmread('congruency.txt');
con_size = size(congruent);
opposite = dlmread('congruency_opp.txt');
opp_size = size(opposite);
for i = 1 : 3
    aa(i) = rocN(congruent(:,i)', opposite(:,i)', 100);
end

% permutation test significance
perm_num=1000;
for n = 1: 1000
    for j = 1 : 3
        aa_temp(1,1:con_size(1)) = congruent(:,j)';
        aa_temp(1,(con_size(1)+1) : (con_size(1)+opp_size(1))) = opposite(:,j)';
        aa_perm = aa_temp( randperm( con_size(1)+opp_size(1) ));
        congruent_perm(j, 1:con_size(1) )= aa_perm(1:con_size(1));
        opposite_perm(j, 1: opp_size(1) )= aa_perm( (con_size(1)+1):(con_size(1)+opp_size(1)) );
        aa_permutation(n,j) = rocN(congruent_perm(j,:), opposite_perm(j,:), 100);
    end
end

% significant test
bin = 0.005;
x_bin = 0 : bin : 1;
for k=1:3
    hist_perm(k,:) = hist( aa_permutation(:,k), x_bin );  % for permutation
    bin_sum = 0;
    n = 0;
    while ( n < (aa(k)/bin) )
         n = n+1;
         bin_sum = bin_sum + hist_perm(k, n);
         if aa(i) > 0.5                  % note it's two tail test
            p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
         else
            p{k} = 2* bin_sum / perm_num;
         end
    end
end
aa
p
