% % do bootstrap to compute the range of threshold or bias
% % -GY 09/21/07
% 
% aa1 = dlmread('data1.txt'); % aa(:,1) is x, aa(:,2) is y, aa(:,3) is repetition
% %aa2 = dlmread('data2.txt'); % the data that need to be compared with later
% aa2=[];
% 
% bootnum = 1000;
% for i = 1 : bootnum
%     i
%     for j = 1 : length(aa1(:,1)) % bootstrap for each heading
%         right_choice1 = round( aa1(j,2)*aa1(j,3) ); 
%         choice1(1:right_choice1) = 1; % right heading labelled by 1
%         choice1(right_choice1+1:aa1(j,3) ) = 2; % left heading labelled by 2
%         
%         for pp = 1 : aa1(j,3)
%             choice_perm1 = choice1( randperm(aa1(j,3)) );
%             choice_heading1(pp) = choice_perm1(1);
%         end
%         temp = find(choice_heading1==1);
%         aa1_perm(j,1) = aa1(j,1);
%         aa1_perm(j,2) = length(temp) / aa1(j,3); 
%         aa1_perm(j,3) = aa1(j,3);        
%     end
%     
%     [bb,tt] = cum_gaussfit_max1(aa1_perm);
%     bias1(i) =  bb;
%     threshold1(i) = tt;
%     
%     % now for the second one
%     if length(aa2)>1  % if aa2 exist
%         for j = 1 : length(aa2(:,1)) % bootstrap for each heading
%             right_choice2 = round( aa2(j,2)*aa2(j,3) ); 
%             choice2(1:right_choice2) = 1; % right heading labelled by 1
%             choice2(right_choice2+1:aa2(j,3) ) = 2; % left heading labelled by 2
% 
%             for pp = 1 : aa2(j,3)
%                 choice_perm2 = choice2( randperm(aa2(j,3)) );
%                 choice_heading2(pp) = choice_perm2(1);
%             end
%             temp = find(choice_heading2==1);
%             aa2_perm(j,1) = aa2(j,1);
%             aa2_perm(j,2) = length(temp) / aa2(j,3); 
%             aa2_perm(j,3) = aa2(j,3);        
%         end
% 
%         [bb,tt] = cum_gaussfit_max1(aa2_perm);
%         bias2(i) =  bb;
%         threshold2(i) = tt;
%     end
% end

% % now compare the confidence interval
aa = dlmread('aa.txt');
if mean(aa(:,1)) > mean(aa(:,2))
    p = length( find(aa(:,1)<mean(aa(:,2))) );    
    p = p/1000
else
    p = length( find(aa(:,1)>mean(aa(:,2))) );
    p = p/1000
end
