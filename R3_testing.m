% addpath ./inexact_alm_mc/
% addpath ./inexact_alm_mc/pro/
% 
% % % count = 0;
% % % A = R ; 
% % % [x , y] = size(R) ; 
% % % num_above = zeros(x,1);
% % % for i = 1:x
% % %     maxi = max(A(i,:));
% % %     mini = min(A(i,:));
% % %     A(i,:) = 10*(A(i,:) - mini)/(maxi - mini) ; 
% % % end
% % 
% % 
% % 
% % for i = 1:x 
% %     count = 0;
% %     for j = 1:y
% %         if(A(i,j) >= 5.4)
% %             count = count + 1;
% %         end
% %     end
% %     num_above(i) = count ; 
% % end
% % 
% % 
% % % A = Amazonratingreduced ; 
% % % 
% % % index = [];
% % % 
% % % [x ,y] = size(A);
% % % 
% % % for i = 1:y
% % %     if(sum(A(:,i)) ~= 0)
% % %         index = [index , i] ; 
% % %     end
% % % end
% % % 
% % % 
% % % B = A(:,index) ; 
% % 
% %%
%

% %%
% load mlens_reduced ; 
% 
% [x , y] = size(D);
% 
% [U , S , V] = OptSpace(full(D),7 , 20 , 1e-4);
% Acomp = U*S*V' ;  
%%

[x, y] = size(Bcomp);
Acomp = zeros(x,y);
for i = 1:x
    row = Bcomp(i,:);
    mm = max(row);
    mi = min(row);
    Acomp(i,:) = 10*(row - mi)/(mm - mi) ;
end

%%

% %%
% % 
% % eta = 0; 
% % fl = round(10*log(y));
% % 
% % for i = 1:x
% %     eta = eta + kthvalue(Acomp(i,:), y - fl) ; 
% % end
% % 
% % eta = eta/x ; 
% % 
% % 
% % %%
% % [x , y] = size(A);
% % eta = 6 ; 
% % count = zeros(1,x);
% % 
% % for i = 1:x
% %     count(i) = length(find(A(i,:) >= eta)) ; 
% % end
% 

% 
% [x,y] = size(M);
% [row, col] = find(M);
% 
% countr = zeros(1,x);
% countc = zeros(1,y);
% 
% for i = 1:x
%     countr(i) = length(find(M(i,:)));
% end
% 
% 
% 
% for i = 1:y
%     countc(i) = length(find(M(:,i)));
% end

% 
% [x, y] = size(M);
% 
% cc = zeros(1,y);
% 
% for i = 1:y
%     cc(i) = length(find(M(:,i)));
% end


%%
% count = zeros(1,x);
% eta = 8.6 ; 
% for i = 1:x
%     count(i) = length(find(Acomp(i,:) > eta)) ; 
% end


 %%
% nm = 100:100:y ; 
% Mall = zeros(x,y, length(nm));
% 
% for i = 1:length(nm)
%     Ma(: , 1:nm(i) , i) = Acomp(: , randsample(1:y , nm(i)));
% end

%%
% load netflix_reduced.mat ; 
eta = 0;
[x , y] = size(Acomp);

for i = 1:x
    eta = eta + kthvalue(Acomp(i,:), y - 150);
end

eta = eta/x ; 


%eta = 8 ; 
count = zeros(1,x);

for i = 1:x
    count(i) = length(find(Acomp(i,:) > eta));
end





    

