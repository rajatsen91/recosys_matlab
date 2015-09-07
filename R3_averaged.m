%Author : Rajat Sen
%for running recosys repeatedly

addpath ./inexact_alm_mc/
addpath ./inexact_alm_mc/pro/
addpath ./CF-toolkit/

%load synthetic_data.mat ;
load netflix_reduced.mat ; 
%load synthetic_data.mat
%load amazon.mat
%load amazon_corr
%load mlens_new
rng(7);

A = Acomp ; %synthetic dataset
B = Bcomp; 

[x , y ] =size(A) ;

eta = eta ; %value of eta(different for different datasets)

%ads = randsample(1:y , round(log(y))); %log(y) number of ads selected

%S_u = randsample(1:x , 100) ; %users selected

biasp = 0:0.05:0.5 ; %different bias probabilities

num_user = 30:10:140 ; %different user numbers

num_ads = 50:10:80 ; 

%num_m = [200 , 300 , 500 , 700 , 1000 , 2000 , 3000 , 4000 , 5000 , 6000 , 7000 , 8000 , 9000];
num_m = [1000 , 1200 , 1400 , 1700 , 2000 , 2200 , 2400 , 2500 , 3000];

pexp = 0:0.025:0.3 ; 
f_all = 50:50:450;

%psucc = zeros(1,length(biasp)) ; %success probability averaged
%psucc1 = zeros(1,length(psucc));
%psucc2 = zeros(1,length(psucc));
%pnsucc = zeros(1,length(psucc));

%naive_sum = zeros(1, length(psucc));

f_lo = 250;

thresh = 1 ; 

%% setting parameters(refer to recosys file)
bias = 1;

method = 0 ;

explorep = 0.1 ;

S_m = 1:y ;

exp_mode = 1 ;

ad_method = 0; 

learning = 0;

naive_thresh1 = 6.5;
naive_thresh2 = 6.6 ; 
% Suu = zeros(length(num_user),170);
% 
% for i = 1:length(num_user)
%     if(i == 1)
%         temp = randsample(1:x , num_user(i));
%         Suu(i,1:num_user(i)) = temp ; 
%     else
%         temp = Suu(i-1 , 1:num_user(i-1)) ;
%         while(1)
%             temp2 = randsample(1:x , num_user(i) - num_user(i-1));
%             tempf = union(temp2 , temp);
%             if(length(tempf) == num_user(i))
%                 break;
%             end
%         end
%         Suu(i,1:num_user(i)) = tempf ; 
%     end
% end




num_sim = 100 ; % number of simulations in averaging one point


%parpool(2);
for j = 7:length(psucc)
    
    warning('off' , 'MATLAB:dispatcher:pathWarning') ;
    
    avgs0 = 0;
    avgs1 = 0;
    avgs2 = 0;
    avg2 = 0;
    avg3 = 0 ; 
    %ads = randsample(1:y , num_ads(j));
%     AA = A(:,randsample(1:y , num_m(j)));
%     [~,q] = size(AA);
%     ads = randsample(1:q , round(log(q)));
%     f_l =(f_lo/y)*q ;
    for i = 1:num_sim
        
        
        [a1_sum , s0 , s1 , s2 ,  T , ad_given , A1 , Rec , naive_avg , naive_success , nss] = recosys2(A ,B , ads, S_u, S_m , biasp(j) , eta , bias , method , 0.1 , exp_mode , ad_method , f_lo , thresh , naive_thresh1 ,naive_thresh2 ,  learning) ;
        
        avgs0 = avgs0 + s0 ;
        avgs1 = avgs1 + s1 ;
        avgs2 = avgs2 + s2 ;
        avg2 = avg2 + nss ; 
        avg3 = avg3 + naive_success;
        if(mod(i,10) == 0)
            i
            j
        end
        
    end
    
    avg2 = avg2/num_sim ;
    avg3 = avg3/num_sim;
    psucc(j) = avgs0/num_sim ;
    psucc1(j) = avgs1/num_sim ;
    psucc2(j) = avgs2/num_sim ;
    naive_sum(j) = avg2 ;
    pnsucc(j) = avg3 ; 
    
end





