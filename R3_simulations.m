%simulations on R3 dataset
%bias vs non-bias

%load R3_comp.mat ;
addpath OptSpace_matlab/ ;
%load mlens_completed.mat ; 

load synthetic_data.mat ; 




c = clock;
rng(c(6));

A = M_final ;    % actual dataset completed


[x , y] = size(A);

degree = randraw('pareto' , [4,3] , x);

degree = round(degree);

A1 = zeros(x,y);

method = 2;



f_l = round(10*log(y)) ; 

%rng('default');

for i = 1:x
    frac = 0.2;
    from_rating = round(frac*degree(i));
    l1 = random_rating_sample(A(i,:) , from_rating);
    l2 = randsample(1:y , degree(i) - from_rating);
    l = [l1, l2] ;
    for j = 1:length(l)
        A1(i,l(j)) = A(i,l(j));
    end
end

num1 = 500;

if(method == 0)
    
    [U1 , S1 , X1] = OptSpace(A1 , [] ,[] , []);
    
    
    A1_comp = U1*S1*X1' ;
    num1 = x; 
end

if(method == 1)
    A1_comp = A1;
    for i = 1:num1
        [recos , act_rec , user_col , mmm] = get_collrec(A1 , i , max(Q_m , 40));
        for j = 1:length(recos)
            A1_comp(i,recos(j)) = act_rec(j);
        end
        if(mod(i,50) == 0)
            i
        end
    end
end

if (method == 2)
    A1_comp = A ; 
    num1 = x ; 
end




% bias = 0;
% 
% explore = 1;
% 
% num_ads = round(log(y)) ;
% 
% Q_m = round(2*log(y)) ; 
% 
% ads = randsample(1:y , num_ads);
% 
% 
% 
% Rec = zeros(x , Q_m);
% 
% ad_given = zeros(x,Q_m);
% 
% T = zeros(1,Q_m) ;
% 
% n_u = 170;

%S_u = randsample(1:num1 , n_u) ;  %set of selected users



eta = zeros(1, length(S_u));

eta_fixed = 5 ;

for u = 1:n_u
    %eta(u) = 1.2*mean(A(u,:));
    eta(u) = eta_fixed ;
    %eta(u) = max(1.1*mean(A(S_u(u),:)) , 0.5*max(A(S_u(u),:)));
end

n_u




%% Biased form of the algorithm



if(bias == 1)
    display('Running Biased algorithm')
    for u = 1:x
        [u_sort , I ] = sort(A1_comp(u,:) , 'descend');
        ad_rat = A1_comp(u,ads);
        %eta(u) = u_sort(Q_m+5);
        u_ind = zeros(1,y);
        cols = find(A1(u,:));
        u_ind(cols) = 1;
        ad_ind = zeros(1, num_ads) ;
        [ad_sort , Ia] = sort(ad_rat , 'descend');
        adc = 1;
        mc = 1;
        rc = 1;
        
        bit = 0;
        
        while(1)
            %                 c = clock;
            %                 rng(c(6));
            %                 rng('shuffle');
            bit = rand(1);
            if(rc > Q_m)
                break;
            end
            if(bit >= 0.5 && adc <= num_ads)
                while(adc <= num_ads)
                    if(u_ind(ads(Ia(adc))) ~= 1)
                        break;
                    else
                        adc = adc + 1;
                    end
                end
                
                if(adc <= num_ads)
                    Rec(u,rc) = ads(Ia(adc));
                    ad_given(u,rc) = 1;
                    u_ind(ads(Ia(adc))) = 1;
                    rc = rc + 1;
                    adc = adc + 1 ; 
                    %A1_comp(u,ads(Ia(adc))) = A(u,ads(Ia(adc)));
                end
            elseif(bit < 0.5 && bit >= 0.1)
                while(mc <= y)
                    if(u_ind(I(mc)) ~= 1)
                        break;
                    else
                        mc = mc + 1;
                    end
                end
                if(mc <= y)
                    Rec(u,rc) = I(mc);
                    u_ind(I(mc)) = 1;
                    %A1_comp(u,I(mc)) = A(u , I(mc));
                    mc = mc + 1;
                    rc = rc + 1;
                end
            else
                while(1)
                    index = randsample(1:y , 1);
                    if(u_ind(index) == 0)
                        Rec(u,rc) = index;
                        u_ind(index) = 1;
                        %A1_comp(u,index) = A(i,index);
                        rc = rc + 1;
                        break;
                    end
                end
                
                
            end
        end
    end
end




if(bias == 0)
    display('Running Unbiased algorithm')
    for u = 1:x
        [u_sort , I ] = sort(A1_comp(u,:) , 'descend');
        %eta(u) = u_sort(Q_m+10);
        u_ind = zeros(1,y);
        cols = find(A1(u,:));
        u_ind(cols) = 1;
        mc = 1;
        rc = 1;
        while(1)
            if(rc > Q_m)
                break;
            end
            bit = rand(1);
            if(bit > 0.1)
                while(mc <= y)
                    if(u_ind(I(mc)) ~= 1)
                        break;
                    else
                        mc = mc + 1;
                    end
                end
                if(mc <= y)
                    Rec(u,rc) = I(mc);
                    u_ind(I(mc)) = 1;
                    %A1_comp(u,I(mc)) = A(u , I(mc));
                    mc = mc + 1;
                    rc = rc + 1;
                end
            else
                while(1)
                    
                    index = randsample(1:y , 1);
                    if(u_ind(index) == 0)
                        Rec(u,rc) = index;
                        u_ind(index) = 1;
                        rc = rc + 1;
                        %display('explore'); 
                        break;
                    end
                end
                
            end
        end
        
    end
end

%% Detection algorithm

%eta = 1.5 times the average rating for an user

display('Running detection algorithm....');


M_hc = zeros(y,Q_m);

n_u = length(S_u);

for t = 1:Q_m
    count = zeros(1,y);
    for i = 1:n_u
        if(A(S_u(i),Rec(S_u(i),t)) < eta(i))
            count(Rec(S_u(i),t)) = count(Rec(S_u(i),t)) + 1 ; 
        end
    end
    
    for i = 1:y
        if(t == 1)
            M_hc(i,t) = count(i);
        else
            M_hc(i,t) = M_hc(i,t-1) + count(i);
        end
    end
end


T = zeros(1,Q_m);
p = zeros(1,Q_m);
a1_sum = zeros(1,Q_m);

for i = 1:Q_m
    a = 2*i - 1 ; 
    T(i) = (2*n_u*i*a/f_l) + 2*sqrt(n_u*i*a*log(y));
    p(i) = T(i)/a;
    [xsort, isort] = sort(M_hc(:,i)' , 'descend');
    for j = 1:a
        a1_sum(i) = a1_sum(i) + xsort(j);
    end
end
success = 0;

for i = 1:Q_m
    if (a1_sum(i) >= T(i)/2)
        success = 1;
        break ; 
    end
end







            
        






















