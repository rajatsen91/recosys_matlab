function [ a1_sum , success , T ,ad_given , A1 ,Rec] = recosys(A , ads , S_u , S_m , biasp , eta , bias , method , explorep , exp_mode , ad_method)
%Author : Rajat Sen
%Function for running recommendation and Bias detection
%ads : set of ads
%A : full matrix
%S_u : users in question
%S_m : movies in question
%biasp : bias probability
%eta : threshold
%bias : bias of no bias
%method : method of recommendation
%explorep : explore probability
%exp_mode : mode of recomendation (top-k or relaxed)


%a1-sum : counts of dislikes for top movies(our criterion)
%T : threshold changing with time
%success : 1 if correctly detected, 0 otherwise
%ad_given : lets you know the slot of ads given
%A1 : the completed matrix, according to the recommendation system





[x , y] = size(A);

degree = randraw('pareto' , [4,3] , x); % the degree distribution of the entrier rated by users (closer to real world data)

degree = round(degree);

A1 = zeros(x,y); %to be filled in with the rated entries of the users, the rest are 0

f_l = round(10*log(y)) ; % f(m) in our model

Q_m = 50 ; % length of detection process



%% creates A1, the incomplete matrix
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

num1 = x;
%% method = 0 : matrix completion
if(method == 0)
    
    D = sparse(A1);
    [Mcomp  iter  svp] = inexact_alm_mc(D , 1e-4 , 15) ;
    A1_comp = Mcomp.U*(Mcomp.V') ;  %% completed matrix
    
    
end

%% method = 1 : collaborative filtering
if(method == 1)
    A1_comp = A1;
    for i = 1:length(S_u)
        [user_col , mmm , act_rec] = get_collrec(A1 , S_u(i) , Q_m+20 , 60);
        A1_comp(S_u(i) , mmm) = act_rec;
    end
end

%% method = 2 : true ideal recommendation
if (method == 2)
    A1_comp = A ;
    
end

%% some more parameters
num_ads = length(ads) ;



n_u = length(S_u); % number of users with us

Rec = zeros(n_u , Q_m); %matrix consisting of recommendations to users over Q_m slots

ad_given = zeros(n_u , Q_m) ;






%% recosys with bias

if(bias == 1)
    %%recommendation algorithm with bias
    for u = 1:n_u
        [u_sort , I ] = sort(A1_comp(S_u(u),:) , 'descend');
        ad_rat = A1_comp(S_u(u),ads);
        u_ind = zeros(1,y);
        cols = find(A1(S_u(u),:));
        u_ind(cols) = 1;
        ad_ind = zeros(1, num_ads) ;
        [ad_sort , Ia] = sort(ad_rat , 'descend');
        adc = 1;
        mc = 1;
        rc = 1;
        adcount = 0;
        cout = 0;
        
        while(cout <= 500)
            cout = cout + 1;
            bit = rand(1);
            if(rc > Q_m)
                break;
            end
            if(bit >= (1-biasp) && adc <= num_ads) %biased ads
                if(ad_method == 1)
                    while(adc <= num_ads)
                        if(u_ind(ads(Ia(adc))) ~= 1)
                            break;
                        else
                            adc = adc + 1;
                        end
                    end
                    
                    if(adc <= num_ads)   %showing ads
                        Rec(u,rc) = ads(Ia(adc));
                        if(A(S_u(u) , ads(Ia(adc))) < eta)
                            ad_given(u,rc) = 2;
                        else
                            ad_given(u,rc) = 1;
                        end
                        u_ind(ads(Ia(adc))) = 1;
                        rc = rc + 1;
                        adc = adc + 1 ;
                    end
                else
                    cin = 0;
                    while(cin <= 4*Q_m)
                        cin = cin + 1;
                        if(length(find(ad_ind)) >= num_ads)
                            break;
                        end
                        new_ad = randsample(setdiff(1:num_ads,find(ad_ind)) , 1);
                        if(ad_ind(new_ad) == 0)
                            ad_ind(new_ad) = 1;
                            u_ind(ads(new_ad)) = 1;
                            Rec(u,rc) = ads(new_ad);
                            rc = rc + 1;
                            if(A(S_u(u) , ads(new_ad)) < eta)
                                ad_given(u,rc) = 2;
                            else
                                ad_given(u,rc) = 1;
                            end
                            break;
                        end
                        adcount = adcount + 1;
                        if(adcount > num_ads)
                            break;
                        end
                        
                    end
                end
            elseif(bit < (1-biasp) && bit >= explorep) %usual recosys exploit
                if(exp_mode == 1)
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
                        mc = mc + 1;
                        rc = rc + 1;
                    end
                else
                    cin = 0;
                    while(cin <= 4*Q_m)
                        if(sum(u_ind(I(1:f_l+10))) ~= f_l + 10) ;
                            rindex = randsample(setdiff(I(1:f_l+10),find(u_ind)) , 1);
                            if(u_ind(rindex) ~= 1)
                                break ;
                            end
                        end
                    end
                    u_ind(rindex) = 1 ;
                    Rec(u,rc) = rindex ;
                    rc = rc + 1 ;
                end
                
            else
                cin = 0;
                while(cin <= 4*Q_m)  %explore with small probability
                    cin = cin + 1;
                    index = randsample(setdiff(1:y,find(u_ind)) , 1);
                    if(u_ind(index) == 0)
                        Rec(u,rc) = index;
                        u_ind(index) = 1;
                        rc = rc + 1;
                        break;
                    end
                end
                
                
            end
        end
    end
end

%% unbiased algorithm
if(bias == 0)
    
    for u = 1:n_u
        [u_sort , I ] = sort(A1_comp(S_u(u),:) , 'descend');
        u_ind = zeros(1,y);
        cols = find(A1(S_u(u),:));
        u_ind(cols) = 1;
        mc = 1;
        rc = 1;
        cout = 1;
        while(cout <= 500)
            cout = cout + 1 ;
            if(rc > Q_m)
                break;
            end
            bit = rand(1);
            if(bit > (1-explorep))   %%usdual recosys exploit
                if(exp_mode == 1)
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
                        mc = mc + 1;
                        rc = rc + 1;
                    end
                else
                    cin = 0;
                    while(cin <= 4*Q_m)
                        cin = cin+1;
                        rindex = randsample(setdiff(I(1:f_l+10),find(u_ind)) , 1);
                        if(u_ind(rindex) ~= 1)
                            break ;
                        end
                    end
                    u_ind(rindex) = 1 ;
                    Rec(u,rc) = rindex ;
                    rc = rc + 1 ;
                end
                
            else
                cin = 0;
                while(cin <= 5*Q_m)
                    cin = cin + 1;
                    
                    index = randsample(setdiff(1:y,find(u_ind)) , 1); %random explore
                    if(u_ind(index) == 0)
                        Rec(u,rc) = index;
                        u_ind(index) = 1;
                        rc = rc + 1;
                        break;
                    end
                end
                
            end
        end
        
    end
end


%% detection algorithm (naive part not coded)


M_hc = zeros(y,Q_m);

n_u = length(S_u);

for t = 1:Q_m
    count = zeros(1,y);
    for i = 1:n_u
        if(A(S_u(i),Rec(i,t)) < eta)
            count(Rec(i,t)) = count(Rec(i,t)) + 1 ;
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
a1_sum = zeros(1,Q_m);

prob = 0;
for i = 1:Q_m
    a = i ;
    prob = prob + a/(f_l - i + 1) ;
    
    %     a = 2*i - 1 ;
    
    %T(i) = (n_u*i*a/(f_l-i)) + 2*sqrt(n_u*i*a*log(y));
    %T(i) = min(max(2*a*log(y)/log(f_l/(2*n_u*i*a)) , (2*a*n_u*i/f_l)) , n_u*i) ;
    %     T(i) = 2*a*log(y)/log(f_l/(2*n_u*i*a)) ;
    %     Tu(i) = n_u*i ;
    %     Tl(i) = (2*a*n_u*i/f_l);
    beta = (a+0.5)*log(y)/(n_u*prob)- 1 ;
    hatp = exp(1+Lambert_W(beta/exp(1)))*n_u*prob ;
    hatbeta = (a+0.5)*log(y)/(hatp)- 1 ;
    T(i) = exp(1+Lambert_W(beta/exp(1)))*n_u*prob ;
    %T(i) = 3*n_u*prob ;
    %T(i) = exp(1+Lambert_W(hatbeta/exp(1)))*hatp ;
    [xsort, isort] = sort(M_hc(:,i)' , 'descend');
    for j = 1:a
        a1_sum(i) = a1_sum(i) + xsort(j);
    end
end

success = 0;

for i = 1:Q_m
    if (a1_sum(i) >= T(i)) %bias detected
        success = 1 ;
        break ;
    end
    
end


if(bias == 0)
    success = mod(1+success , 2) ; %wrong if bias detected in unbiased case
end





end

