function [Rec] = linear_recos(D, Q_m , S_u , A)
%Author : Rajat Sen
%Linear Bandit recommendation scheme
tic ;
n_u = length(S_u);

Rec = zeros(n_u , Q_m);

[U , S , V] = OptSpace(full(D),30 , 10 , []);

V = S*V' ;

V = V' ;
% MC = inexact_alm_mc(D , 1e-4 , 30);
% 
% U = MC.U;
% V = MC.V ; 

[x , r] = size(U);

[~ , l] = size(A);

for u = 1:n_u
    for t = 1:Q_m
        mt = find(full(D(S_u(u),:)));
        rt = full(D(S_u(u) , mt));
        
        %         cvx_begin quiet
        %         variable y(r,1)
        %         minimize (m_regress(V , rt , r , mt , y))
        %         subject to
        %         cvx_end
        
        
        rt = rt' ; 
        M = V(mt,:);
        
        y = M\rt ; 
        
        beta = sqrt(2/3)*(1/(2*t))^(0.25);
        
        
        pseudo = y'*normr(V)' ;
        [pmax , imax] = max(pseudo);
        
        ps2 = (V(imax,:)/norm(V(imax,:)',2))*normr(V)' ;
        
        movies = find(abs(ps2) > (1 - (beta^2)/2));
        
        if(~isempty(movies))
            Rec(u,t) = randsample(movies,1);
        else
            Rec(u,t) = imax ;
        end
        D(S_u(u) , Rec(u,t)) = A(S_u(u) , Rec(u,t));
        
        
        
        
    end
    
end


toc;




end

