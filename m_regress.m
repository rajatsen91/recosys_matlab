function [ out ] = m_regress(V , rt , r , mt , y)
%to be used in the regression step of linear bandits

out = 0;

for i = 1:length(rt)
    out = out + (r) * abs(rt(i) - V(mt(i) , :)*y) ;
end


out = out + norm(y,2) ; 




end

