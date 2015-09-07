function [ user_col , movie ,act_rec] = get_collrec(D , u , num , k  )

%collaborative filtering recommendation
%D : sparse matrix
%u : active user
%num : number of movies to be recommended
%k : top k most matched users contribute to average rating

[x,~] = size(D);

user_col = zeros(1,x);

m_u = find(full(D(u,:)));
for i = 1:x
    
    m_i = find(full(D(i,:)));
    inter = intersect(m_u, m_i);
    
    R_u = full(D(u,inter));
    R_i = full(D(i,inter));
    uni = union(m_i,m_u) ;
    
   
    
    if(isempty(inter))
        user_col(i) = 0;
    else
        user_col(i) = (R_u*R_i'/(norm(R_u,2)*norm(R_i,2)))*(length(inter)/length(uni));
    end
    
end



[xsort, isort] = sort(user_col , 'descend');

movie = [];

for i = 1:x
    m_i = find(full(D(isort(i),:)));
    
    movie = union(movie, m_i);
    if(length(movie) > num)
        break;
    end
end

act_rec = zeros(1 , length(movie));

for i = 1:length(movie)
    sum = 0;
    for j = 1:k
        if(full(D(isort(j), movie(i))) > 0)
            sum = sum + xsort(j);
            act_rec(i) = act_rec(i) + xsort(j)*full(D(isort(j), movie(i)));
        end
        
    end
    if(sum > 0)
        act_rec(i) = act_rec(i)/sum ;
    end
end



























end

