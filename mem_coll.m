function [M_comp]= mem_coll(D , S_u ,type)
if(type == 0)
    D = round(D);
    M_comp = full(D);
    [x,~] = size(D);
    for i = 1:length(S_u)
        [~,meanPred]=predictPreferencePD(full(D(S_u(i),:)), D(setdiff(1:x , S_u(i)),:),10,1);
        M_comp(S_u(i),:) = meanPred;
    end
end

if(type == 1)
    [x,~] = size(D);
    M_comp = full(D);
    user_col = zeros(length(S_u),x);
    
    
    for j = 1:length(S_u)
        m_u = find(full(D(S_u(j),:)));
        for i = 1:x
            
            m_i = find(full(D(i,:)));
            inter = intersect(m_u, m_i);
            
            R_u = full(D(S_u(j),inter));
            R_i = full(D(i,inter));
            uni = union(m_i,m_u) ;
            
            
            
            if(isempty(inter))
                user_col(j,i) = 0;
            else
                user_col(j,i) = (R_u*R_i'/(norm(R_u,2)*norm(R_i,2)))*(length(inter)/length(uni));
            end
            
        end
    end
    for i = 1:length(S_u)
        nz = find(D(S_u(i),:));
        avg = mean(full(D(S_u(i),nz)));
        prefpred = predictPreferenceMemBased(user_col(i,:) , D , avg);
        M_comp(S_u(i),:) = prefpred ;
    end
    
end