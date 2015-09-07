function [ I ] = random_rating_sample(ent , num)



I = [];

entries = exp(ent);


l = length(entries);

dist = zeros(1,l+1);

sum = 0;
for i = 1:l
    sum = sum + entries(i);
end

x = entries/sum;

sum = 0;
dist(1) = 0;
for i = 1:l
    sum = sum + x(i);
    dist(i+1) = sum;
end


for i = 1:num
    bit = rand(1);
    for j = 1:l
        if(bit >= dist(j) && bit < dist(j+1))
            I = [I , j];
            break;
        end
    end

end



    

