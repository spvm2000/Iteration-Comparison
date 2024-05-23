function X = initialization(pop,ub,lb,dim)
    X = zeros(pop,dim); 
    for i = 1:pop
       for j = 1:dim
           X(i,j) = (ub(j) - lb(j))*rand() + lb(j);  
       end
    end
end
