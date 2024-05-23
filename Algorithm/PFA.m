function Data_o = PFA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
pop = Data_i.X;                        
pop_fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ObjectiveFunction=Data_i.fobj;      pop_size = Data_i.pop;      problem_size = Data_i.dim;
N = pop_size;       dim = problem_size;         max_iter = Data_i.maxIter;
Lb=Data_i.lb;       Ub=Data_i.ub;      
fitness=pop_fitness;
[fit_global, index]=min(pop_fitness);
pop_old=pop;

global_pop=pop(index,:); path_p=global_pop;
path_old=path_p;
ub=Ub; lb=Lb;

IterCurve(1)=fit_global;
meanc(1)=mean(pop_fitness);

iter=0;
while iter < max_iter
    iter = iter+1;
    u1=-1+(1-(-1)).*rand(N,1);
    u2=-1+(1-(-1)).*rand(1,dim);
    r3=rand(1,dim);
    eps=(1 - iter/(max_iter)).*u1;
    A=u2.*exp(-2.*iter/max_iter);
    path_p(1,:) = path_p(1,:) + ...
                ((2).*r3(1,:)).*(path_old(1,:) - path_p(1,:)) + 1.*A(1,:);                                                    % check upper bound for pathfinder
    path_p(1,:)=BoundaryCheck(path_p(1,:),Data_i.ub,Data_i.lb,Data_i.dim);
    path_old=path_p;
    path_fitness=ObjectiveFunction(path_p(1,:));
    [fitness, index] = sort(fitness);                                                 
    pop=pop(index,:);
    if path_fitness<fitness(1)                                                        % compare new fitness with old
        fitness(1)=path_fitness;
        pop(1,:)=path_p(1,:);    
    else
        fitness(1)=fitness(1);
        pop(1,:)=pop(1,:);
    end
    for ii=2:pop_size
%            d=sqrt(abs(((pop(ii,:)).^2 - (pop(ii-1,:)).^2)));
       d=sqrt(sum(((pop(ii,:)) - (pop(ii-1,:))).^2));                        % calculate D(i,j)
       alpha=1+rand(1,dim);                                                            % assign alpha
       beta=1+rand(1,dim);                                                             % and beta
       pop(ii,:) = pop(ii,:) ...
            + ((alpha(1,:)).*rand).*(global_pop(1,:) - pop(ii,:)) ...
            + ((beta(1,:)).*rand).*(pop_old(ii-1,:) - pop(ii,:)) + 1.*(eps(ii,1).*d);         % update position of followers

       pop(ii,:)=BoundaryCheck(pop(ii,:),Data_i.ub,Data_i.lb,Data_i.dim);
    end
    for i=1:pop_size
        fitness(i)=ObjectiveFunction(pop(i,:));                               	% calculate new fitness of followers
        if fitness(i) < pop_fitness(i)                                          % update position via new fitness
            pop_fitness(i) = fitness(i);
            pop_old(i,:) = pop(i, :);
            bsf_index = i;
            
        else
            pop_fitness(i) = pop_fitness(i);
            pop_old(i,:) = pop_old(i, :);
        end
    end
    [fit_global, index] = min(pop_fitness);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fit_global
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fit_global;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    end
    global_pop = pop_old(index,:);
    path_p=global_pop;
    IterCurve(iter)=fit_global;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iter;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end