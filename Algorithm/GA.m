function Data_o = GA(Data_i,Data_o)
t=clock;                             
cnt=1;                               
eps = 1e-2;                          
pc = 0.9;                            
pm = 0.01;                           
global sn                            
X = Data_i.X;
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
for gen=1:Data_i.maxIter
    for b = 1:Data_i.pop
        
        fitness(b)=Data_i.fobj(X(b,:));
    end

    [dad] = selection(X,fitness);

    newpop = crossover(dad,pc,Data_i.fobj);
    newpop = mutation(newpop,pm,Data_i.lb,Data_i.ub);
    for c = 1:Data_i.pop
        fitness_temp(c) = Data_i.fobj(newpop(c,:));
        if fitness_temp(c)<fitness(c)
            fitness(c)=fitness_temp(c);
            X(c,:)=newpop(c,:);
            if fitness(c)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness(c);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=c;
            end    
        end
    end
    [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)] = min(fitness_temp);              
    IterCurve(gen)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter) = etime(clock,t);               
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=gen;
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end