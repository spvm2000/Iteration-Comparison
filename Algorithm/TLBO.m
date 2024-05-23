function Data_o = TLBO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve = zeros(1,Data_i.maxIter);
CostFunction=Data_i.fobj;       nVar = Data_i.dim;      VarSize = [1 nVar];   VarMin=Data_i.lb(1);
VarMax=Data_i.ub(1);        MaxIt=Data_i.maxIter;   nPop=Data_i.pop;
empty_individual.Position = [];
empty_individual.Cost = [];
pop = repmat(empty_individual, nPop, 1);
BestSol.Cost = inf;
for i=1:nPop
    pop(i).Position = X(i,:);
    pop(i).Cost = Cost(i);
    
    if pop(i).Cost < BestSol.Cost
        BestSol = pop(i);
    end
end

for it=1:MaxIt
    Mean = 0;
    for i=1:nPop
        Mean = Mean + pop(i).Position;
    end
    Mean = Mean/nPop;
    Teacher = pop(1);
    for i=2:nPop
        if pop(i).Cost < Teacher.Cost
            Teacher = pop(i);
        end
    end
    for i=1:nPop
        % Create Empty Solution
        newsol = empty_individual;
        
        % Teaching Factor
        TF = randi([1 2]);
        
        % Teaching (moving towards teacher)
        newsol.Position = pop(i).Position ...
            + rand(VarSize).*(Teacher.Position - TF*Mean);
        
        % Clipping
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        
        % Evaluation
        newsol.Cost = CostFunction(newsol.Position);
        
        % Comparision
        if newsol.Cost<pop(i).Cost
            pop(i) = newsol;
            if pop(i).Cost < BestSol.Cost
                BestSol = pop(i);
            end
        end
    end
    for i=1:nPop
        
        A = 1:nPop;
        A(i)=[];
        j = A(randi(nPop-1));
        
        Step = pop(i).Position - pop(j).Position;
        if pop(j).Cost < pop(i).Cost
            Step = -Step;
        end
        
        % Create Empty Solution
        newsol = empty_individual;
        
        % Teaching (moving towards teacher)
        newsol.Position = pop(i).Position + rand(VarSize).*Step;
        
        % Clipping
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        
        % Evaluation
        newsol.Cost = CostFunction(newsol.Position);
        
        % Comparision
        if newsol.Cost<pop(i).Cost
            pop(i) = newsol;
            if pop(i).Cost < BestSol.Cost
                BestSol = pop(i);
                if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>BestSol.Cost
                    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestSol.Cost;
                    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                end
            end
        end
    end
    IterCurve(it) = BestSol.Cost;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end