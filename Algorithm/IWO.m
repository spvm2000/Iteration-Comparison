function Data_o = IWO(Data_i,Data_o)                         
ti=clock;                               
pop.Position = Data_i.X;                        
pop.Cost=Data_i.F_value;                  
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);

CostFunction=Data_i.fobj;   nVar = Data_i.dim;  VarSize = [1 nVar];
VarMin =Data_i.lb(1);   VarMax =Data_i.ub(1);    
empty_plant.Position = [];
empty_plant.Cost = [];
%% IWO Parameters
MaxIt = Data_i.maxIter;            
nPop0 = round(Data_i.pop*0.5);             
nPop = Data_i.pop;              
Smin = 0;               
Smax = 5;               
Exponent = 2;           
sigma_initial = 0.5;   
sigma_final = 0.001;	

pop = repmat(empty_plant, nPop0, 1);
for i = 1:nPop0
    pop(i).Position = Data_i.X(i,:);
    pop(i).Cost = Data_i.F_value(i);
end

IterCurve=zeros(1,Data_i.maxIter);

for it = 1:Data_i.maxIter
    sigma = ((MaxIt - it)/(MaxIt - 1))^Exponent * (sigma_initial - sigma_final) + sigma_final;
    Costs = [pop.Cost];
    BestCost = min(Costs);
    WorstCost = max(Costs);
    newpop = [];
    for i = 1:numel(pop)
        ratio = (pop(i).Cost - WorstCost)/(BestCost - WorstCost);
        S = floor(Smin + (Smax - Smin)*ratio);
        for j = 1:S
            newsol = empty_plant;
            newsol.Position = pop(i).Position + sigma * randn(VarSize);
            newsol.Position = max(newsol.Position, VarMin);
            newsol.Position = min(newsol.Position, VarMax);
            newsol.Cost = CostFunction(newsol.Position);
            newpop = [newpop
                      newsol];
        end    
    end
    pop = [pop
           newpop];
    [~, SortOrder]=sort([pop.Cost]);
    pop = pop(SortOrder);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>pop(SortOrder(1)).Cost
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=pop(SortOrder(1)).Cost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=SortOrder(1);
    end    
    if numel(pop)>nPop
        pop = pop(1:nPop);
    end
    BestSol = pop(1);
    IterCurve(it) = BestSol.Cost;           
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end