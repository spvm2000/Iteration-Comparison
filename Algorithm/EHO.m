function Data_o = EHO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ProblemFunction=Data_i.fobj;        DisplayFlag = true;     RandSeed = round(sum(100*clock));
OPTIONS.popsize = Data_i.pop;       OPTIONS.Maxgen = Data_i.maxIter;    OPTIONS.numVar =  Data_i.dim;
OPTIONS.numClan = 5;               OPTIONS.fobj=Data_i.fobj;
OPTIONS.numElephantInEachClan = OPTIONS.popsize/OPTIONS.numClan;        OPTIONS.OrderDependent = false;
OPTIONS.MaxFEs = 0.3E4;     nEvaluations = OPTIONS.popsize;
MinParValue=Data_i.lb;      MaxParValue=Data_i.ub;      
for popindex = 1 : OPTIONS.popsize
    Population(popindex).chrom = X(popindex,:);
    Population(popindex).cost =Cost(popindex);
end
Population = ClearDups(Population, MaxParValue, MinParValue);
for i=1:Data_i.pop
    Population(popindex).cost=Data_i.fobj(Population(popindex).chrom);
end
Population = PopSort(Population);
AverageCost = ComputeAveCost(Population);
MinCost = [Population(1).cost];
AvgCost = [AverageCost];

Keep = 2;
numElephantInEachClan = OPTIONS.numElephantInEachClan*ones(1, OPTIONS.numClan);
dim = OPTIONS.numVar;
alpha = 0.5;
beta = 0.1;
GenIndex = 1;
while GenIndex<= Data_i.maxIter
    for j = 1 : Keep
        chromKeep(j,:) = Population(j).chrom;
        costKeep(j) = Population(j).cost;
    end

    j =1;
    popindex = 1;
    while popindex <= OPTIONS.popsize
        for cindex = 1 : OPTIONS.numClan
            Clan{cindex}(j) = Population(popindex);
            popindex = popindex + 1;
        end 
        j = j+1;
    end
    j = 1;
    popindex =1;
    while popindex <= OPTIONS.popsize
        for cindex = 1 : OPTIONS.numClan
            ClanCenter =  CaculateClanCenter(OPTIONS, Clan, cindex);
            NewClan{cindex}(j).chrom = Clan{cindex}(j).chrom ...
                + alpha*(Clan{cindex}(1).chrom - Clan{cindex}(j).chrom).* rand(1,dim);
            if sum( NewClan{cindex}(j).chrom - Clan{cindex}(j).chrom) == 0
                NewClan{cindex}(j).chrom = beta * ClanCenter ;
            end
            popindex = popindex + 1;
        end % end for cindex
        j = j+1;
    end
    for cindex = 1 : OPTIONS.numClan
        NewClan{cindex}(end).chrom = MinParValue...
            + (MaxParValue - MinParValue + 1) .* rand(1,OPTIONS.numVar);
    end
    SavePopSize = OPTIONS.popsize;
    for i=1:OPTIONS.numClan
        OPTIONS.popsize = numElephantInEachClan(i);
        % Make sure the population does not have duplicates.
        NewClan{i} = ClearDups(NewClan{i}, MaxParValue, MinParValue);
        % Make sure each individual is legal.        
        NewClan{i} = FeasibleFunction(OPTIONS, NewClan{i},Data_i);
        % Calculate cost
        NewClan{i} = CostFunction(OPTIONS, NewClan{i});
        % the number of fitness evaluations
        nEvaluations = nEvaluations +  OPTIONS.popsize;
        % Sort from best to worst
        NewClan{i} = PopSort(NewClan{i});
    end
    OPTIONS.popsize = SavePopSize;
    Population = CombineClan(OPTIONS, NewClan);
    Population = PopSort(Population);
    n = length(Population);
    for k3 = 1 : Keep
        Population(n-k3+1).chrom = chromKeep(k3,:);
        Population(n-k3+1).cost = costKeep(k3);
    end
    [Population,index] = PopSort(Population);
    [AverageCost, nLegal] = ComputeAveCost(Population);
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    IterCurve(GenIndex)=Population(1).cost;
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Population(1).cost
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Population(1).cost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
    end    
    GenIndex = GenIndex+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=GenIndex;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [Population] = ClearDups(Population, MaxParValue, MinParValue)
    if length(MaxParValue) == 1
        for i = 1 : length(Population)
            Chrom1 = sort(Population(i).chrom);
            for j = i+1 : length(Population)
                Chrom2 = sort(Population(j).chrom);
                if isequal(Chrom1, Chrom2)
                    parnum = ceil(length(Population(j).chrom) * rand);
                    Population(j).chrom(parnum) = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand);
                end
            end
        end
    else
        for i = 1 : length(Population)
            Chrom1 = sort(Population(i).chrom);
            for j = i+1 : length(Population)
                Chrom2 = sort(Population(j).chrom);
                if isequal(Chrom1, Chrom2)
                    parnum = ceil(length(Population(j).chrom) * rand);
                    Population(j).chrom(parnum) = floor(MinParValue(parnum) ...
                        + (MaxParValue(parnum) - MinParValue(parnum) + 1) * rand);
                end
            end
        end
    end
end

function [Population, indices] = PopSort(Population)
    popsize = length(Population);
    Cost = zeros(1, popsize);
    indices = zeros(1, popsize);
    for i = 1 : popsize
        Cost(i) = Population(i).cost;
    end
    [Cost, indices] = sort(Cost, 2, 'ascend');
    Chroms = zeros(popsize, length(Population(1).chrom));
    for i = 1 : popsize
        Chroms(i, :) = Population(indices(i)).chrom;
    end
    for i = 1 : popsize
        Population(i).chrom = Chroms(i, :);
        Population(i).cost = Cost(i);
    end
end

function [AveCost, nLegal] = ComputeAveCost(Population)
    Cost = [];
    nLegal = 0;
    for i = 1 : length(Population)
        if Population(i).cost < inf
            Cost = [Cost Population(i).cost];
            nLegal = nLegal + 1;
        end
    end
    % Compute average cost.
    AveCost = mean(Cost);
end

function Population = CombineClan(OPTIONS, NewClan)
    j =1;
    popindex = 1;
    while popindex <= OPTIONS.popsize
        for clanindex = 1 : OPTIONS.numClan
            Population(popindex)   = NewClan{clanindex}(j);
            popindex = popindex + 1;
        end
        j = j+1;
    end
end

function [Population] = FeasibleFunction(OPTIONS, Population,Data_i)
    for i = 1 : OPTIONS.popsize
        for k = 1 : OPTIONS.numVar
            Population(i).chrom(k) = max(Population(i).chrom(k), Data_i.lb(k));
            Population(i).chrom(k) = min(Population(i).chrom(k), Data_i.ub(k));
        end
    end
end

function [Population] = CostFunction(OPTIONS, Population)
    for popindex = 1 : OPTIONS.popsize
        Population(popindex).cost = OPTIONS.fobj(Population(popindex).chrom);    
    end
end

function ClanCenter =  CaculateClanCenter( OPTIONS, Clan, cindex)

ClanCenter = zeros(1, OPTIONS.numVar);

for Elephantindex = 1 : OPTIONS.numElephantInEachClan
    ClanCenter = ClanCenter + Clan{cindex}(Elephantindex).chrom;
end

ClanCenter = (1/OPTIONS.numElephantInEachClan)* ClanCenter;
end