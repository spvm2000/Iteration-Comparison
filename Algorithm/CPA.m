function Data_o = CPA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
d=Data_i.dim;       Lb=Data_i.lb;       Ub=Data_i.ub;       Dim=Data_i.dim;

group_iter = 2;
attraction_rate = 0.8;
growth_rate = 2;
reproduction_rate = 1.8;
nPop=Data_i.pop;
nCP = round(nPop*(1/3),0);
nPrey = nPop-round(nPop*(1/3),0);

ini.position=[];
ini.cost=[];
life=repmat(ini,nPop,1);
NewCP=repmat(ini,nCP*group_iter+nCP,1);

for i=1:Data_i.pop
    life(i).position=X(i,:);
    life(i).cost=Cost(i);
end
costs=Cost;
[fmin,I]=min(costs);
best=life(I).position;
N_iter=1;

while N_iter<=Data_i.maxIter
    [CP,life]=CPA_Grouping(life,nCP,nPrey,group_iter,N_iter);
    [NewCP,a]=CPA_Growth(CP,NewCP,Dim,attraction_rate,growth_rate,group_iter);
    NewCP=CPA_Reproduction(CP,NewCP,Dim,reproduction_rate,a);
    NewCP=CPA_UpdateCost(NewCP,Dim,Lb,Ub,Data_i);
    [life,BestIndex]=CPA_Combine(NewCP,life);
    BestSol=life(BestIndex);
    fmin=BestSol.cost;
    best=BestSol.position;
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fmin
       Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fmin;
       Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=BestIndex;
    end
    IterCurve(N_iter)=fmin;
    N_iter=N_iter+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=N_iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end


function [CP,pop]=CPA_Grouping(pop,nCP,nPrey,Group_Iter,it)

empty_Prey.position=[];
empty_Prey.cost=[];
Select=reshape(1:nPrey,nCP,[]);

% First, sort the population
if it==1
    for x=1:nCP+nPrey
        costs(x,:)=pop(x).cost;
    end
else
    for x=1:(2+Group_Iter)*nCP+nPrey
        costs(x,:)=pop(x).cost;
    end
end

for i=1:nCP+nPrey
    [~, I]=min(costs);
    index(i)=I;
    costs(I)=10^30;
end
pop=pop(index);

Plant=pop(1:nCP);
Prey=pop(nCP+1:end);

empty_CP.Plant=[];
empty_CP.Prey=repmat(empty_Prey,0,1);
empty_CP.nPrey=0;
CP=repmat(empty_CP,nCP,1);

% Group the Carnivorous Plants with their nearby Preys
for q=1:nCP
    CP(q).Plant=Plant(q);
    for r=1:nPrey/nCP
        CP(q).Prey=[CP(q).Prey
            Prey(Select(q,r))];
        CP(q).nPrey=CP(q).nPrey+1;
    end
end
end

function [NewCP,a]=CPA_Growth(CP,NewCP,Dim,HuntingChance,alpha,Group_Iter)

a=1;
nCP=numel(CP);

for Grp=1:nCP
    for Group_cycle=1:Group_Iter
        v=randi(CP(Grp).nPrey);
        if HuntingChance>rand   %Growth of Carnivorous Plant by hunting prey
 
            Step=alpha.*rand(1,Dim);
            NewCP(a).position=Step.*CP(Grp).Plant.position...
                +(1-Step).*CP(Grp).Prey(v).position;
            
            a=a+1;
        else                    %Mating of Prey
            %To ensure prey v and prey j are different
            u=v;

            while v==u
                u=randi(CP(Grp).nPrey);
            end
            
            Step=alpha.*rand(1,Dim);
            if CP(Grp).Prey(v).cost<CP(Grp).Prey(u).cost
                Step=1-Step;   %So that it moves to good prey
            end
            
            NewCP(a).position=Step.*CP(Grp).Prey(u).position...
                +(1-Step).*CP(Grp).Prey(v).position;
            
            a=a+1;
        end
    end
end
end

function NewCP=CPA_Reproduction(CP,NewCP,Dim,alpha,a)
%Mating of Carnivorous Plant
for i=1:numel(CP)
    for j=1:Dim
        
        %To ensure plant j and plant v are different
        v=i;
        while v==i
            v=randi(numel(CP));
        end
        
        Step=CP(i).Plant.position(j)-CP(v).Plant.position(j);
        if CP(v).Plant.cost<CP(i).Plant.cost
            Step=-Step;   %So that it moves to good plant
        end
        
        NewCP(a).position(j)=CP(1).Plant.position(j)+alpha.*rand.*(Step);
    end
    a=a+1;
end
end

function NewCP=CPA_UpdateCost(NewCP,Dim,Lb,Ub,Data_i)
n=numel(NewCP);
for i=1:n
    Flag4ub=NewCP(i).position>Ub;
    Flag4lb=NewCP(i).position<Lb;
    NewCP(i).position=(NewCP(i).position.*(~(Flag4ub+Flag4lb)))...
        +(rand(1,Dim).*(Ub-Lb)+Lb).*(Flag4ub+Flag4lb);
    NewCP(i).cost=Data_i.fobj(NewCP(i).position);
end
end

function [pop,BestIndex]=CPA_Combine(NewCP,pop)
pop=[pop
     NewCP];
for i=1:size(pop,1)
    costs(i,:)=pop(i).cost;
end
[~,BestIndex]=min(costs);
end