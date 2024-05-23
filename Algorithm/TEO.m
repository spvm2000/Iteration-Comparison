function Data_o = TEO(Data_i,Data_o)
rand_num=[];                        
ti=clock;                               
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
best_value=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);    
index=Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter);
CostFunction=Data_i.fobj;       nVar=Data_i.dim;        VarSize=[1 nVar];   VarMin=Data_i.lb(1);
VarMax=Data_i.ub(1);            MaxIt=Data_i.maxIter;       nPop=Data_i.pop;    
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=1;                                                    %全局最小适应度值
c1=1.1;               % Personal Learning Coefficient
c2=c1*2;
pCrossover=0.7;
nCrossover=2*round(pCrossover*nPop/2);
pMutation=0.4;
nMutation=round(pMutation*nPop);
mu=0.02;
sigma=0.1*(VarMax-VarMin); 
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    pop(i).Position=X(i,:);
    pop(i).Cost=Cost(i);
end
[pop, F]=NonDominatedSorting(pop);
pop=CalcCrowdingDistance(pop,F);
[pop, F]=SortPopulation(pop);
for it=1:Data_i.maxIter
    popc=repmat(empty_individual,nCrossover/2,2);
    for k=1:nCrossover/2
        i1=randi([1 nPop]);
        p1=pop(i1);
        
        i2=randi([1 nPop]);
        p2=pop(i2);
        
        MaxRank=max([pop.Rank]);
        c=c1+c2*(MaxIt-it)/MaxIt;
        ratio=it/MaxIt;
        
        [popc(k,1), popc(k,2)]=Crossover(p1,p2,MaxRank,c,ratio);
        
        popc(k,1).Position=max(min(popc(k,1).Position,VarMax),VarMin);
        popc(k,2).Position=max(min(popc(k,2).Position,VarMax),VarMin);

        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        i=randi([1 nPop]);
        p=pop(i);
        popm(k).Position=Mutate(p.Position,mu,sigma);
        popm(k).Position=max(min(popm(k).Position,VarMax),VarMin);
        popm(k).Cost=CostFunction(popm(k).Position);
    end
    pop=[pop
         popc
         popm];
    [pop, F]=NonDominatedSorting(pop);
    pop=CalcCrowdingDistance(pop,F);
    pop=SortPopulation(pop);
    pop=pop(1:nPop);
    [pop, F]=NonDominatedSorting(pop);
    pop=CalcCrowdingDistance(pop,F);
    [pop, F]=SortPopulation(pop);
    F1=pop(F{1});
    IterCurve(it)=F1.Cost;
    for i=1:Data_i.pop              %求全局最优
        if pop(i).Cost<best_value
            best_value=pop(i).Cost;
            index=i;
        end    
    end    
    if best_value >= IterCurve(it)
        best_value = IterCurve(it);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=best_value;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    end    
end    
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [pop, F]=NonDominatedSorting(pop)
    nPop=numel(pop);
    for i=1:nPop
        pop(i).DominationSet=[];
        pop(i).DominatedCount=0;
    end
    F{1}=[];
    for i=1:nPop
        for j=i+1:nPop
            p=pop(i);
            q=pop(j);
            if Dominates(p,q)
                p.DominationSet=[p.DominationSet j];
                q.DominatedCount=q.DominatedCount+1;
            end
            if Dominates(q.Cost,p.Cost)
                q.DominationSet=[q.DominationSet i];
                p.DominatedCount=p.DominatedCount+1;
            end
            pop(i)=p;
            pop(j)=q;
        end
        if pop(i).DominatedCount==0
            F{1}=[F{1} i];
            pop(i).Rank=1;
        end
    end
    k=1;
    while true
        Q=[];
        for i=F{k}
            p=pop(i);
            for j=p.DominationSet
                q=pop(j);
                q.DominatedCount=q.DominatedCount-1;
                if q.DominatedCount==0
                    Q=[Q j]; %#ok
                    q.Rank=k+1;
                end   
                pop(j)=q;
            end
        end
        if isempty(Q)
            break;
        end
        F{k+1}=Q; %#ok
        k=k+1; 
    end
end

function pop=CalcCrowdingDistance(pop,F)
    nF=numel(F);
    for k=1:nF 
        Costs=[pop(F{k}).Cost];
        nObj=size(Costs,1);
        n=numel(F{k});
        d=zeros(n,nObj);
        for j=1:nObj
            [cj, so]=sort(Costs(j,:));
            d(so(1),j)=inf;
            for i=2:n-1
                d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
            end
            d(so(end),j)=inf;
        end
        for i=1:n
            pop(F{k}(i)).CrowdingDistance=sum(d(i,:));
        end 
    end
end

function [pop, F]=SortPopulation(pop)
    % Sort Based on Crowding Distance
    [~, CDSO]=sort([pop.CrowdingDistance],'descend');
    pop=pop(CDSO);
    % Sort Based on Rank
    [~, RSO]=sort([pop.Rank]);
    pop=pop(RSO);
    % Update Fronts
    Ranks=[pop.Rank];
    MaxRank=max(Ranks);
    F=cell(MaxRank,1);
    for r=1:MaxRank
        F{r}=find(Ranks==r);
    end
end

function [P1, P2]=Crossover(P1,P2,MaxRank,c,ratio)
    x1=P1.Position;
    r1=P1.Rank;
    x2=P2.Position;
    r2=P2.Rank;
    nVar=numel(x1);
    Tenv2=(ones(1,nVar)-c*rand(1,nVar)).*x2;
    y1=Tenv2+(x1-Tenv2)*exp(-1*ratio*(r1/MaxRank));
    Tenv1=(ones(1,nVar)-c*rand(1,nVar)).*x1;
    y2=Tenv1+(x2-Tenv1)*exp(-1*ratio*(r2/MaxRank));
    P1.Position=y1;
    P2.Position=y2;
end

function b=Dominates(x,y)
    if isstruct(x)
        x=x.Cost;
    end
    if isstruct(y)
        y=y.Cost;
    end
    b=all(x<=y) && any(x<y);
end

function y=Mutate(x,mu,sigma)
    nVar=numel(x);
    nMu=ceil(mu*nVar);
    j=randsample(nVar,nMu);
    if numel(sigma)>1
        sigma = sigma(j);
    end
    y=x;
    y(j)=x(j)+sigma.*randn(size(j));
end