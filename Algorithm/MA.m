function Data_o = MA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
ObjectiveFunction=Data_i.fobj;    ProblemSize=[1 Data_i.dim];
LowerBound=Data_i.lb(1);          UpperBound=Data_i.ub(1);
GlobalBest.Cost=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
GlobalBest.Position=X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
%% Mayfly Parameters
MaxIt=Data_i.maxIter;                 % Maximum Number of Iterations
nPop=round(Data_i.pop*0.5); nPopf=Data_i.pop-nPop;          % Population Size (males and females)
g=0.8;                      % Inertia Weight
gdamp=1;                    % Inertia Weight Damping Ratio
a1=1.0;                     % Personal Learning Coefficient
a2=1.5; a3=1.5;             % Global Learning Coefficient
beta=2;                     % Distance sight Coefficient
dance=5;                    % Nuptial Dance
fl=1;                       % Random flight
dance_damp=0.8;             % Damping Ratio
fl_damp=0.99;
nc=20;                      % Number of Offsprings (also Parnets)
nm=round(0.05*nPop);        % Number of Mutants
mu=0.01;
VelMax=0.1*(UpperBound-LowerBound);     VelMin=-VelMax;
empty_mayfly.Position=[];
empty_mayfly.Cost=[];
empty_mayfly.Velocity=[];
empty_mayfly.Best.Position=[];
empty_mayfly.Best.Cost=[];

Mayfly=repmat(empty_mayfly,nPop,1);         % Males
Mayflyf=repmat(empty_mayfly,nPopf,1);       % Females

for i=1:nPop
    Mayfly(i).Position=X(i,:);
    Mayfly(i).Velocity=zeros(ProblemSize);
    Mayfly(i).Cost=fitness(i);

    Mayfly(i).Best.Position=Mayfly(i).Position;
    Mayfly(i).Best.Cost=Mayfly(i).Cost;
    if Mayfly(i).Best.Cost<GlobalBest.Cost
        GlobalBest=Mayfly(i).Best;
    end
end    

for i=1:nPop
    Mayflyf(i).Position=X(nPop+i,:);
    Mayflyf(i).Velocity=zeros(ProblemSize);
    Mayflyf(i).Cost=fitness(nPop+i);
    if Mayflyf(i).Best.Cost<GlobalBest.Cost
        GlobalBest=Mayflyf(i).Best;
    end
end    

IterCurve=zeros(1,MaxIt);
for it=1:MaxIt
    for i=1:nPop
        % Update Females
        e=unifrnd(-1,+1,ProblemSize);
        rmf=(Mayfly(i).Position-Mayflyf(i).Position);
        if Mayflyf(i).Cost>Mayfly(i).Cost
            Mayflyf(i).Velocity = g*Mayflyf(i).Velocity ...
                +a3*exp(-beta.*rmf.^2).*(Mayfly(i).Position-Mayflyf(i).Position);
        else
            Mayflyf(i).Velocity = g*Mayflyf(i).Velocity+fl*(e);
        end
        % Apply Velocity Limits
        Mayflyf(i).Velocity = max(Mayflyf(i).Velocity,VelMin);
        Mayflyf(i).Velocity = min(Mayflyf(i).Velocity,VelMax);
        % Update Position
        Mayflyf(i).Position = Mayflyf(i).Position + Mayflyf(i).Velocity;
        % Velocity Mirror Effect
        %IsOutside=(Mayflyf(i).Position<LowerBound | Mayflyf(i).Position>UpperBound);
        %Mayflyf(i).Velocity(IsOutside)=-Mayflyf(i).Velocity(IsOutside);
        % Position Limits
        Mayflyf(i).Position = max(Mayflyf(i).Position,LowerBound);
        Mayflyf(i).Position = min(Mayflyf(i).Position,UpperBound);
        % Evaluation
        Mayflyf(i).Cost = ObjectiveFunction(Mayflyf(i).Position);

        % Update Global Best (Uncomment if you use the PGB-IMA version)
        %if Mayflyf(i).Best.Cost<GlobalBest.Cost
        %    GlobalBest=Mayflyf(i).Best;
        %end
    end

    for i=1:nPop
        % Update Males
        rpbest=(Mayfly(i).Best.Position-Mayfly(i).Position);
        rgbest=(GlobalBest.Position-Mayfly(i).Position);
        e=unifrnd(-1,+1,ProblemSize);
        % Update Velocity
        if Mayfly(i).Cost>GlobalBest.Cost
            Mayfly(i).Velocity = g*Mayfly(i).Velocity ...
                +a1*exp(-beta.*rpbest.^2).*(Mayfly(i).Best.Position-Mayfly(i).Position) ...
                +a2*exp(-beta.*rgbest.^2).*(GlobalBest.Position-Mayfly(i).Position);
        else
            Mayfly(i).Velocity = g*Mayfly(i).Velocity+dance*(e);
        end
        % Apply Velocity Limits
        Mayfly(i).Velocity = max(Mayfly(i).Velocity,VelMin);
        Mayfly(i).Velocity = min(Mayfly(i).Velocity,VelMax);
        % Update Position
        Mayfly(i).Position = Mayfly(i).Position + Mayfly(i).Velocity;
        % Velocity Mirror Effect
        %IsOutside=(Mayfly(i).Position<LowerBound | Mayfly(i).Position>UpperBound);
        %Mayfly(i).Velocity(IsOutside)=-Mayfly(i).Velocity(IsOutside);
        % Position Limits
        Mayfly(i).Position = max(Mayfly(i).Position,LowerBound);
        Mayfly(i).Position = min(Mayfly(i).Position,UpperBound);
        % Evaluation
        Mayfly(i).Cost = ObjectiveFunction(Mayfly(i).Position);

        % Update Personal Best
        if Mayfly(i).Cost<Mayfly(i).Best.Cost
            Mayfly(i).Best.Position=Mayfly(i).Position;
            Mayfly(i).Best.Cost=Mayfly(i).Cost;
            % Update Global Best
            if Mayfly(i).Best.Cost<GlobalBest.Cost
                GlobalBest=Mayfly(i).Best;
                if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>GlobalBest.Cost
                    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GlobalBest.Cost;
                    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                end
            end
        end
    end
    [~, SortMayflies]=sort([Mayfly.Cost]);
    Mayfly=Mayfly(SortMayflies);
    [~, SortMayflies]=sort([Mayflyf.Cost]);
    Mayflyf=Mayflyf(SortMayflies);

    MayflyOffspring=repmat(empty_mayfly,nc/2,2);
    for k=1:nc/2
        % Select Parents
        i1=k;
        i2=k;
        p1=Mayfly(i1);
        p2=Mayflyf(i2);
        % Apply Crossover
        [MayflyOffspring(k,1).Position, MayflyOffspring(k,2).Position]=Crossover(p1.Position,p2.Position,LowerBound,UpperBound);
        % Evaluate Offsprings
        MayflyOffspring(k,1).Cost=ObjectiveFunction(MayflyOffspring(k,1).Position);
        if MayflyOffspring(k,1).Cost<GlobalBest.Cost
            GlobalBest=MayflyOffspring(k,1);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>GlobalBest.Cost
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GlobalBest.Cost;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=k;
            end
        end

        MayflyOffspring(k,2).Cost=ObjectiveFunction(MayflyOffspring(k,2).Position);
        if MayflyOffspring(k,2).Cost<GlobalBest.Cost
            GlobalBest=MayflyOffspring(k,2);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>GlobalBest.Cost
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GlobalBest.Cost;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=k;
            end
        end

        MayflyOffspring(k,1).Best.Position = MayflyOffspring(k,1).Position;
        MayflyOffspring(k,1).Best.Cost = MayflyOffspring(k,1).Cost;
        MayflyOffspring(k,1).Velocity= zeros(ProblemSize);
        MayflyOffspring(k,2).Best.Position = MayflyOffspring(k,2).Position;
        MayflyOffspring(k,2).Best.Cost = MayflyOffspring(k,2).Cost;
        MayflyOffspring(k,2).Velocity= zeros(ProblemSize);
    end
    MayflyOffspring=MayflyOffspring(:);
    % Mutation
    MutMayflies=repmat(empty_mayfly,nm,1);
    for k=1:nm
        % Select Parent
        i=randi([1 nc]);
        p=MayflyOffspring(i);
        %p=Mayfly(i);
        MutMayflies(k).Position=Mutate(p.Position,mu,LowerBound,UpperBound);
        % Evaluate Mutant
        MutMayflies(k).Cost=ObjectiveFunction(MutMayflies(k).Position);
        if MutMayflies(k).Cost<GlobalBest.Cost
            GlobalBest=MutMayflies(k);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>GlobalBest.Cost
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GlobalBest.Cost;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=k;
            end
        end
        MutMayflies(k).Best.Position = MutMayflies(k).Position;
        MutMayflies(k).Best.Cost = MutMayflies(k).Cost;
        MutMayflies(k).Velocity= zeros(ProblemSize);
    end
    % Create Merged Population
    MayflyOffspring=[MayflyOffspring
        MutMayflies]; %#ok
    split=round((size(MayflyOffspring,1))/2);
    newmayflies=MayflyOffspring(1:split);
    Mayfly=[Mayfly
        newmayflies]; %#ok
    newmayflies=MayflyOffspring(split+1:size(MayflyOffspring,1));
    Mayflyf=[Mayflyf
        newmayflies]; %#ok
    [~, SortMayflies]=sort([Mayfly.Cost]);
    Mayfly=Mayfly(SortMayflies);
    Mayfly=Mayfly(1:nPop); % Keep best males
    [~, SortMayflies]=sort([Mayflyf.Cost]);
    Mayflyf=Mayflyf(SortMayflies);
    Mayflyf=Mayflyf(1:nPopf); % Keep best females
    IterCurve(it)=GlobalBest.Cost;
    g=g*gdamp;
    dance = dance*dance_damp;
    fl = fl*fl_damp;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);            
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [off1, off2]=Crossover(x1,x2,LowerBound,UpperBound)
    L=unifrnd(0,1,size(x1));
    off1=L.*x1+(1-L).*x2;
    off2=L.*x2+(1-L).*x1;
    % Position Limits
    off1=max(off1,LowerBound); off1=min(off1,UpperBound);
    off2=max(off2,LowerBound); off2=min(off2,UpperBound);
end

function y=Mutate(x,mu,LowerBound,UpperBound)
    nVar=numel(x);
    nmu=ceil(mu*nVar);
    j=randsample(nVar,nmu);
    sigma(1:nVar)=0.1*(UpperBound-LowerBound);
    y=x;
    y(j)=x(j)+sigma(j).*(randn(size(j))');
    y=max(y,LowerBound); y=min(y,UpperBound);
end

