function Data_o = BBO(Data_i,Data_o)
rand_num=[];                           
ti=clock;                                
pop.Position = Data_i.X;                        
pop.Cost=Data_i.F_value;                   
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
CostFunction=Data_i.fobj;        % Cost Function

nVar=Data_i.dim;             % Number of Decision Variables


VarMin=Data_i.lb(1);         % Decision Variables Lower Bound
VarMax=Data_i.ub(1);         % Decision Variables Upper Bound

%% BBO Parameters

MaxIt=Data_i.maxIter;          % Maximum Number of Iterations

nPop=Data_i.pop;               % Number of Habitats (Population Size)

KeepRate=0.2;                   % Keep Rate
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats

nNew=nPop-nKeep;                % Number of New Habitats

% Migration Rates
mu=linspace(1,0,nPop);          % Emmigration Rates
lambda=1-mu;                    % Immigration Rates

alpha=0.9;

pMutation=0.1;

sigma=0.02*(VarMax-VarMin);

%% Initialization

% Empty Habitat
habitat.Position=[];
habitat.Cost=[];

% Create Habitats Array
pop=repmat(habitat,nPop,1);

% Initialize Habitats
for i=1:nPop
    pop(i).Position=Data_i.X(i,:);
    pop(i).Cost=Data_i.F_value(i);
end

% Sort Population
[~, SortOrder]=sort([pop.Cost]);
pop=pop(SortOrder);

% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop
for it=1:MaxIt
    newpop=pop;
    for i=1:nPop
        for k=1:nVar
            % Migration
            if rand<=lambda(i)
                % Emmigration Probabilities
                EP=mu;
                EP(i)=0;
                EP=EP/sum(EP);
                
                % Select Source Habitat
                j=RouletteWheelSelection(EP);
                
                % Migration
                newpop(i).Position(k)=pop(i).Position(k) ...
                    +alpha*(pop(j).Position(k)-pop(i).Position(k));
            end
            
            % Mutation
            if rand<=pMutation
                newpop(i).Position(k)=newpop(i).Position(k)+sigma*randn;
            end
        end
        
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
        % Evaluation
        newpop(i).Cost=CostFunction(newpop(i).Position);
    end
    
    % Sort New Population
    [~, SortOrder]=sort([newpop.Cost]);
    newpop=newpop(SortOrder);
    
    % Select Next Iteration Population
    pop=[pop(1:nKeep)
         newpop(1:nNew)];
     
    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Update Best Solution Ever Found
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>pop(1).Cost
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=pop(1).Cost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=SortOrder(1);
    end    
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    IterCurve(it)=BestSol.Cost;
    
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function j=RouletteWheelSelection(P)
    r=rand;
    C=cumsum(P);
    j=find(r<=C,1,'first');
end