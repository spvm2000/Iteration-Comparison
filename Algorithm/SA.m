function Data_o = SA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                           
LightRays = Data_i.X;                
fitness=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
nVar=Data_i.dim;
VarMax=Data_i.ub(1);
VarMin=Data_i.lb(1);
T0=0.1;         % Initial Temp.
alpha=0.99;     % Temp. Reduction Rate
nMove=5;        % Number of Neighbors per Individual
mu = 0.5;       % Mutation Rate
sigma = 0.1*(Data_i.ub(1)-Data_i.lb(1));    % Mutation Range (Standard Deviation)
MaxIt=Data_i.maxIter;
MaxSubIt=20;                    % Maximum Number of Sub-iterations
sigma = 0.1*(VarMax-VarMin);    % Mutation Range (Standard Deviation)
nPop=Data_i.pop;                % Population Size
% Create Structure for Individuals
empty_individual.Position=[];
empty_individual.Cost=[];
% Create Population Array
pop=repmat(empty_individual,nPop,1);
% Initialize Best Solution
BestSol.Cost=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
for i=1:nPop
    % Initialize Position
    pop(i).Position=Data_i.X(i,:);     
    % Evaluation
    pop(i).Cost=Data_i.F_value(i);           
    if pop(i).Cost<BestSol.Cost
        BestSol.Cost=pop(i).Cost;                                   
    end
end

T=T0;
for it=1:MaxIt 
     for subit=1:MaxSubIt
        newpop=repmat(empty_individual,nPop,nMove);
        for i=1:nPop
            for j=1:nMove
                % Create Neighbor
                newpop(i,j).Position=Mutate(pop(i).Position,mu,sigma,VarMin,VarMax);
                % Evaluation
                newpop(i,j).Cost=Data_i.fobj(newpop(i,j).Position);
            end
        end
        newpop=newpop(:);
        % Sort Neighbors
        [~, SortOrder]=sort([newpop.Cost]);
        newpop=newpop(SortOrder);
        for i=1:nPop
            if newpop(i).Cost<=pop(i).Cost
                pop(i)=newpop(i);
            else
                DELTA=(newpop(i).Cost-pop(i).Cost)/pop(i).Cost;
                P=exp(-DELTA/T);
                if rand<=P
                    pop(i)=newpop(i);
                end
            end
            % Update Best Solution Ever Found
            if pop(i).Cost<BestSol.Cost
                BestSol.Cost=pop(i).Cost;
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestSol.Cost;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
     end   
     IterCurve(it)=BestSol.Cost;
     T=alpha*T;
     sigma = 0.98*sigma;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function y=Mutate(x,mu,sigma,VarMin,VarMax)
    A=(rand(size(x))<=mu);
    J=find(A==1);
    y=x;
    y(J)=x(J)+sigma*randn(size(J));
    % Clipping
    y=max(y,VarMin);
    y=min(y,VarMax);  
end