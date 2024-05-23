function Data_o = ABC(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
CostFunction=Data_i.fobj;        % Cost Function

nVar=Data_i.dim;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin= Data_i.lb;         % Decision Variables Lower Bound
VarMax= Data_i.ub;         % Decision Variables Upper Bound

%% ABC Settings

MaxIt=Data_i.maxIter;              % Maximum Number of Iterations

nPop=Data_i.pop;               % Population Size (Colony Size)

nOnlooker=nPop;         % Number of Onlooker Bees

L=round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)

a=1;                    % Acceleration Coefficient Upper Bound

%% Initialization
% Create Initial Population
Position=Data_i.X;
Cost=Data_i.F_value;
[BestSol,index]=min(Cost);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestSol;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;

% Abandonment Counter
C=zeros(nPop,1);

% Array to Hold Best Cost Values
BestCost=zeros(1,Data_i.maxIter);

%% ABC Main Loop

for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee_Position=Position(i,:)+phi.*(Position(i,:)-Position(k,:));
        
        % Evaluation
        newbee_Cost=CostFunction(newbee_Position);
        
        % Comparision
        if newbee_Cost<=Cost(i)
            Position(i,:)=newbee_Position;
            Cost(i)=newbee_Cost;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    MeanCost = mean(Cost);
    for i=1:nPop
        F(i) = exp(-Cost(i)/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P,Data_i.pop);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        temp=numel(K);
        if temp<=0
            temp=1;
        end
        k=K(randi([1 temp]));
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee_Position=Position(i,:)+phi.*(Position(i,:)-Position(k,:));
        
        % Evaluation
        newbee_Cost=CostFunction(newbee_Position);
        
        % Comparision
        if newbee_Cost<=Cost(i)
            Position(i,:)=newbee_Position;
            Cost(i)=newbee_Cost;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            Position(i,:)=unifrnd(VarMin,VarMax,VarSize);
            Cost(i)=CostFunction(Position(i,:));
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if Cost(i)<=BestSol
            BestSol=Cost(i);
        end
    end
    
    BestCost(it)=BestSol;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = BestSol;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=BestCost;
end   

function i=RouletteWheelSelection(P,pop)
    r=rand;
    C=cumsum(P);
    i=find(r<=C,1,'first');
    if i<1 
        i=randi(pop);
    end
    if length(i)<1
        i=randi(pop);
    end
end