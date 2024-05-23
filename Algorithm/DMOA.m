function Data_o = DMOA(Data_i,Data_o)
t_c=clock;                            
cnt=1;                              
%% Problem Definition
nVar=Data_i.dim;             
VarSize=[1 nVar];   
VarMin=Data_i.lb;         % Decision Variables Lower Bound
VarMax=Data_i.ub;         % Decision Variables Upper Bound
MaxIt=Data_i.maxIter;              % Maximum Number of Iterations
nPop=Data_i.pop;               % Population Size (Family Size)

nBabysitter= 3;         % Number of babysitters
nAlphaGroup=nPop-nBabysitter;         % Number of Alpha group
nScout=nAlphaGroup;         % Number of Scouts
L=round(0.6*nVar*nBabysitter); % Babysitter Exchange Parameter 
peep=2;             % Alpha female vocalization 
% Empty Mongoose Structure
empty_mongoose.Position=[];
empty_mongoose.Cost=[];

% Initialize Population Array
pop=repmat(empty_mongoose,nAlphaGroup,1);

% Initialize Best Solution Ever Found
BestSol=inf;
tau=inf;
Iter=1;
sm=inf(nAlphaGroup,1);

% Create Initial Population
X=Data_i.X;
fit=Data_i.F_value;
for i=1:nAlphaGroup
    pop(i).Position=X(i,:);
    pop(i).Cost=fit(i);
    if pop(i).Cost<=BestSol
        BestSol=pop(i).Cost;
    end
end
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);
% Abandonment Counter
C=zeros(nAlphaGroup,1);
CF=(1-Iter/MaxIt)^(2*Iter/MaxIt);

% Array to Hold Best Cost Values
IterCurve=zeros(1,Data_i.maxIter);

%% DMOA Main Loop

for it=1:MaxIt
    % Alpha group
     F=zeros(nAlphaGroup,1);
     MeanCost = mean([pop.Cost]);
    for i=1:nAlphaGroup
        % Calculate Fitness Values and Selection of Alpha
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
      % Foraging led by Alpha female
    for m=1:nAlphaGroup
        
        % Select Alpha female
        i=RouletteWheelSelection(P,nAlphaGroup);
        
        % Choose k randomly, not equal to Alpha
        K=[1:i-1 i+1:nAlphaGroup]; 
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);
        
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newpop.Cost=Data_i.fobj(newpop.Position);
        
        % Comparision
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
        else
            C(i)=C(i)+1;
        end
        
    end   
    
    % Scout group
    for i=1:nScout
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);
        
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newpop.Cost=Data_i.fobj(newpop.Position);
        
        % Sleeping mould
        sm(i)=(newpop.Cost-pop(i).Cost)/max(newpop.Cost,pop(i).Cost);
        
        % Comparision
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
        else
            C(i)=C(i)+1;
        end
        
    end    
    % Babysitters
    for i=1:nBabysitter
        if C(i)>=L
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            pop(i).Cost=Data_i.fobj(pop(i).Position);
            C(i)=0;
        end
    end    
     % Update Best Solution Ever Found
    for i=1:nAlphaGroup
        if pop(i).Cost<=BestSol
            BestSol=pop(i).Cost;
        end
    end    
        
   % Next Mongoose Position
   newtau=mean(sm);
   for i=1:nScout
        M=(pop(i).Position.*sm(i))/pop(i).Position;
        if newtau>tau
           newpop.Position=pop(i).Position-CF*phi*rand.*(pop(i).Position-M);
        else
           newpop.Position=pop(i).Position+CF*phi*rand.*(pop(i).Position-M);
        end
        tau=newtau;
   end
       
   % Update Best Solution Ever Found
    for i=1:nAlphaGroup
        if pop(i).Cost<=BestSol
            BestSol=pop(i).Cost;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    
    IterCurve(it)=BestSol;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = BestSol;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                          
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
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