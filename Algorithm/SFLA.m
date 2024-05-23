function Data_o = SFLA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                            
X = Data_i.X;                        
Ffun=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
CostFunction =Data_i.fobj;      nVar = Data_i.dim;      VarSize = [1 nVar];     VarMin = Data_i.lb(1);
VarMax = Data_i.ub(1);      MaxIt = Data_i.maxIter;     

nMemeplex = 5;                              % Number of Memeplexes
nPopMemeplex = Data_i.pop/nMemeplex;        % Memeplex Size
nPop=Data_i.pop;

I = reshape(1:nPop, nMemeplex, []);
% FLA Parameters
fla_params.q = max(round(0.3*nPopMemeplex),2);   % Number of Parents
fla_params.alpha = 3;   % Number of Offsprings
fla_params.beta = 5;    % Maximum Number of Iterations
fla_params.sigma = 2;   % Step Size
fla_params.CostFunction = CostFunction;
fla_params.VarMin = VarMin;
fla_params.VarMax = VarMax;

% Empty Individual Template
empty_individual.Position = [];
empty_individual.Cost = [];
pop = repmat(empty_individual, nPop, 1);
for i=1:nPop
    pop(i).Position = X(i,:);
    pop(i).Cost = Ffun(i);
end

pop = SortPopulation(pop);
BestSol = pop(1);
for it = 1:MaxIt
    fla_params.BestSol = BestSol;

    % Initialize Memeplexes Array
    Memeplex = cell(nMemeplex, 1);
    
    % Form Memeplexes and Run FLA
    for j = 1:nMemeplex
        % Memeplex Formation
        Memeplex{j} = pop(I(j,:));
        % Run FLA
        Memeplex{j} = RunFLA(Memeplex{j}, fla_params);
        
        % Insert Updated Memeplex into Population
        pop(I(j,:)) = Memeplex{j};
    end
    
    % Sort Population
    [pop,sortt] = SortPopulation(pop);
    
    % Update Best Solution Ever Found
    BestSol = pop(1);
    
    % Store Best Cost Ever Found
    IterCurve(it) = BestSol.Cost;
    if BestSol.Cost<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestSol.Cost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=sortt(1);
    end    
end
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function pop = RunFLA(pop, params)
    %% FLA Parameters
    q = params.q;           % Number of Parents
    alpha = params.alpha;   % Number of Offsprings
    beta = params.beta;     % Maximum Number of Iterations
    sigma = params.sigma;
    CostFunction = params.CostFunction;
    VarMin = params.VarMin;
    VarMax = params.VarMax;
    VarSize = size(pop(1).Position);
    BestSol = params.BestSol;
    
    nPop = numel(pop);      % Population Size
    P = 2*(nPop+1-(1:nPop))/(nPop*(nPop+1));    % Selection Probabilities
    
    % Calculate Population Range (Smallest Hypercube)
    LowerBound = pop(1).Position;
    UpperBound = pop(1).Position;
    for i = 2:nPop
        LowerBound = min(LowerBound, pop(i).Position);
        UpperBound = max(UpperBound, pop(i).Position);
    end
    
    %% FLA Main Loop

    for it = 1:beta
        
        % Select Parents
        L = RandSample(P,q);
        B = pop(L);
        
        % Generate Offsprings
        for k=1:alpha
            
            % Sort Population
            [B, SortOrder] = SortPopulation(B);
            L = L(SortOrder);
            
            % Flags
            ImprovementStep2 = false;
            Censorship = false;
            
            % Improvement Step 1
            NewSol1 = B(end);
            Step = sigma*rand(VarSize).*(B(1).Position-B(end).Position);
            NewSol1.Position = B(end).Position + Step;
            if IsInRange(NewSol1.Position, VarMin, VarMax)
                NewSol1.Cost = CostFunction(NewSol1.Position);
                if NewSol1.Cost<B(end).Cost
                    B(end) = NewSol1;
                else
                    ImprovementStep2 = true;
                end
            else
                ImprovementStep2 = true;
            end
            
            % Improvement Step 2
            if ImprovementStep2
                NewSol2 = B(end);
                Step = sigma*rand(VarSize).*(BestSol.Position-B(end).Position);
                NewSol2.Position = B(end).Position + Step;
                if IsInRange(NewSol2.Position, VarMin, VarMax)
                    NewSol2.Cost = CostFunction(NewSol2.Position);
                    if NewSol2.Cost<B(end).Cost
                        B(end) = NewSol2;
                    else
                        Censorship = true;
                    end
                else
                    Censorship = true;
                end
            end
                
            % Censorship
            if Censorship
                B(end).Position = unifrnd(LowerBound, UpperBound);
                B(end).Cost = CostFunction(B(end).Position);
            end
            
        end
        
        % Return Back Subcomplex to Main Complex
        pop(L) = B;
        
    end
end

function L = RandSample(P, q, replacement)

    if ~exist('replacement','var')
        replacement = false;
    end

    L = zeros(q,1);
    for i=1:q
        L(i) = randsample(numel(P), 1, true, P);
        if ~replacement
            P(L(i)) = 0;
        end
    end

end

function b = IsInRange(x, VarMin, VarMax)
    b = all(x>=VarMin) && all(x<=VarMax);
end

function [pop, SortOrder] = SortPopulation(pop)

    % Get Costs
    Costs = [pop.Cost];
    % Sort the Costs Vector
    [~, SortOrder]=sort(Costs);
    % Apply the Sort Order to Population
    pop = pop(SortOrder);
end