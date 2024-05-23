function Data_o = GEO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  

%% initialization
PopulationSize = Data_i.pop;
MaxIterations = Data_i.maxIter;
nvars=Data_i.dim;
x = Data_i.X;
FitnessScores = Data_i.F_value;
IterCurve  = zeros(1,Data_i.maxIter);

% solver-specific initialization
FlockMemoryX = x;
[FlockMemoryF,~] = min(FitnessScores);


options.AttackPropensity = [0.5 ,   2];
options.CruisePropensity = [1   , 0.5];
AttackPropensity = linspace (options.AttackPropensity(1), options.AttackPropensity(2), MaxIterations);
CruisePropensity = linspace (options.CruisePropensity(1), options.CruisePropensity(2), MaxIterations);
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(FitnessScores);

%% main loop

for CurrentIteration = 1 : MaxIterations
	
	% prey selection (one-to-one mapping)
	DestinationEagle = randperm (PopulationSize)';
	
	% calculate AttackVectorInitial (Eq. 1 in paper)
	AttackVectorInitial = FlockMemoryX (DestinationEagle,:) - x;
	
	% calculate Radius
	Radius = VecNorm (AttackVectorInitial, 2, 2);
	
	% determine converged and unconverged eagles
	ConvergedEagles = sum (Radius,2) == 0;
	UnconvergedEagles = ~ ConvergedEagles;
	
	% initialize CruiseVectorInitial
	CruiseVectorInitial = 2 .* rand (PopulationSize, nvars) - 1; % [-1,1]
	
	% correct vectors for converged eagles
	AttackVectorInitial (ConvergedEagles, :) = 0;
	CruiseVectorInitial (ConvergedEagles, :) = 0;
	
	% determine constrained and free variables
	for i1 = 1 : PopulationSize
		if UnconvergedEagles (i1)
			vConstrained = false ([1, nvars]); % mask
			idx = datasample (find(AttackVectorInitial(i1,:)), 1, 2);
			vConstrained (idx) = 1;
			vFree = ~vConstrained;
			CruiseVectorInitial (i1,idx) = - sum(AttackVectorInitial(i1,vFree).*CruiseVectorInitial(i1,vFree),2) ./ (AttackVectorInitial(i1,vConstrained)); % (Eq. 4 in paper)
		end
	end
	
	% calculate unit vectors
	AttackVectorUnit = AttackVectorInitial ./ VecNorm (AttackVectorInitial, 2, 2);
	CruiseVectorUnit = CruiseVectorInitial ./ VecNorm (CruiseVectorInitial, 2, 2);
	
	% correct vectors for converged eagles
	AttackVectorUnit(ConvergedEagles,:) = 0;
	CruiseVectorUnit(ConvergedEagles,:) = 0;
	
	% calculate movement vectors
	AttackVector = rand (PopulationSize, 1) .* AttackPropensity(CurrentIteration) .* Radius .* AttackVectorUnit; % (first term of Eq. 6 in paper)
	CruiseVector = rand (PopulationSize, 1) .* CruisePropensity(CurrentIteration) .* Radius .* CruiseVectorUnit; % (second term of Eq. 6 in paper)
	StepVector = AttackVector + CruiseVector;
	
	% calculate new x
	x = x + StepVector;
	
	% enforce bounds
	lbExtended = repmat (Data_i.lb,[PopulationSize,1]);
	ubExtended = repmat (Data_i.ub,[PopulationSize,1]);
	
	lbViolated = x < lbExtended;
	ubViolated = x > ubExtended;
	
	x (lbViolated) = lbExtended (lbViolated);
	x (ubViolated) = ubExtended (ubViolated);

	% calculate fitness
    for j=1:PopulationSize
        FitnessScores = Data_i.fobj (x(j,:));
        % update memory
        if  FitnessScores < FlockMemoryF
            FlockMemoryF  = FitnessScores ;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=j;
        end
    end
	IterCurve (CurrentIteration) =FlockMemoryF;
	Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = FlockMemoryF;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=CurrentIteration;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function N = VecNorm (A, p, dim)  
    N = nthroot(sum(A.^p,dim),p); 
end

