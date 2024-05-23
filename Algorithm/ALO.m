function Data_o = ALO(Data_i,Data_o)                        
ti=clock;                             
antlion_position = Data_i.X;                        
antlions_fitness=Data_i.F_value;               
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
% 专属变量定义
N=Data_i.pop;     fobj=Data_i.fobj;       lb=Data_i.lb; ub=Data_i.ub;  dim=Data_i.dim;
ant_position=ALO_init(N,dim,ub,lb);     Max_iter=Data_i.maxIter;

Sorted_antlions=zeros(N,dim);
Elite_antlion_position=zeros(1,dim);
Elite_antlion_fitness=inf;
antlions_fitness=zeros(1,N);
ants_fitness=zeros(1,N);

[sorted_antlion_fitness,sorted_indexes]=sort(antlions_fitness);
for newindex=1:N
     Sorted_antlions(newindex,:)=antlion_position(sorted_indexes(newindex),:);
end

Elite_antlion_position=Sorted_antlions(1,:);
Elite_antlion_fitness=sorted_antlion_fitness(1);
Current_iter=1; 
while Current_iter<Max_iter+1
    % This for loop simulate random walks
    for i=1:size(ant_position,1)
        % Select ant lions based on their fitness (the better anlion the higher chance of catching ant)
        Rolette_index=RouletteWheelSelection(1./sorted_antlion_fitness);
        if Rolette_index==-1  
            Rolette_index=1;
        end
      
        % RA is the random walk around the selected antlion by rolette wheel
        RA=Random_walk_around_antlion(dim,Max_iter,lb,ub, Sorted_antlions(Rolette_index,:),Current_iter);
        
        % RA is the random walk around the elite (best antlion so far)
        [RE]=Random_walk_around_antlion(dim,Max_iter,lb,ub, Elite_antlion_position(1,:),Current_iter);
        
        ant_position(i,:)= (RA(Current_iter,:)+RE(Current_iter,:))/2; % Equation (2.13) in the paper          
    end
    
    for i=1:size(ant_position,1)  
        
        % Boundar checking (bring back the antlions of ants inside search
        % space if they go beyoud the boundaries
        Flag4ub=ant_position(i,:)>ub;
        Flag4lb=ant_position(i,:)<lb;
        ant_position(i,:)=(ant_position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        ants_fitness(1,i)=fobj(ant_position(i,:));        
       
    end
    
    % Update antlion positions and fitnesses based of the ants (if an ant 
    % becomes fitter than an antlion we assume it was cought by the antlion  
    % and the antlion update goes to its position to build the trap)
    double_population=[Sorted_antlions;ant_position];
    double_fitness=[sorted_antlion_fitness ants_fitness];
        
    [double_fitness_sorted I]=sort(double_fitness);
    double_sorted_population=double_population(I,:);
        
    antlions_fitness=double_fitness_sorted(1:N);
    Sorted_antlions=double_sorted_population(1:N,:);
        
    % Update the position of elite if any antlinons becomes fitter than it
    if antlions_fitness(1)<Elite_antlion_fitness 
        Elite_antlion_position=Sorted_antlions(1,:);
        Elite_antlion_fitness=antlions_fitness(1);
        if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>antlions_fitness(1)
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=antlions_fitness(1);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=I(1);
        end
    end
      
    % Keep the elite in the population
    Sorted_antlions(1,:)=Elite_antlion_position;
    antlions_fitness(1)=Elite_antlion_fitness;
  
    % Update the convergence curve
    IterCurve(Current_iter)=Elite_antlion_fitness;
    Current_iter=Current_iter+1; 
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Current_iter;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function choice = RouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand() * accumulation(end);
  chosen_index = -1;
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end  

function [RWs]=Random_walk_around_antlion(Dim,max_iter,lb, ub,antlion,current_iter)
    if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
        lb=ones(1,Dim)*lb;
        ub=ones(1,Dim)*ub;
    end
    
    if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
        lb=lb';
        ub=ub';
    end
    
    I=1; % I is the ratio in Equations (2.10) and (2.11)
    
    if current_iter>max_iter/10
        I=1+100*(current_iter/max_iter);
    end
    
    if current_iter>max_iter/2
        I=1+1000*(current_iter/max_iter);
    end
    
    if current_iter>max_iter*(3/4)
        I=1+10000*(current_iter/max_iter);
    end
    
    if current_iter>max_iter*(0.9)
        I=1+100000*(current_iter/max_iter);
    end
    
    if current_iter>max_iter*(0.95)
        I=1+1000000*(current_iter/max_iter);
    end
    
    
    % Dicrease boundaries to converge towards antlion
    lb=lb/(I); % Equation (2.10) in the paper
    ub=ub/(I); % Equation (2.11) in the paper
    
    % Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
    if rand<0.5
        lb=lb+antlion; % Equation (2.8) in the paper
    else
        lb=-lb+antlion;
    end
    
    if rand>=0.5
        ub=ub+antlion; % Equation (2.9) in the paper
    else
        ub=-ub+antlion;
    end
    
    % This function creates n random walks and normalize accroding to lb and ub
    % vectors 
    for i=1:Dim
        X = [0 cumsum(2*(rand(max_iter,1)>0.5)-1)']; % Equation (2.1) in the paper
        %[a b]--->[c d]
        a=min(X);
        b=max(X);
        c=lb(i);
        d=ub(i);      
        X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
        RWs(:,i)=X_norm;
    end
end