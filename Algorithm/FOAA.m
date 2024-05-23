function Data_o = FOAA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
cg_curve=zeros(1,Data_i.maxIter);
%% Problem Definition
maxt = Data_i.maxIter;      
n = Data_i.pop;            
dim=Data_i.dim;
Sol = Data_i.X;               
Fitness = Data_i.F_value;       
lb = Data_i.lb;
ub = Data_i.ub;

% Initialize the original position
% for i = 1:n
%     X(i,:) = lb+(ub-lb).*rand(1,dim); % the position of X axis
%     Y(i,:) = lb+(ub-lb).*rand(1,dim); % the position of Y axis
%     D(i,:) = (X(i,:).^2 + Y(i,:).^2).^0.5; % Caculate the distance
%     Sol(i,:) = 1./D(i,:); % the solution set
%     Fitness(i) = fun(Sol(i,:)); % Caculate the fitness
% end

[bestSmell,index] = min(Fitness); % Get the min fitness and its index
% new_X = X(index,:); % the X axis of min fitness
% new_Y = Y(index,:); % the Y axis of min fitness
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=bestSmell;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
Smellbest = bestSmell;

% Start main loop
for t = 1:maxt
    for i = 1:n
        % Refer to the process of initializing
        if t==1
            X(i,:) = lb + (ub - lb).*rand(1,dim);
            Y(i,:) = lb + (ub - lb).*rand(1,dim);
        else
            X(i,:) = new_X + (ub - lb).*rand(1,dim);
            Y(i,:) = new_Y + (ub - lb).*rand(1,dim);
        end
        D(i,:) = (X(i,:).^2 + Y(i,:).^2).^0.5;
        Sol(i,:) = 1./D(i,:);
        Fitness(i) = Data_i.fobj(Sol(i,:));
    end
    [bestSmell,index] = min(Fitness);
    new_X = X(index,:);
    new_Y = Y(index,:);
    % If the new value is smaller than the best value,update the best value
    if (bestSmell < Smellbest)
        Smellbest = bestSmell;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    end
    cg_curve(t) = Smellbest;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Smellbest;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=cg_curve;
end     