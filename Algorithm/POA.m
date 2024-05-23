function Data_o = POA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 
%% Problem Definition
Max_iterations = Data_i.maxIter;    
SearchAgents = Data_i.pop;          
dimension = Data_i.dim;             
X = Data_i.X;                       
fit = Data_i.F_value;               
lowerbound=Data_i.lb;                             % Lower limit for variables
upperbound=Data_i.ub;                             % Upper limit for variables
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);
Curve=zeros(1,Data_i.maxIter);
for t=1:Max_iterations
    
    %% UPDATE location of food
    
%     X_FOOD=[];
    k=randperm(SearchAgents,1);
    X_FOOD=X(k,:);
    F_FOOD=fit(k);
    
    %%
    for i=1:SearchAgents
        
        %% PHASE 1: Moving towards prey (exploration phase)
        I=round(1+rand(1,1));
        if fit(i)> F_FOOD
            X_new=X(i,:)+ rand(1,1).*(X_FOOD-I.* X(i,:)); %Eq(4)
        else
            X_new=X(i,:)+ rand(1,1).*(X(i,:)-1.*X_FOOD); %Eq(4)
        end
        X_new = max(X_new,lowerbound);
        X_new = min(X_new,upperbound);
        
        % Updating X_i using (5)
        f_new = Data_i.fobj(X_new);
        if f_new <= fit (i)
            X(i,:) = X_new;
            fit (i)=f_new;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
        %% END PHASE 1: Moving towards prey (exploration phase)
        
        %% PHASE 2: Winging on the water surface (exploitation phase)
        X_new=X(i,:)+0.2*(1-t/Max_iterations).*(2*rand(1,dimension)-1).*X(i,:);% Eq(6)
        X_new= max(X_new,lowerbound);
        X_new= min(X_new,upperbound);
        
        % Updating X_i using (7)
        f_new = Data_i.fobj(X_new);
        if f_new <= fit (i)
            X(i,:) = X_new;
            fit (i)=f_new;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
        %% END PHASE 2: Winging on the water surface (exploitation phase)
    end
    
    %% update the best condidate solution
    [best , ~]=min(fit);
    if t==1                            
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
    end
    Curve(t) = fbest;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fbest;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Curve;
end

