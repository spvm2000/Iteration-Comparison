function Data_o = TSOO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
voltages = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
best_score=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
best_voltage=voltages(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
ub=Data_i.ub;           lb=Data_i.lb;
fobj=Data_i.fobj;
l=0;  L=Data_i.maxIter;
while l<L
    for i=1:size(voltages,1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=voltages(i,:)>ub;
        Flag4lb=voltages(i,:)<lb;
        voltages(i,:)=(voltages(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
               
        fitness=fobj(voltages(i,:));% Calculate objective function for each search agent
        
        % Update the best
        if fitness<best_score 
            best_score=fitness; 
            best_voltage=voltages(i,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
        
    end
    
    t=2-l*(2/L); 
    K=1;   % K is a real can be 0, 1, 2,....
    % Update the voltage of search agents 
    for i=1:size(voltages,1)
        r1=rand();r2=rand(); 
        r3 = rand();
        T=2*t*r1-t;  
        C1=K*r2*t+1;              
        for j=1:size(voltages,2)
            if r3<0.5
                voltages(i,j)=best_voltage(j)+exp(-T)*(voltages(i,j)-C1*best_voltage(j));
            elseif r3>=0.5
              voltages(i,j)=best_voltage(j)+exp(-T)*(cos(T*2*pi)+sin(T*2*pi))*abs(voltages(i,j)-C1*best_voltage(j));
            end           
        end
    end
    l=l+1;
    IterCurve(l)=best_score;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end