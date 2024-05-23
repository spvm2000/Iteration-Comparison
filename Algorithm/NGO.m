function Data_o = NGO(Data_i,Data_o)
rand_num=[];                        
t=clock;                             
cnt=1;                               
rand_num=[];                         
X = Data_i.X;
x_temp=X;
fitness=Data_i.F_value;                         
IterCurve=zeros(1,Data_i.maxIter);
waring_bird_per=0.2;                            
fitness_temp=fitness;
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
for gen=1:Data_i.maxIter
    [SortFitness,indexSort] = sort(fitness_temp);               
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=indexSort(1);
    gBest = x_temp(indexSort(1),:);                             
    [gworstFitness,gworstFitness_index]=max(fitness_temp);      
    gworst=x_temp(gworstFitness_index,:);                       
    for i=1:Data_i.pop
        rand_index=randperm(Data_i.pop,1);                                    
        p=mean(X(rand_index,:));                                              
        F_p=mean(fitness(rand_index));                                        
        I=randi(2);
        rand_num(2,gen)=I;
        if fitness(i) > F_p
            x_temp(i,:)=X(i,:)+rand(1,Data_i.dim).*(p-I*X(i,:));
        else
            x_temp(i,:)=X(i,:)+rand(1,Data_i.dim).*(X(i,:)-p);
        end
        
        x_temp(i,:)=max(x_temp(i,:),Data_i.lb.*ones(1,Data_i.dim));
        x_temp(i,:)=min(x_temp(i,:),Data_i.ub.*ones(1,Data_i.dim));
        
        fitness_temp(i)=Data_i.fobj(x_temp(i,:));
        if fitness_temp(i)<fitness(i)
            X(i,:)=x_temp(i,:);
            fitness(i)=fitness_temp(i);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fitness(i)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness(i);   
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                gBest=X(i,:);
            end    
        end
        R=0.02;
        x_temp(i,:)=X(i,:)+(-R+2*R*rand(1,Data_i.dim)).*X(i,:);
        x_temp(i,:)=max(x_temp(i,:),Data_i.lb.*ones(1,Data_i.dim));
        x_temp(i,:)=min(x_temp(i,:),Data_i.ub.*ones(1,Data_i.dim));
        
        fitness_temp(i)=Data_i.fobj(x_temp(i,:));
        if fitness_temp(i)<fitness(i)
            X(i,:)=x_temp(i,:);
            fitness(i)=fitness_temp(i); 
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fitness(i)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness(i);   
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;  
                gBest=X(i,:);
            end
        end
    end

    [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
    IterCurve(gen)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=gen;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end