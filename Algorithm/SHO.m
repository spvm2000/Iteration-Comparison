function Data_o = SHO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
hyena_pos = Data_i.X;                        
hyena_fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
Iteration=1;        Max_iterations=Data_i.maxIter;    fitness=Data_i.fobj;  N=Data_i.pop;
while Iteration<=Max_iterations
    for i=1:size(hyena_pos,1)
        hyena_pos(i,:)=BoundaryCheck(hyena_pos(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        hyena_fitness(1,i)=fitness(hyena_pos(i,:));  
    end
    
    if Iteration==1
        [fitness_sorted FS]=sort(hyena_fitness);
        sorted_population=hyena_pos(FS,:);
        best_hyenas=sorted_population;
        best_hyena_fitness=fitness_sorted;
        
    else
        double_population=[pre_population;best_hyenas];
        double_fitness=[pre_fitness best_hyena_fitness];
        [double_fitness_sorted FS]=sort(double_fitness);
        double_sorted_population=double_population(FS,:);
        fitness_sorted=double_fitness_sorted(1:N);
        sorted_population=double_sorted_population(1:N,:);
        best_hyenas=sorted_population;
        best_hyena_fitness=fitness_sorted;
    end
    NOH=noh(best_hyena_fitness);
    Best_hyena_score=fitness_sorted(1);
    Best_hyena_pos=sorted_population(1,:);
    pre_population=hyena_pos;
    pre_fitness=hyena_fitness;
    
    a=5-Iteration*((5)/Max_iterations);
    HYE=0;
    CV=0;
    for i=1:size(hyena_pos,1)
        for j=1:size(hyena_pos,2)
             for k=1:NOH
                HYE=0;
                r1=rand();
                r2=rand(); 
                Var1=2*a*r1-a; 
                Var2=2*r2; 
                distance_to_hyena=abs(Var2*sorted_population(k)-hyena_pos(i,j));
                HYE=sorted_population(k)-Var1*distance_to_hyena;
                CV=CV+HYE;        
                distance_to_hyena=0;
             end
            hyena_pos(i,j)=(CV/(NOH+1));
            CV=0;
        end
    end
    IterCurve(Iteration)=Best_hyena_score;
    Iteration=Iteration+1; 
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iteration;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function X=noh(best_hyena_fitness)
    min = 0.5;
    max = 1;
    count=0;
    M=(max-min).*rand(1,1) + min;
    M=M+best_hyena_fitness(1);
    for i=2:numel(best_hyena_fitness)
        if M<=best_hyena_fitness(i)
            count=count+1;
        end
    end
     
    X=count;
    clear count
    clear M
    clear px
end


