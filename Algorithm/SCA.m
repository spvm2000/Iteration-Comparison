function Data_o = SCA(Data_i,Data_o)   
    a = 2; 
    t_c=clock;                             
    cnt=1;                               
    X = Data_i.X;
    fitness = Data_i.F_value;
    IterCurve=zeros(1,Data_i.maxIter);
    [SortFitness,indexSort] = sort(fitness);
    Best_Pos=indexSort(1);
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Pos;
    gBest = X(indexSort(1),:);      
    gBestFitness = SortFitness(1);  
    for t = 1:Data_i.maxIter
           r1 = a - t*(a/Data_i.maxIter);
        for i = 1:Data_i.pop
            for j = 1:Data_i.dim
                r2 = rand()*(2*pi);
                r3 = 2*rand();
                r4 = rand();
                if r4<0.5
                    X(i,j) = X(i,j)+(r1*sin(r2)*abs(r3*gBest(j)-X(i,j))); 
                else
                    X(i,j) = X(i,j)+(r1*cos(r2)*abs(r3*gBest(j)-X(i,j))); 
                end
            end
             X(i,:) = BoundaryCheck(X(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        end
         for i = 1:Data_i.pop
             fitness(i) = Data_i.fobj(X(i,:));
            if  fitness(i)<gBestFitness 
                gBestFitness= fitness(i); 
                gBest=X(i,:);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
         end
         Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=gBestFitness;  
         IterCurve(t)=gBestFitness;
end

Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t; 
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end