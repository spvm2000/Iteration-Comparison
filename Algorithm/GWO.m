function Data_o = GWO(Data_i,Data_o)
    t_c=clock;                             
    cnt=1;                               
    IterCurve=zeros(1,Data_i.maxIter);
    Alpha_pos=zeros(1,Data_i.dim);
    Alpha_score=inf; 

    Beta_pos=zeros(1,Data_i.dim);
    Beta_score=inf; 

    Delta_pos=zeros(1,Data_i.dim);
    Delta_score=inf; 
    Positions = Data_i.X;
    fitness = Data_i.F_value;
    [SortFitness,indexSort] = sort(fitness);
    Best_Pos=indexSort(1);
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Pos;
    Alpha_pos = Positions(indexSort(1),:);
    Alpha_score = SortFitness(1);
    Beta_pos = Positions(indexSort(2),:);
    Beta_score = SortFitness(2);
    Delta_pos = Positions(indexSort(3),:);
    Delta_score = SortFitness(3);
    gBest = Alpha_pos;
    gBestFitness = Alpha_score;  
    for t = 1:Data_i.maxIter
        a=2-t*((2)/Data_i.maxIter);      
        for i = 1:Data_i.pop
            for j = 1:Data_i.dim
                r1=rand();   
                r2=rand();               
                A1=2*a*r1-a; 
                C1=2*r2;     
                D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); 
                X1=Alpha_pos(j)-A1*D_alpha; 
              
                r1=rand();
                r2=rand();      
                A2=2*a*r1-a;
                C2=2*r2;           
                D_beta=abs(C2*Beta_pos(j)-Positions(i,j));  
                X2=Beta_pos(j)-A2*D_beta; 
                
                r1=rand();
                r2=rand();
                A3=2*a*r1-a; 
                C3=2*r2; 
                D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); 
                X3=Delta_pos(j)-A3*D_delta; 
                Positions(i,j)=(X1+X2+X3)/3;
            end
             Positions(i,:) = BoundaryCheck(Positions(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        end
         for i = 1:Data_i.pop
             fitness(i) = Data_i.fobj(Positions(i,:));
             
            if  fitness(i)<Alpha_score 
                Alpha_score= fitness(i); 
                Alpha_pos=Positions(i,:);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
            if  fitness(i)>Alpha_score &&  fitness(i)<Beta_score 
                Beta_score= fitness(i); 
                Beta_pos=Positions(i,:);
            end
            if  fitness(i)>Alpha_score &&  fitness(i)>Beta_score &&  fitness(i)<Delta_score 
                Delta_score= fitness(i); 
                Delta_pos=Positions(i,:);
            end
         end
        gBest = Alpha_pos;
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Alpha_score;  
        Best_fitness = Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
        IterCurve(t)=Best_fitness;
    end
        Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);         
        Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end