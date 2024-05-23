function Data_o = SSA(Data_i,Data_o)
rand_num=[];                         
t=clock;                                                            
X = Data_i.X;
fitness = Data_i.F_value;
X_temp=X;
fitness_temp=fitness;
IterCurve=zeros(1,Data_i.maxIter);
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
for gen=1:Data_i.maxIter
    [SortFitness,indexSort] = sort(fitness_temp);               
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=indexSort(1);
    gBest = X_temp(indexSort(1),:);                            
    [gworstFitness,gworstFitness_index]=max(fitness_temp);      
    gworst=X_temp(gworstFitness_index,:);                       
    [~,whole_fit_index]=sort(fitness);                          
    
    r2=rand(1);
    rand_num(1,gen)=r2;
    st=0.5+(1-0.5)*rand(1,1);
    rand_num(2,gen)=st;
    for a=1:round(Data_i.pop*0.2)
        if r2<st
            X_temp(a,:)=X_temp(a,:) * exp(-a/(randn(1)*Data_i.maxIter));
        else
            X_temp(a,:)=X_temp(a,:)+randn(1).*ones(1,Data_i.dim);
        end
        X_temp(a,:)=BoundaryCheck(X_temp(a,:),Data_i.ub,Data_i.lb,Data_i.dim);
        fitness_temp(a) = Data_i.fobj(X_temp(a,:));
    end
    for a=round(Data_i.pop*0.2)+1:Data_i.pop
        if a>Data_i.pop/2
            X_temp(a,:)=randn(1)*exp((gworst-X_temp(a,:))/a^2);
        else
            A=floor(rand(1,Data_i.dim)*2)*2-1;
            A_jia=A'*(A*A')^-1;
            X_temp(a,:)=gBest+abs(X_temp(a,:)-gBest)*A_jia.*ones(1,Data_i.dim);
        end
        X_temp(a,:)=BoundaryCheck(X_temp(a,:),Data_i.ub,Data_i.lb,Data_i.dim);
        fitness_temp(a) = Data_i.fobj(X_temp(a,:));
    end

    [SortFitness,indexSort] = sort(fitness_temp);
    gBest = X_temp(indexSort(1),:);      
    gBestFitness = SortFitness(1);       
    gworst=X_temp(indexSort(Data_i.pop),:);     
    gworstFitness=SortFitness(end);      
    gBestFitness = SortFitness(1);       

    c=randperm(numel(fitness_temp));       
    b=c(Data_i.pop-round(Data_i.pop*0.2,0):Data_i.pop);         
    for a=1:length(b)
        if fitness_temp(b(a))>gBestFitness
            X_temp(b(a),:)=gBest+randn(1)*abs(X_temp(b(a),:)-gBest);
        elseif fitness_temp(b(a))==gworstFitness
            X_temp(b(a),:)=X_temp(b(a),:)+(2*rand(1)-1)*((abs(X_temp(b(a),:)-gworst)/(fitness_temp(b(a))-gworstFitness)+0.000001));
        end    
        X_temp(a,:)=BoundaryCheck(X_temp(a,:),Data_i.ub,Data_i.lb,Data_i.dim);
        fitness_temp(a) = Data_i.fobj(X_temp(a,:));
    end
    for Q_J=1:Data_i.pop
        if fitness(Q_J)>fitness_temp(Q_J)
            fitness(Q_J)=fitness_temp(Q_J);
            X(Q_J,:)=X_temp(Q_J,:);
            if fitness(Q_J)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness(Q_J);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Q_J;
            end
        end
    end
    IterCurve(gen)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=gen;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end