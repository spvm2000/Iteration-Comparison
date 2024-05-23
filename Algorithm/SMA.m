function Data_o = SMA(Data_i,Data_o)
rand_num=[];                         
t=clock;                               
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
AllFitness=fitness;
z=0.03; % parameter
it=0;   %Number of iterations
bestPositions=zeros(1,Data_i.dim);

while  it <= Data_i.maxIter
    for i=1:Data_i.pop
        Flag4ub=X(i,:)>Data_i.ub;
        Flag4lb=X(i,:)<Data_i.lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+Data_i.ub.*Flag4ub+Data_i.lb.*Flag4lb;
        AllFitness(i)=Data_i.fobj(X(i,:));
        %update the best fitness value and best position
        if AllFitness(i) < Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
            bestPositions=X(i,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = AllFitness(i);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
       it=it+1;
       if (it > Data_i.maxIter)
           break;
       else
           IterCurve(it)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
       end
    end
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(Data_i.pop);
    bestFitness = SmellOrder(1);
    
    S=bestFitness-worstFitness+eps;  
    
    for i=1:Data_i.pop
        for j=1:Data_i.dim
            if i<=(Data_i.pop/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        bestPositions=X(SmellIndex(1),:);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = bestFitness;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
    end

    a = atanh(-(it/Data_i.maxIter)+1);  
    b = 1-it/Data_i.maxIter;
    % Update the Position of search agents
    for i=1:Data_i.pop
        if rand<z     %Eq.(2.7)
            X(i,:) = (Data_i.ub-Data_i.lb)*rand+Data_i.lb;
        else
            p =tanh(abs(AllFitness(i)-Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)));  %Eq.(2.2)
            vb = unifrnd(-a,a,1,Data_i.dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,Data_i.dim);
            for j=1:Data_i.dim
                r = rand();
                A = randi([1,Data_i.pop]);  
                B = randi([1,Data_i.pop]);
                if r<p    %Eq.(2.1)
                    X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
end

%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end        