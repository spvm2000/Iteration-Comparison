function Data_o = PSO(Data_i,Data_o)
    c_t=clock;                             
    cnt=1;                                 
    c1 = 2.0;
    c2 = 2.0;
    vmax = 2.*ones(1,Data_i.dim);
    vmin = -2.*ones(1,Data_i.dim);
    V = initialization(Data_i.pop,vmax,vmin,Data_i.dim);
    X = Data_i.X;
    fitness = Data_i.F_value;
    pBest = X;
    pBestFitness = fitness;
    [~,index] = min(fitness);
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
    gBestFitness = fitness(index);
    gBest = X(index,:);

    Xnew = X;            
    fitnessNew = fitness;
    IterCurve = zeros(1,Data_i.maxIter);
    for t = 1:Data_i.maxIter 
       for i = 1:Data_i.pop
          r1 = rand(1,Data_i.dim);
          r2 = rand(1,Data_i.dim);
          V(i,:) = V(i,:) + c1.*r1.*(pBest(i,:) - X(i,:)) + c2.*r2.*(gBest - X(i,:));
          V(i,:) = BoundaryCheck(V(i,:),vmax,vmin,Data_i.dim);
          Xnew(i,:) = X(i,:) + V(i,:);
          Xnew(i,:) = BoundaryCheck(Xnew(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
          fitnessNew(i) = Data_i.fobj(Xnew(i,:));
          if fitnessNew(i) < pBestFitness(i)
              pBest(i,:) = Xnew(i,:);
              pBestFitness(i) = fitnessNew(i);    
          end
          if fitnessNew(i)<gBestFitness
              gBestFitness = fitnessNew(i);
              gBest = Xnew(i,:);
              Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
          end   
       end
       X = Xnew;
       fitness = fitnessNew;
       Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = gBestFitness;
       IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
    end
    Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,c_t);              
    Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                            
    Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end