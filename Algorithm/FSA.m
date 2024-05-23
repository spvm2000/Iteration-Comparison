function Data_o = FSA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
S = Data_i.X;                        
Fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
iteration=Data_i.maxIter;             d=Data_i.dim;     Lb=Data_i.lb;   Ub=Data_i.ub;       n=Data_i.pop;
[fmin,I]=min(Fitness);      Fun=Data_i.fobj;
best=S(I,:);
Lbe=S;
Lbest=Fitness;
iter=0; 
for t=1:iteration
    for i=1:n
        S(i,:)=S(i,:)+(-S(i,:)+best)*rand+(-S(i,:)+Lbe(i,:))*rand;
        BoundaryCheck(S(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        Fnew(i)=Fun(S(i,:));
        if (Fnew(i)<=Lbest(i)) 
            Lbe(i,:)=S(i,:);
            Lbest(i)=Fnew(i);
        end
        if Fnew(i)<=fmin
            best=S(i,:);
            fmin=Fnew(i);
        end
    end
    for i=1:n
        Si(i,:)=best+(best-Lbe(i,:)).*rand;
        Fitness(i)=Fun(Si(i,:));
        if (Fitness(i)<=Fnew(i)) 
            S(i,:)=Si(i,:);
            Lbe(i,:)=Si(i);
        end
    end
    [fmini,I]=min(Fitness);
    besti=Si(I,:);
    if fmini<=fmin
        best=besti;
        fmin=fmini;
        if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fmini
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fmini;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=I;
        end
    end
    iter=iter+1;
    IterCurve(iter)=fmin;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iter;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end