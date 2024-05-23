function Data_o = EO(Data_i,Data_o)
ti=clock;                               
pos = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
t=0; V=1;
a1=2;
a2=1;
GP=0.5;

Ceq1=zeros(1,Data_i.dim);   Ceq1_fit=inf; 
Ceq2=zeros(1,Data_i.dim);   Ceq2_fit=inf; 
Ceq3=zeros(1,Data_i.dim);   Ceq3_fit=inf; 
Ceq4=zeros(1,Data_i.dim);   Ceq4_fit=inf;
fit_old=fitness;
Max_NFE=Data_i.maxIter;
while t<Max_NFE
    for i=1:Data_i.pop
        for j=1:Data_i.dim
            if  pos(i,j)>Data_i.ub(j)
                pos(i,j)=(Data_i.ub(j));
            elseif  pos(i,j)<Data_i.lb(j)
                pos(i,j)=Data_i.lb(j);
            end
        end
        fitness(i)=Data_i.fobj(pos(i,:));
        if fitness(i)<Ceq1_fit 
            Ceq1_fit=fitness(i);  Ceq1=pos(i,:);
        elseif fitness(i)>Ceq1_fit && fitness(i)<Ceq2_fit  
            Ceq2_fit=fitness(i);  Ceq2=pos(i,:);              
        elseif fitness(i)>Ceq1_fit && fitness(i)>Ceq2_fit && fitness(i)<Ceq3_fit
            Ceq3_fit=fitness(i);  Ceq3=pos(i,:);
        elseif fitness(i)>Ceq1_fit && fitness(i)>Ceq2_fit && fitness(i)>Ceq3_fit && fitness(i)<Ceq4_fit
            Ceq4_fit=fitness(i);  Ceq4=pos(i,:);
        end

        t=t+1;  
        if t>Max_NFE
            break;
        else
            IterCurve(t)=Ceq1_fit;
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Ceq1_fit
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Ceq1_fit;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end    
        end
    end
    %---------------- Memory saving-------------------   
    if t==Data_i.pop
        fit_old=fitness; 
        C_old=pos;
    end
    
    for i=1:Data_i.pop
        if fit_old(i)<fitness(i)                %Local optimal fitness value
            fitness(i)=fit_old(i); 
            pos(i,:)=C_old(i,:);
        end
    end

    C_old=pos;  fit_old=fitness;
    Ceq_ave=(Ceq1+Ceq2+Ceq3+Ceq4)/4;                              % averaged candidate 
    C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];                     % Equilibrium pool
    t1=(1-t/Max_NFE)^(a2*t/Max_NFE);

    for i=1:Data_i.pop
       lambda=rand(1,Data_i.dim);                                % lambda in Eq(11)
       r=rand(1,Data_i.dim);                                     % r in Eq(11)  
       Ceq=C_pool(randi(size(C_pool,1)),:);               % random selection of one candidate from the pool
       F=a1*sign(r-0.5).*(exp(-lambda.*t1)-1);             % Eq(11)
       r1=rand(); r2=rand();                                % r1 and r2 in Eq(15)
       GCP=0.5*r1*ones(1,Data_i.dim)*(r2>=GP);                   % Eq(15)
       G0=GCP.*(Ceq-lambda.*pos(i,:));                      % Eq(14)
       G=G0.*F;                                           % Eq(13)
       pos(i,:)=Ceq+(pos(i,:)-Ceq).*F+(G./lambda*V).*(1-F);   % Eq(16)  
       if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Data_i.fobj(pos(i,:))
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Data_i.fobj(pos(i,:));
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end 
    end
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end