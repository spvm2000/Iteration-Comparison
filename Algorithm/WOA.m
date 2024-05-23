function Data_o = WOA(Data_i,Data_o)
rand_nun=[];
ti=clock;
X=Data_i.X;
fitness=Data_i.F_value;
X_temp=X;
fitness_temp=fitness;
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
Best_Pos=X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
IterCurve=zeros(1,Data_i.maxIter);
for t=1:Data_i.maxIter
    for i=1:Data_i.pop
        F_UB=X_temp(i,:)>Data_i.ub;
        F_LB=X_temp(i,:)<Data_i.lb;
        X_temp(i,:)=(X_temp(i,:).*(~(F_UB+F_LB)))+Data_i.ub.*F_UB+Data_i.lb.*F_LB; 
        fitness(i)=Data_i.fobj(X_temp(i,:));
        if fitness(i)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness(i);            
            Best_Pos=X_temp(i,:);                         
            Data_o.Best_Pos=i;                         
        end
    end
    
    a=2-t*(2/Data_i.maxIter);
    a2=-1+t*((-1)/Data_i.maxIter);
    for i=1:Data_i.pop
        rand_nun(1,t)=rand();
        r1=rand_nun(1,t);

        rand_nun(2,t)=rand();
        r2=rand_nun(2,t);

        A=2*a*r1-a;
        C=2*r2;

        b=1;
        l=(a2-1)*rand-1;
        
        rand_nun(3,t)=rand();
        p=rand_nun(3,t);
        for j=1:Data_i.dim
            if p<0.5
                if abs(A)>=1
                    rand_leader_index=floor(Data_i.pop*rand()+1);
                    X_rand=X_temp(rand_leader_index,:);
                    D_X_rand=abs(C*X_rand(j)-X_temp(i,j));
                    X_temp(i,j)=X_rand(j)-A*D_X_rand;
                elseif abs(A)<1
                    D_leader=abs(C*Best_Pos(j)-X_temp(i,j));
                    X_temp(i,j)=Best_Pos(j)-A*D_leader;
                end
            elseif p>=0.5
                distance2Leader=abs(Best_Pos(j)-X_temp(i,j));
                X_temp(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Best_Pos(j);
            end    
        end 
        if fitness(i)>fitness_temp(i)
            fitness(i)=fitness_temp(i);
            X(i,:)=X_temp(i,:);
            if fitness(i)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness(i);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                Best_Pos=X_temp(i,:);
            end
        end
    end
      IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                           
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end