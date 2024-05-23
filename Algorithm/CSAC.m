function Data_o = CSAC(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X_t = Data_i.X;                        
fitness=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ub=Data_i.ub;       lb=Data_i.lb;       fobj=Data_i.fobj;
Max_iter=Data_i.maxIter;          [best,indx]=min(fitness); X_c(1,:)=X_t(indx,:);
l=0;        SearchAgents_no=Data_i.pop;     dim=Data_i.dim;
while l<Max_iter
    a=pi-pi*(l/Max_iter)^2; 
    for i=1:SearchAgents_no
        w=a*rand-a;
        p=1-0.9*(l/Max_iter)^0.5;
        c=0.8;
        for j=1:dim
            if l>(c*Max_iter) % make c>=0.75 for functions F1-F13 and make it c<=0.4 for other functions
                X_t(i,j)=X_c(1,j)+(X_c(1,j)-X_t(i,j))*tan(w*rand);% escape local stagnation
            else
                X_t(i,j)=X_c(1,j)-(X_c(1,j)-X_t(i,j))*tan(w*p);
            end
    
            if X_t(i,j)>ub(j)
                X_t(i,j)=ub(j);
            elseif X_t(i,j)<lb(j)
                X_t(i,j)=lb(j);
            end
        end
        fitness(i)=fobj(X_t(i,:));
        if fitness(i)<best
            best=fitness(i); % Update best
            X_c(1,:)=X_t(i,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>best
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=best;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end    
        end
    end
    l=l+1;    
    IterCurve(l)=best;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end