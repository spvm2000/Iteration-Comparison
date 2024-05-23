function Data_o = SPBO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Objective_values=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
sol=X;      ans=Objective_values;       student=Data_i.pop;         variable=Data_i.dim;  maxi=Data_i.ub(1);   mini=Data_i.lb(1);      
Best_fitness = min(ans);  fobj=Data_i.fobj;
Best_student=sol(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);

for t=1:Data_i.maxIter
    for do=1:1:variable
        sum=zeros(1,variable);
        for gw=1:1:variable
            for fi=1:1:student
                sum(1,gw)=sum(1,gw)+sol(fi,gw);
            end
            mean(1,gw)=sum(1,gw)/student;
        end
        par=sol;
        par1=sol;
        check=rand(student,1);
        mid=rand(student,1);
        for dw=1:student
            if Best_fitness==ans(1,dw)
                jg=ans(randperm(numel(ans),1));
                for oi=1:1:student
                     if jg==ans(1,oi)
                         lk=oi;
                     end
                end
                par1(dw,do)=par(dw,do)+(((-1)^(round(1+rand)))*rand*(par(dw,do)-par(lk,do)));
            elseif check(dw,1)<mid(dw,1)
                rta=rand;
                if rta>rand
                    par1(dw,do)=Best_student(1,do)+(rand*(Best_student(1,do)-par(dw,do)));      % Equation (2a)
                else
                    par1(dw,do)=par(dw,do)+(rand*(Best_student(1,do)-par(dw,do)))+((rand*(par(dw,do)-mean(1,do))));         % Equation (2b)
                end
            else
                an=rand;
                if rand>an
                    par1(dw,do)=par(dw,do)+(rand*(mean(1,do)-par(dw,do)));      % Equation (3)
                else
                    par1(dw,do)=mini+(rand*(maxi-mini));
                end    
            end    
        end
        for z=1:1:student
            if par1(z,do)>maxi
                 par1(z,do)=maxi;
            elseif par1(z,do)<mini
                 par1(z,do)=mini;
            end
        end
        X=par1;
        for i=1:1:size(X,1)
            % Calculate the objective values
            Objective_values(i,1)=fobj(X(i,:));
        end
        fun1=Objective_values;
        for vt=1:1:student
            if ans(1,vt)>fun1(1,vt)
                ans(1,vt)=fun1(1,vt);
                sol(vt,:)=par1(vt,:);
            end
        end
        [Best_fitness1,index]=min(ans);
        Best_student1=sol(index,:);
        if Best_fitness>Best_fitness1
            Best_fitness=Best_fitness1;
            Best_student=Best_student1;
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Best_fitness1
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_fitness1;
            end
        end
    end 
    IterCurve(t)=Best_fitness;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end