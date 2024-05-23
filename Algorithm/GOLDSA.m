function Data_o = GOLDSA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
X=Data_i.X;
GoldSA_value=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
GoldSA_position=X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);

Destination_values=Data_i.F_value;

Ub=Data_i.ub;
Lb=Data_i.lb;
a=-pi;                           
b=pi;

gold=0.618;      % golden proportion coefficient, around 0.618
x1=a+(1-gold)*(b-a);          
x2=a+gold*(b-a);
t=1;

while t<Data_i.maxIter
    for i=1:size(X,1) % in i-th solution
        r=rand;
        r1=(2*pi)*r;
        r2=r*pi; 
        for j=1:size(X,2) % in j-th dimension
            X(i,j)= X(i,j)*abs(sin(r1)) - r2*sin(r1)*abs(x1*GoldSA_position(j)-x2*X(i,j));
        end
    end

    for i=1:size(X,1)
        % Check if solutions go outside the search spaceand bring them back
        Boundary_Ub=X(i,:)>Ub;
        Boundary_Lb=X(i,:)<Lb;
        X(i,:)=(X(i,:).*(~(Boundary_Ub+Boundary_Lb)))+Ub.*Boundary_Ub+Lb.*Boundary_Lb;
        % Calculate the objective values
        Destination_values(1,i)=Data_i.fobj(X(i,:));
        
        % Update the destination if there is a better solution
        if Destination_values(1,i)<GoldSA_value
            GoldSA_position=X(i,:);
            GoldSA_value=Destination_values(1,i);
            b=x2;
            x2=x1;
            x1=a+(1-gold)*(b-a);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GoldSA_value;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        else
            a=x1;
            x1=x2;
            x2=a+gold*(b-a);
        end
        
        if x1==x2
            a=-pi*rand; 
            b=pi*rand;
            x1=a+(1-gold)*(b-a); 
            x2=a+gold*(b-a);  
        end 
    end
    t=t+1;
    IterCurve(t)=GoldSA_value;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end