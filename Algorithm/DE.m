function Data_o = DE(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Cost);
IterCurve=zeros(1,Data_i.maxIter);
[Best_Cost,ind] = min(Cost);
Best_X = X(ind,:);

F=0.5;         
Cr=0.5;

it=0;
while it<Data_i.maxIter
    for i=1:Data_i.pop
        Rand_ind=randperm(Data_i.pop);
        Rand_ind(Rand_ind==i)=[];
        
        a = Rand_ind(1);
        b = Rand_ind(2);
        c = Rand_ind(3);

        % Mutation
        y = X(a,:)+F.*(X(b,:)-X(c,:));

        % Crossover
        z=zeros(1,Data_i.dim);
        j0=randi([1 Data_i.dim]);

        for j=1:Data_i.dim
            if j==j0 || rand<=Cr
                z(j)=y(j);
            else
                z(j)=X(i,j);
            end
        end

        New_X = min(max(z,Data_i.lb),Data_i.ub);

        New_Cost=Data_i.fobj(New_X);
        it=it+1;
        if it>Data_i.maxIter
            break;
        end

        if New_Cost<Cost(i)
            X(i,:) = New_X;
            Cost(i) = New_Cost;            
            if Cost(i)<Best_Cost
               Best_X = X(i,:);
               Best_Cost = Cost(i);
               Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Cost;
               Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
        IterCurve(it)=Best_Cost;
    end
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end