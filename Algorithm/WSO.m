function Data_o = WSO(Data_i,Data_o)
rand_num=[];                        
ti=clock;                          
Positions = Data_i.X;                        
fitness_old=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
 Max_iter=Data_i.maxIter;            dim=Data_i.dim;     fobj=Data_i.fobj;           pop_size=Data_i.pop;
King=zeros(1,dim);      King_fit=inf;
fitness_old=inf*ones(1,pop_size);            fitness_new=inf*ones(1,pop_size);
l=1;
W1=2*ones(1,pop_size);
Wg=zeros(1,pop_size);     R=0.1;
King_fit=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
King=Positions(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
while l<=Max_iter
    [~,tindex]=sort(fitness_old);
    Co=Positions(tindex(2),:);
    iter =l;
    com=randperm(pop_size);
    for i=1:pop_size
        RR=rand;
        if RR<R
            D_V(i,:)=2*RR*(King-Positions(com(i),:))+1*W1(i)*rand*(Co-Positions(i,:)); 
        else
            D_V(i,:)=2*RR*(Co-King)+1*rand*(W1(i)*King-Positions(i,:)); 
        end
        Positions_new(i,:)=Positions(i,:)+D_V(i,:);
        Positions_new(i,:)=BoundaryCheck(Positions_new(i,:),Data_i.ub,Data_i.lb,Data_i.dim);  
        fitness=fobj(Positions_new(i,:));
        fitness_new(i) = fitness;

        if fitness<King_fit
            King_fit=fitness;
            King=Positions_new(i,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>King_fit
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=King_fit;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end

        if fitness<fitness_old(i)
            Positions(i,:)=Positions_new(i,:);
            fitness_old(i)=fitness;
            Wg(i)=Wg(i)+1;
            W1(i)=1*W1(i)*(1-Wg(i)/Max_iter)^2;
        end   
    end
    IterCurve(l)=King_fit;
    l=l+1;    
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end