function Data_o = RSO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
Positions= Data_i.X;                       
Ffun=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
Max_iterations=Data_i.maxIter;      objective=Data_i.fobj;
Lower_bound=Data_i.lb;       Upper_bound=Data_i.ub;       dimension=Data_i.dim;
l=0;
x = 1;
y = 5;
R = floor((y-x).*rand(1,1) + x);
Score=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
Position=Positions(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
IterCurve=zeros(1,Data_i.maxIter);
while l<Max_iterations
    for i=1:size(Positions,1)  
        Flag4Upper_bound=Positions(i,:)>Upper_bound;
        Flag4Lower_bound=Positions(i,:)<Lower_bound;
        Positions(i,:)=(Positions(i,:).*(~(Flag4Upper_bound+Flag4Lower_bound)))+Upper_bound.*Flag4Upper_bound+Lower_bound.*Flag4Lower_bound;               
        fitness=objective(Positions(i,:));
        if fitness<Score 
            Score=fitness; 
            Position=Positions(i,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fitness
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end    
        end
        
    end
   
    A=R-l*((R)/Max_iterations); 
    
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
            C=2*rand();          
            P_vec=A*Positions(i,j)+abs(C*((Position(j)-Positions(i,j))));                   
            P_final=Position(j)-P_vec;
            Positions(i,j)=P_final;
        end
    end
    l=l+1;    
    IterCurve(l)=Score;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end