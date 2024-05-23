function Data_o = SOA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  

%% Problem Definition
Max_iterations = Data_i.maxIter;      
Search_Agents = Data_i.pop;      
Positions = Data_i.X;                   
fit = Data_i.F_value;     

[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);

[Score,index]=min(fit);
Position=Positions(index,:); 

Convergence=zeros(1,Max_iterations);

l=0;

while l<Max_iterations
    
    Fc=2-l*((2)/Max_iterations); 
    
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); 
            r2=rand(); 
            
            A1=2*Fc*r1-Fc; 
            C1=2*r2; 
            b=1;             
            ll=(Fc-1)*rand()+1;  
       
            D_alphs=Fc*Positions(i,j)+A1*((Position(j)-Positions(i,j)));                   
            X1=D_alphs*exp(b.*ll).*cos(ll.*2*pi)+Position(j);
            Positions(i,j)=X1;
            
        end
    end

    for i=1:Search_Agents 
        Positions=BoundaryCheck(Positions,Data_i.ub,Data_i.lb,Data_i.dim);               
        fitness=Data_i.fobj(Positions(i,:));
        if fitness<Score 
            Score=fitness; 
            Position=Positions(i,:);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end  
    end
    l=l+1; 
    Convergence(l) = Score;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Score;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence;
end


