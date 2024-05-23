function Data_o = TSAA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
Positions = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
Score=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
SearchAgents=Data_i.pop;        Max_iterations=Data_i.maxIter;      
Position=zeros(1,Data_i.dim);   objective=Data_i.fobj;
t=1;
while t<=Max_iterations
    xmin=1;
    xmax=4;
    xr=xmin+rand()*(xmax-xmin);
    xr=fix(xr);
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
            A1=((rand()+rand())-(2*rand()))/xr;     
            c2=rand();
            if(i==1)
                c3=rand(); 
                if(c3>=0.5)
                    d_pos=abs(Position(j)-c2*Positions(i,j));
                    Positions(i,j)=Position(j)+A1*d_pos;
                else
                    d_pos=abs(Position(j)-c2*Positions(i,j));
                    Positions(i,j)=Position(j)-A1*d_pos;
                end
            else
                c3=rand(); 
                if(c3>=0.5)
                    d_pos=abs(Position(j)-c2*Positions(i,j));
                    Pos(i,j)=Position(j)+A1*d_pos;
                else                          
                    Pos(i,j)=Position(j)-A1*d_pos;
                end
                 Positions(i,j)=(Pos(i,j)+Positions(i-1,j))/2;
            end
        end
        Positions(i,:)=BoundaryCheck(Positions(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
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
    IterCurve(t)=Score;
    t=t+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end