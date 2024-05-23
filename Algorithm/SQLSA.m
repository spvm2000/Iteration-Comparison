function Data_o = SQLSA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
x = Data_i.X;                        
fitness=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
Max_iter=Data_i.maxIter;    N=Data_i.pop;       lb=Data_i.lb;   ub=Data_i.ub;     dim=Data_i.dim;

nfs=4;                  
hnt=1;                  
ant=3;                 
noft=46;               
Fmax=1.11;                 
Fmin=0.5;                 %minimum gliding distance
A=rand(N,1);           
r=rand(N,1);            %pulse flying rate for each SSA
Gc=1.9; %Gliding constant
% Initializing arrays
F=zeros(N,1);           % Frequency
v=zeros(N,dim);         % Velocities

for inum=1:length(fitness)
    update_ssa(inum)=randi(3);
end
[fmin,index]=min(fitness);          %find the initial best fitness value,
bestsol=x(index,:);
iter=1;
while iter<=Max_iter
    for ii=1:Data_i.pop
        if update_ssa(ii) == 1
            F(ii)=Fmin+(Fmax-Fmin)*rand;              %randomly chose the gliding distance
            v(ii,:)=v(ii,:)+F(ii)*Gc*(x(ii,:)-bestsol)*1;  %update the velocity for acorn tree squrrels
            x(ii,:)=x(ii,:)+v(ii,:);
        elseif update_ssa(ii) == 2
            F(ii)=Fmin+(Fmax-Fmin)*rand;              %randomly chose the gliding distance
            v(ii,:)=v(ii,:)+F(ii)*Gc*(x(ii,:)-bestsol)*2;  %update the for normal tree squrrels
            x(ii,:)=x(ii,:)+v(ii,:);
        else
            F(ii)=Fmin+(Fmax-Fmin)*rand;              %randomly chose the gliding distance
            v(ii,:)=v(ii,:)+F(ii)*Gc*(x(ii,:)-bestsol)*3;  %update the velocity for hickory tree squrrels
            x(ii,:)=x(ii,:)+v(ii,:); 
        end
        x(ii,:)=BoundaryCheck(x(ii,:),Data_i.ub,Data_i.lb,Data_i.dim);
        if rand>r(ii)
            % The factor 0.001 limits the step sizes of random flyes
            eps=-1+(1-(-1))*rand;
            x(ii,:)=bestsol+eps*mean(A);
        end
        fitnessnew=Data_i.fobj(x(ii,:));
        if fitnessnew<fitness(ii)
            fitness(ii)=fitnessnew;
        end    
        if fitnessnew<=fmin
            bestsol=x(ii,:);
            fmin=fitnessnew;
        end
    end
    IterCurve(iter)=fmin;
    iter=iter+1; 
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iter;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end
