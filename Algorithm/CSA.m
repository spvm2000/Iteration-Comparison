function Data_o = CSA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
tmax = Data_i.maxIter;             
N = Data_i.pop;                 
pd = Data_i.dim;               
x = Data_i.X;                   
ft = Data_i.F_value;     
l=Data_i.lb;
u=Data_i.ub;
ffit=zeros(1,Data_i.maxIter);
AP=0.1; % Awareness probability
fl=2; % Flight length (fl)

mem=x; % Memory initialization
fit_mem=ft; % Fitness of memory positions

[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(ft);

for t=1:tmax

    num=ceil(N*rand(1,N)); % Generation of random candidate crows for following (chasing)
    for i=1:N
        if rand>AP
            xnew(i,:)= x(i,:)+fl*rand*(mem(num(i),:)-x(i,:)); % Generation of a new position for crow i (state 1)
        else
            for j=1:pd
                xnew(i,j)=l(1,j)-(l(1,j)-u(1,j))*rand; % Generation of a new position for crow i (state 2)
            end
        end
    end

    for i=1:N % Update position and memory
        if xnew(i,:)>=l & xnew(i,:)<=u
            x(i,:)=xnew(i,:); % Update position
            tempft=Data_i.fobj(x(i,:));
            if tempft<fit_mem(i)
                mem(i,:)=x(i,:); % Update memory
                fit_mem(i)=tempft;
            end
        end
    end
    [bestfit,index]=min(fit_mem); % Best found value until iteration t
    ffit(t)=bestfit;
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = bestfit;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=ffit;
end



