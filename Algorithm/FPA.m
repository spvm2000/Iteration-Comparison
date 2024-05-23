function Data_o = FPA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
Sol = Data_i.X;                        
Fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
para=[Data_i.pop 0.8];
n=Data_i.pop;           % Population size, typically 20 to 40
p=para(2);           % The probabibility of switch
% Iteration parameters
N_iter=Data_i.maxIter;         % Total number of iterations
d=Data_i.dim;   Lb=Data_i.lb;          Ub=Data_i.ub;        Fun=Data_i.fobj;
[fmin,I]=min(Fitness);
best=Sol(I,:);
S=Sol;

for t=1:N_iter
    for i=1:n
        if rand>p
            L=Levy(d);
            dS=L.*(Sol(i,:)-best);      % Caclulate the step increments
            S(i,:)=Sol(i,:)+dS;
            S(i,:)=BoundaryCheck(S(i,:),Ub,Lb,Data_i.dim);
        else
            epsilon=rand;
            JK=randperm(n);
            S(i,:)=S(i,:)+epsilon*(Sol(JK(1),:)-Sol(JK(2),:));
            S(i,:)=BoundaryCheck(S(i,:),Ub,Lb,Data_i.dim);
        end
        Fnew=Fun(S(i,:));
        if (Fnew<=Fitness(i))
            Sol(i,:)=S(i,:);
            Fitness(i)=Fnew;
        end
    end
    [value,index]=min(Fitness);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>value
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=value;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    end
    IterCurve(t)=value;
end    
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function L=Levy(d)
% For details of the Levy flights, see Chapter 11 of the following book:
% Xin-She Yang, Nature-Inspired Optimization Algorithms, Elsevier, (2014).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% Mantegna's algorithm for Levy random numbers
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;
end
