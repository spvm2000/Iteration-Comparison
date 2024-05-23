function Data_o = JAYA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
x = Data_i.X;                        
f=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
pop = Data_i.pop;                   % Population size
var = Data_i.dim;                   % Number of design variables
maxGen = Data_i.maxIter;            % Maximum number of iterations
mini = Data_i.lb;                   % Lower Bound of Variables
maxi = Data_i.ub;                   % Upper Bound of Variables
objective = Data_i.fobj;            % Cost Function

[row,var] = size(mini);
fnew = zeros(pop,1);
fopt= zeros(pop,1);
xopt=zeros(1,var);
gen=1;
while(gen <= maxGen)
    [row,col]=size(x);
    [t,tindex]=min(f);
    Best=x(tindex,:);
    [w,windex]=max(f);
    worst=x(windex,:);
    xnew=zeros(row,col);
    for i=1:row
        for j=1:col
            xnew(i,j)=(x(i,j))+rand*(Best(j)-abs(x(i,j))) - (worst(j)-abs(x(i,j)));  % 
        end
    end

    for i=1:row
        xnew(i,:) = max(min(xnew(i,:),maxi),mini);   
        fnew(i,:) = objective(xnew(i,:));
    end
    for i=1:pop
        if(fnew(i)<f(i))                
            x(i,:) = xnew(i,:);
            f(i) = fnew(i);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>f(i)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=f(i);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end    
        end
    end

    fnew = []; xnew = [];
    [IterCurve(gen),ind] = min(f);
    xopt(gen,:)= x(ind,:);
    gen = gen+1;
end    
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=gen;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end