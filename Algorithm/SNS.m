function Data_o = SNS(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
x = Data_i.X;                        
f=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
nDim=Data_i.dim;   LB=Data_i.lb;        UB=Data_i.ub;    GloMin=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
nUser=Data_i.pop;       MaxIter=Data_i.maxIter;     Cost=Data_i.fobj;
for Iter = 1:MaxIter
    for i = 1:nUser
        % Select a random mood
        Mood = randi(4);
        
        % Select random user
        Id = [1:i-1 i+1:nUser];
        j = Id(randi(nUser-1));
        Id(Id == j) = [];
        
        % Follow the procedure of the selected mood
        if Mood == 1
            r = x(j,:) - x(i,:);
            R = rand(1,nDim).*r;
            nx = x(j,:) + (1-2*rand(1,nDim)).*R;
        elseif Mood == 2
            k = Id(randi(nUser-2));
            D = sign(f(i) - f(j))*(x(j,:) - x(i,:));
            nx = x(k,:) + rand(1,nDim).*D;
        elseif Mood == 4
            Group = randperm(nUser,randi(nUser));
            M = mean(x(Group,:), 1);
            nx = x(i,:) + rand(1,nDim).*(M - randi(2)*x(i,:));
        else
            d = randi(nDim);
            n = LB(d) + rand*(UB(d) - LB(d));
            nx = x(i,:);
            t = rand;
            nx(d) = t*n + (1-t)*x(j,d);
        end

        % Clamp the new solution Eq. (6)
        nx = min(max(nx, LB), UB);
        
        % Evaluation of new generated solution(nx)
        nf = Cost(nx);
        
        % Greedy selection
        if nf < f(i)
            f(i) = nf;
            x(i, :) = nx;
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>nf
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=nf;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
    end
    
    % Determine the best user
    [fBest, BestIndex] = min(f);
    xBest = x(BestIndex,:);
    IterCurve(Iter)=fBest;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end