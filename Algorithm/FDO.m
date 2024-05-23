function Data_o = FDO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
dimensions=Data_i.dim;      fitness=Data_i.fobj;        upper_bound=Data_i.ub;      lower_bound=Data_i.lb;
scout_bee_number=Data_i.pop;

weightFactor = 0.0;
for i=1: scout_bee_number
    scouts(i).xs = X(i,:);
    scouts(i).pace = Cost(i);
end
for iterate = 1 : Data_i.maxIter
    for s=1:scout_bee_number
        current_bee = scouts(s);
        best_bee = getBestScoutBee(scouts,Data_i);
        if fitness(best_bee.xs) ~= 0
            bbf = fitness(best_bee.xs);
            cbf = fitness(current_bee.xs);
            if bbf < (0.05 * cbf)
                fitness_weight = 0.2;
            else
                fitness_weight = fitness(best_bee.xs)/fitness( current_bee.xs) - weightFactor;
            end
        end
        for d=1:dimensions
            x = current_bee.xs(d);
            pace = 0.0;
            random = Levy(1);
            distance_from_best_bee = best_bee.xs(d)- x;
            if fitness_weight == 1
                pace = x * random;
            elseif fitness_weight ==0
                pace = distance_from_best_bee * random;
            else
                pace = (distance_from_best_bee * fitness_weight);
                if random < 0
                    pace = pace * -1;
                end
            end
            x = x +pace;
            tempBee.xs(d) = x;
        end
        tempBee.xs=BoundaryCheck(tempBee.xs,Data_i.ub,Data_i.lb,Data_i.dim);
        tempBee.pace = Data_i.fobj(tempBee.xs);
        if fitness(tempBee.xs) < fitness( current_bee.xs)           %更新局部最优
            scouts(s) = tempBee;
        elseif  size(current_bee.pace, 2) > 1
            for m=1 : dimensions
                x = current_bee.xs(m);
                distance_from_best_bee = best_bee.xs(m)- x;
                x = current_bee.xs(m) + (distance_from_best_bee*fitness_weight)+current_bee.pace(m);
                 x = BoundaryCheck(x.Data_i,ub,Data_i.lb,Data_i.dim);
                tempBee.xs(m)= x;
            end
            if fitness(tempBee.xs) < fitness( current_bee.xs)
                current_bee = tempBee;
            else
                for k=1 : dimensions
                    x = current_bee.xs(k);
                    random = Levy(1);
                    x= x + x*random;
                     x = BoundaryCheck(x.Data_i,ub,Data_i.lb,Data_i.dim);
                    tempBee.xs(k)= x;
                end
            end
        
            if fitness(tempBee.xs) < fitness( current_bee.xs)
                current_bee = tempBee;
            end
        end
        if scouts(i).pace<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=scouts(i).pace;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end    
    end
    IterCurve(iterate)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
end    

Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iterate;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function max_bee = getBestScoutBee(scouts,Data_i)
    max_bee = scouts(1);
    for n=2: Data_i.pop
        if Data_i.fobj(max_bee.xs) >  Data_i.fobj(scouts(n).xs)
            max_bee =  scouts(n);
        end
    end
end

function o=Levy(d)
%% levy函数
    beta=1.5;
    %使用标量和向量计算 gamma 函数
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;
    v=randn(1,d);
    step=u./abs(v).^(1/beta);
    o=step;
end