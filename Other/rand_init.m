function Data_i = rand_init(Data_i)
%% rand初始化
for i = 1:Data_i.pop
       for j = 1:Data_i.dim
           if  Data_i.ub(j)==Data_i.lb(j)
               disp(['ub=',num2str(Data_i.ub(j))]);
               disp(['lb=',num2str(Data_i.lb(j))]);
               error('上下界不能相同!');
           elseif (Data_i.ub(j)>0 && Data_i.lb(j)>=0) || (Data_i.ub(j)>0 && Data_i.lb(j)<=0)
                Data_i.X(i,j) = rand*(Data_i.ub(j)-Data_i.lb(j))+Data_i.lb(j);       
           else
               Data_i.X(i,j) = rand*(Data_i.ub(j))+rand*(Data_i.lb(j))
           end
       end
    end
end