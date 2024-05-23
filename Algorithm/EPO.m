function Data_o = EPO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
x = Data_i.X;                        
y=Data_i.F_value;              

[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);

ts = Data_i.maxIter;        n_dim = Data_i.dim;     k_particle = Data_i.pop;   ub=Data_i.ub;   lb=Data_i.lb;
l_scale = rand*Data_i.ub(1);    res = 0.05; 
eta = (res / l_scale)^(1 / ts);    fobj=Data_i.fobj;  

eta_max = eta;
x_center = Data_i.X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
xmin = Data_i.X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
ymin=y (Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter));
for t = 1 : ts
    x_delta = l_scale*rands(k_particle, n_dim);
    x = x + x_delta;
    for k = 1 : k_particle
        x(k,:)=BoundaryCheck(x(k,:),Data_i.ub,Data_i.lb,Data_i.dim);
        y(k) = fobj(x(k, :));
    end
    [y0, m] = min(y);
    if y0 ~= inf
        x_center = x(m, :);
        y_center = y(m);
    end
    l_scale = l_scale*eta;
    eta = eta_max - t .* ( (eta_max - 0.8) / ts);
    if y_center < ymin
        xmin = x_center;
        ymin = y_center;  
    end

    if y0 < Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        x_center = x(m, :);
        y_center = fobj(x_center);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=y0;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=m;
    end
    IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end