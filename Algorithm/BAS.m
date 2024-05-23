function Data_o = BAS(Data_i,Data_o)
rand_num=[];                        
ti=clock;                            
x = Data_i.X;                        
fitness=Data_i.F_value;                 

[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);

d0=0.001;   d1=3;   d=d1;   eta_d=0.95; l0=0.0;   l1=0.0;   l=l1;   eta_l=0.95;     step=0.8;
n=Data_i.maxIter;       k=Data_i.dim;      fun=Data_i.fobj;   eta_step=0.95;
x0=x;           xbest=x(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
fbest=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
fbest_store=fbest;      

for i=1:n
    for j=1:Data_i.pop
        dir=rands(k,1);
        dir=dir/(eps+norm(dir));

        xleft(j,:)=x(j,:)+(dir*d)';
        fleft(j)=fun(xleft(j,:));
    
        xright(j,:)=x(j,:)-(dir*d)';
        fright(j)=fun(xright(j,:));

        w=l*rands(k,1);
        x(j,:)=x(j,:)-(step*dir*sign(fleft(j)-fright(j))+w)';
        f(j)=fun(x(j,:));
        if f(j)<fbest
            xbest=x(j,:);
            fbest=f(j);
            if fbest<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fbest;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=j;
            end    
        end
    end
    d=d*eta_d+d0;
    l=l*eta_l+l0;
    step=step*eta_step;
    IterCurve(i)=fbest;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=i;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end