function Data_o = LAPO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
xold = Data_i.X;                        
fold=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ub=Data_i.ub;       lb=Data_i.lb;   ng=Data_i.pop;   tmax=Data_i.maxIter;   xav=mean(xold,1);
nv=Data_i.dim;      fav=Data_i.fobj(xav);       fobj=Data_i.fobj;

for t=1:Data_i.maxIter
    [a,b]=sort(fold);
    if fav<fold(b(1,ng))
        fold(b(1,ng))=(fav);
        xold(b(1,ng),:)=xav;
    end
    
    [a,b]=sort(fold);
    xav=mean(xold,1);
    fav=fobj(xav);
    for r=1:ng
        u=randperm(ng);
        if fav<fold(u(r))
            for i=1:nv
                xnew(r,i)=xold(r,i)+rand*(+xav(1,i)-rand*xold(u(r),i));
                if  xnew(r,i)>ub(1,i)
                    xnew(r,i)=ub(1,i);
                end
                if  xnew(r,i)<lb(1,i)
                    xnew(r,i)=lb(1,i);
                end
                
            end
            %          xnew(r,2)=xold(r,2)+rand*(+xav1(1,2)-rand*xold(u(r),2));
        else
            for i=1:nv
                xnew(r,i)=xold(r,i)-rand*(+xav(1,i)-rand*xold(u(r),i));
                % cheak uper and lower band
                if  xnew(r,i)>ub(1,i)
                    xnew(r,i)=ub(1,i);
                end
                if  xnew(r,i)<lb(1,i)
                    xnew(r,i)=lb(1,i);
                end
            end
            %          xnew(r,2)=xold(r,2)-rand*(+xav1(1,2)-rand*xold(u(r),2));
        end
        %% better result
        fnew(r)=fobj(xnew(r,:));
        if fnew(r)>fold(r)                  
            fnew(r)=fold(r);
            xnew(r,:)=xold(r,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fnew(r)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fnew(r);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=r;
            end    
        end
    end
    %% Phase 2
    fold=fnew;
    xold=xnew;
    [a,b]=sort(fold);
    
    for r=1:ng
        for i=1:nv
            yu=1-(t*1/tmax)*exp(-(t^1)/tmax);
            y=(+yu*(xold(b(1,ng),i)-xold(b(1,1),i)));
            xnew(r,i)=xold(r,i)-rand*y;
            % cheak uper and lower band
            if  xnew(r,i)>ub(1,i)
                xnew(r,i)=ub(1,i);
            end
            if  xnew(r,i)<lb(1,i)
                xnew(r,i)=lb(1,i);
            end
        end
        fnew(r,1)=fobj(xnew(r,:));
        %% better result
        if fnew(r)>fold(r)
            fnew(r)=fold(r);
            xnew(r,:)=xold(r,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fnew(r)
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fnew(r);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=r;
            end
        end
    end
    
    fold=fnew;
    xold=xnew;
    [a,b]=sort(fold);
    IterCurve(t)=a(1);
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end