function Data_o = CSAF(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
nest = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
n=Data_i.pop;    Lb=Data_i.lb;      Ub=Data_i.ub;    
fmin=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
bestnest=nest(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
pa=0.25;    it=1;
while it<=Data_i.maxIter
    new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
    [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,Data_i);
    if fnew<fmin
        fmin=fnew;
        bestnest=best;
    end
    IterCurve(it)=fmin;
    it=it+1;
end    

Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,Data_i)
    % Evaluating all new solutions
    for j=1:size(nest,1)
        fnew=Data_i.fobj(newnest(j,:));
        if fnew<=fitness(j)
           fitness(j)=fnew;
           nest(j,:)=newnest(j,:);
        end
    end
    % Find the current best
    [fmin,K]=min(fitness) ;
    best=nest(K,:);
end

function nest=get_cuckoos(nest,best,Lb,Ub)
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    for j=1:size(nest,1)
        s=nest(j,:);
        u=randn(size(s))*sigma;
        v=randn(size(s));
        step=u./abs(v).^(1/beta);    
        stepsize=0.01*step.*(s-best);
        s=s+stepsize.*randn(size(s));
        nest(j,:)=BoundaryCheck(s,Ub,Lb,length(Lb));
    end
end