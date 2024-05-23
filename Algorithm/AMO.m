function Data_o = AMO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
p = Data_i.X;                        
fit=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
time=1;
total_time=1;
low=Data_i.lb;
up=Data_i.ub;
GlobalMin=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
while time<=Data_i.maxIter
    for i=1:Data_i.pop
        fit(i)=Data_i.fobj(p(i,:));
    end

    ind=find(fit==min(fit));
    ind=ind(end);
    GlobalParams=p(ind,:);
    for i=1:Data_i.pop
        p(i,:)=BoundaryCheck(p(i,:),up,Data_i.lb,Data_i.dim);
        f=i*ones(1,Data_i.dim);
        FF=normrnd(0,1);
        for d=1:Data_i.dim
            if i==1
                 lseq=[Data_i.pop-1 Data_i.pop i i+1 i+2];
                 j=randperm(5);
                 f(d)=lseq(j(2));
             elseif i==2
                 lseq=[Data_i.pop i-1 i i+1 i+2];
                 j=randperm(5);
                 f(d)=lseq(j(2));
             elseif i==Data_i.pop-1
                 lseq=[i-2 i-1 i Data_i.pop 1];
                 j=randperm(5);
                 f(d)=lseq(j(2));
             elseif i==Data_i.pop
                 lseq=[i-2 i-1 i 1 2];
                 j=randperm(5);
                 f(d)=lseq(j(2));
             else
                 lseq=[i-2 i-1 i i+1 i+2];
                 j=randperm(5);
                 f(d)=lseq(j(2));
            end
            newV(i,d)=p(i,d)+FF.*(p(f(d),d)-p(i,d));
        end
        newV(i,:)=BoundaryCheck(newV(i,:),up,Data_i.lb,Data_i.dim);
        fit_V(i)=Data_i.fobj(newV(i,:));
        if fit_V(i)<=fit(i)                     
            p(i,:)= newV(i,:);
            fit(i)=fit_V(i);
        end
    end
    [sortVal, sortIndex] = sort(fit);
    for i=1:Data_i.pop     
        pro(sortIndex(i))=(Data_i.pop -i+1)./Data_i.pop;
    end
    ind=find(fit==min(fit));
    ind=ind(end);
    if (fit(ind)<GlobalMin)
        GlobalMin=fit(ind);
        GlobalParams=p(ind,:);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GlobalMin;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=ind;
    end
    [r1,r2,r3,r4,r5] =getindex(Data_i.pop);
    for i=1:Data_i.pop
        for j=1:Data_i.dim
            if rand>pro(i) 
                newVV(i,j)= p(r1(i),j)+rand*(GlobalParams(j)-p(i,j))+ rand*(p(r3(i),j)-p(i,j));
            else
                newVV(i,j)= p(i,j);
            end
        end
        newVVnewVV(i,:)=BoundaryCheck(newVV(i,:),up,Data_i.lb,Data_i.dim);
        fit_VV(i)=Data_i.fobj(newVV(i,:));
        if fit_VV(i)<=fit(i)
            p(i,:)=newVV(i,:);
            fit(i)=fit_VV(i);
        end
    end
    ind=find(fit==min(fit));
    ind=ind(end);
    if (fit(ind)<GlobalMin)
        GlobalMin=fit(ind);
        GlobalParams=p(ind,:);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GlobalMin;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=ind;
    end
    IterCurve(time)=GlobalMin;
    time=time+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=time;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [r1,r2,r3,r4,r5] =getindex(popsize)
r1=zeros(1,popsize);
r2=zeros(1,popsize);
r3=zeros(1,popsize);
r4=zeros(1,popsize);
r5=zeros(1,popsize);
for i=1:popsize
    sequence=1:popsize;
    sequence(i)=[];

    temp=floor(rand*(popsize-1))+1;
    r1(i)=sequence(temp);
    sequence(temp)=[];

    temp=floor(rand*(popsize-2))+1;
    r2(i)=sequence(temp);
    sequence(temp)=[];

    temp=floor(rand*(popsize-3))+1;
    r3(i)=sequence(temp);
    sequence(temp)=[];

    temp=floor(rand*(popsize-4))+1;
    r4(i)=sequence(temp);
    sequence(temp)=[];

    temp=floor(rand*(popsize-5))+1;
    r5(i)=sequence(temp);
end
end