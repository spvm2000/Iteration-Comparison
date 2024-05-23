function Data_o = CSOA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
PopPos = Data_i.X;                        
PopFit=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
global nfe;
alg=1;  crs=1;
nPop=Data_i.pop;    MaxIt=Data_i.maxIter;   LF=2.0;     Pv=rand;    maxrn=1;
Low=Data_i.lb;  Up=Data_i.ub;   Dim=Data_i.dim;     
maxnfe=10000*Dim;
costrn=zeros(1,maxrn);
posrn=zeros(maxrn,Dim);
BestF=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
BestX=Data_i.X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
for it=1:MaxIt
    npos=PopPos;
    for i=1:nPop
         L=randperm(nPop,1);
         M=randperm(nPop,1);
         for j=1:Dim
             c1=-1+2*rand;
             r1=rand;
             npos(L,j)=r1*PopPos(L,j)+(1-r1)*PopPos(M,j)+c1*(PopPos(L,j)-PopPos(M,j));
             c2=-1+2*rand;
             r2=rand;
              npos(M,j)=r2*PopPos(M,j)+(1-r2)*PopPos(L,j)+c2*(PopPos(M,j)-PopPos(L,j));
         end
    end
    for i=1:nPop
        npos(i,:)=BoundaryCheck(npos(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        nfit(i)=Data_i.fobj(npos(i,:));
    end
    for i=1:nPop
       if nfit(i)<PopFit(i)
           PopFit(i)=nfit(i);
           for j=1:Dim
               PopPos(i,j)=npos(i,j);
           end
       end
    end
    for i=1:nPop
       for j=1:Dim
           normpos(i,j)=(PopPos(i,j)-Low)/(Up-Low);
       end
    end
    for i=1:nPop
        r3=rand;
        if r3<Pv
            for j=1:Dim
                r4=rand;
                rpos(i,j)=r4*normpos(i,j)+(1-r4)*normpos(i,j);
            end
        else
             for j=1:Dim
                 rpos(i,j)=normpos(i,j);
             end
        end
    end
    for i=1:nPop
        npos(i,:)=rpos(i,:).*(Up-Low)+Low;
        npos(i,:)=BoundaryCheck(npos(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        nfit(i)=Data_i.fobj(npos(i,:));
    end
    for i=1:nPop
        if nfit(i)<PopFit(i)
            PopFit(i)=nfit(i);
            for j=1:Dim
                PopPos(i,j)=npos(i,j);
            end
        end
    end
    for i=1:nPop
        if PopFit(i)<=BestF
            BestF=PopFit(i);
            BestX=PopPos(i,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>BestF
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestF;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
    end
   IterCurve(it)=BestF;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end