function Data_o = EM(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
d=Data_i.dim;       options.lk=Data_i.lb;   options.uk=Data_i.ub;   options.m=Data_i.pop;
options.m0=options.m;   options.MAXITER=Data_i.maxIter;     options.n=length(options.uk);
options.IndexLB=0;      options.LSITER=5;       options.delta=0.1;
options.ObjFunction=Data_i.fobj;    options.Display_Flag=1;     options.run_parallel_index=0;
options.run=5;
nEval=0;
ObjFunction=Data_i.fobj;    n=options.n;    uk=options.uk;  lk=options.lk;  m=options.m;    m0=options.m0;
MAXITER=options.MAXITER;    LSITER=options.LSITER;  delta=options.delta;  

[x,bestX,ibest,ObjFunctionValue,nEval]=EMInitialize(options,nEval,Data_i);
for it=1:MAXITER
    [x,bestX,ibest,ObjFunctionValue,nEval]=EMLocal(x,bestX,ibest,options,ObjFunctionValue,nEval);
    F=EMCalcF(x,options,ObjFunctionValue);
    [x, bestX, ibest,ObjFunctionValue,nEval]=EMMove(F,x,ibest,options,nEval);
    %--------------------------------------------------------------------------------
    [value,index]=min(ObjFunctionValue);
    IterCurve(it)=value(1);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>value(1)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=value(1);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
    end    
    %--------------------------------------------------------------------------------    
end
bestFitness=ObjFunctionValue(ibest);
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [x,xbest,ibest,ObjFunctionValue,nEval]=EMInitialize(options,nEval,Data_i)
    ObjFunction=options.ObjFunction; % the name of the objective function.
    n=options.n;    % n: dimension of the problem.
    uk=options.uk;  % up: upper bound in the kth dimension.
    lk=options.lk;  % lp: lower bound in the kth dimension.
    m=options.m;    % m: number of sample points
    m0=options.m0;   % m: number of initial sample points
    
    for i = 1 : m0
        x(i,:) = Data_i.X(i,:);
        ObjFunctionValue(i)=Data_i.fobj(x(i,:));
    end
    nEval=nEval+size(x,1);
    [index1,index2]=sort(ObjFunctionValue);
    x=x(index2(1:m),:);
    xbest=x(1,:);
    ibest=1;
    ObjFunctionValue=ObjFunctionValue(index2(1:m));
end

function [x,xbest,ibest,ObjFunctionValue,nEval]=EMLocal(x,xbest,ibest,options,ObjFunctionValue,nEval)
    ObjFunction=options.ObjFunction; % the name of the objective function
    n=options.n;    % dimension of the problem.
    uk=options.uk;   % upper bound in the kth dimension.
    lk=options.lk;  % lower bound in the kth dimension.
    m=options.m; % m: number of sample points
    LSITER=options.LSITER; % LSITER: maximum number of local search iterations
    delta=options.delta; % local search parameter, ? ? [0, 1]
    IndexLB=options.IndexLB; % The index for local search on best sample or all oints if IndexLB=1 local search on best only
    if IndexLB==0
        counter =1;
        Length= delta*((uk - lk));
        for i = 1 : m
            for k = 1: n
                landa1=rand;
                while counter <LSITER
                    y=x(i,:);
                    landa2= rand;
                    if landa1 > 0.5
                        y(:,k)=y(:,k)+landa2*Length(k);
                    else
                        y(:,k)=y(:,k)-landa2*Length(k);
                    end
                    ObjFunctiony=feval(ObjFunction,y);
                    nEval=nEval+size(y,1);
                    if ObjFunctiony < feval(ObjFunction,x(i,:))
                        x(i,:)=y;
                        ObjFunctionValue(i)=ObjFunctiony;
                        counter=LSITER-1;
                    end
                    nEval=nEval+1;
                    counter=counter + 1;
                end
            end
            
        end
        [xbest,ibest]=argmin(x,ObjFunctionValue,options);
    else
        counter =1;
        Length= delta*((uk - lk));
        for k = 1: n
            landa1= rand;
            while counter <LSITER
                y=xbest;
                landa2= rand;
                if landa1 > 0.5
                    y(k)=y(k)+landa2*Length(k);
                else
                    y(k)=y(k)-landa2*Length(k);
                end
                ObjFunctiony=feval(ObjFunction,y);
                nEval=nEval+size(y,1);
                if ObjFunctiony < feval(ObjFunction,xbest)
                    xbest=y;
                    ObjFunctionValue(ibest)=ObjFunctiony;
                    counter=LSITER-1;
                end
                nEval=nEval+1;
                counter=counter + 1;
            end
        end
        x(ibest,:)=xbest;
    end
end

function F=EMCalcF(x,options,ObjFunctionValue)
    n=options.n;    % dimension of the problem.
    m=options.m;    % m: number of sample points
    [xbest,ibest,xworst,iworst]=argmin(x,ObjFunctionValue,options);
    for i = 1 : m
        q(i)=real(exp(-2*n*(ObjFunctionValue(i)-ObjFunctionValue(ibest))./(sum(ObjFunctionValue(1:m)-ObjFunctionValue(ibest)))));
    end
    F=zeros(m,n);
    for i = 1 : m
        for j = 1 : m
            if i~=j
                if ObjFunctionValue(j) < ObjFunctionValue(i)
                    F(i,:)= F(i,:) + real(x(j,:)- x(i,:))*q(i)*q(j)./(norm(x(j,:)- x(i,:))).^2; % Attraction
                else
                    F(i,:)= F(i,:) - real(x(j,:)- x(i,:))*q(i)*q(j)./(norm(x(j,:)- x(i,:))).^2;% Repulsion
                end
            end
        end
    end
end

function  [x, xbest, ibest,ObjFunctionValue,nEval]=EMMove(F,x,ibest,options,nEval)
    ObjFunction=options.ObjFunction; % the name of the objective function
    n=options.n;     % dimension of the problem.
    uk=options.uk;   % upper bound in the kth dimension.
    lk=options.lk;   % lower bound in the kth dimension.
    m=options.m;     % m: number of sample points
    for i = 1 : m
        if i ~= ibest
            landa = rand;
            F(i,:)= F(i,:)/norm(F(i,:));
            for k = 1 : n
                if F(i,k)> 0
                    x(i,k)=x(i,k) + landa*F(i,k)*(uk(k)- x(i,k));
                else
                    x(i,k)= x(i,k) + landa*F(i,k)*(x(i,k)-lk(k));
                end
            end
        end
    end
    ObjFunctionValue=feval(ObjFunction,x);
    [xbest,ibest]=argmin(x,ObjFunctionValue,options);
    nEval=nEval+size(x,1);
end

function [xb,ib,xw,iw]=argmin(x,f,options)
    [minf,ib]=min(f);
    xb=x(ib,:);
    [maxf,iw]=max(f);
    xw=x(iw,:);
end