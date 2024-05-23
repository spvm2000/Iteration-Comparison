function Data_o = AEFA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                            
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);

FCheck=1; Rpower=1;         
tag=1;                      
BestValues=[];MeanValues=[];
Rnorm=2;     lb=Data_i.lb; ub=Data_i.ub; D=Data_i.dim;
N=Data_i.pop;
X=Data_i.X;
V=zeros(N,D);
max_it=Data_i.maxIter;
N=Data_i.pop;

for iteration=1:max_it
    for i=1:N
        fitness(i)=Data_i.fobj(X(i,:));
    end

    if tag==1
        [best, best_X]=min(fitness); %minimization.   
    else
        [best, best_X]=max(fitness); %maximization.
    end

    if iteration==1
        Fbest=best;Lbest=X(best_X,:);
    end

    if tag==1
        if best<Fbest  %minimization.
            Fbest=best;
            Lbest=X(best_X,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>best
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=best;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=best_X;
            end 
        end    
    else
        if best>Fbest  %maximization
            Fbest=best;
            Lbest=X(best_X,:);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)<best
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=best;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=best_X;
            end
        end    
    end
    BestValues=[BestValues Fbest];
    MeanValues=[MeanValues mean(fitness)];
    Fmax=max(fitness); Fmin=min(fitness); Fmean=mean(fitness);
    if Fmax==Fmin
        M=ones(N,1);
        Q=ones(N,1);
    else
        if tag==1 %for minimization
            best=Fmin;worst=Fmax; 
        else
            best=Fmax;worst=Fmin; 
        end
        Q=exp((fitness-worst)./(best-worst));
    end
    Q=Q./sum(Q);
    fper=3;
    if FCheck==1
        cbest=fper+(1-iteration/max_it)*(100-fper); 
        cbest=round(N*cbest/100);
    else
        cbest=N; 
    end
    [Qs s]=sort(Q,'descend');
    for i=1:N
        E(i,:)=zeros(1,D);
        for ii=1:cbest
            j=s(ii);
            if j~=i
                R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
                 for k=1:D 
                     E(i,k)=E(i,k)+ rand*(Q(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
                 end
             end
        end    
    end
    alfa=30;K0=500;
    K=K0*exp(-alfa*iteration/max_it);
    %---------------------------------------------------------------------------------- 
    %%%Calculation of accelaration.
    a=E*K;

    %Charge movement
    %----------------------------------------------------------------------------------
    V=rand(N,D).*V+a;
    X=X+V;
    X=max(X,lb);X=min(X,ub);   % Check the bounds of the variables
    IterCurve(iteration)=Fbest;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iteration;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end