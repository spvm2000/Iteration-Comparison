function Data_o = GSA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
min_flag=1;                 
Rnorm=2; low=Data_i.lb; up=Data_i.ub; dim=Data_i.dim;
X=Data_i.X;   BestChart=[];   MeanChart=[];  V=zeros(Data_i.pop,Data_i.dim);
ElitistCheck=1; Rpower=1;  max_it=Data_i.maxIter;  F_index=1;
for iteration=1:max_it
    for i=1:Data_i.pop
        X(i,:) = BoundaryCheck(X(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        fitness(i)=Data_i.fobj(X(i,:));
    end

    if min_flag==1
        [best best_X]=min(fitness); %minimization.
    else
        [best best_X]=max(fitness); %maximization.
    end

    if iteration==1
       Fbest=best; Lbest=X(best_X,:);
    end

    if min_flag==1
        if best<Fbest  %minimization.
            Fbest=best; Lbest=X(best_X,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Fbest;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=best_X;
        end
    else 
        if best>Fbest  %maximization
            Fbest=best; Lbest=X(best_X,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Fbest;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=best_X;
        end
    end
    IterCurve(iteration)=best;
    BestChart=[BestChart Fbest];
    MeanChart=[MeanChart mean(fitness)];

    %eq.14-20
    [M]=massCalculation(fitness,min_flag);
    %eq.13.
    G=Gconstant(iteration,max_it);
    %eq.7-10,21.
    a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);
    %eq.11-12
    [X,V]=move(X,a,V);
    
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iteration;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function G=Gconstant(iteration,max_it)
  alfa=20;G0=100;
  G=G0*exp(-alfa*iteration/max_it); %eq. 28.
end

function a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it)
    [N,dim]=size(X);
     final_per=2; 
     if ElitistCheck==1
         kbest=final_per+(1-iteration/max_it)*(100-final_per); %kbest in eq. 21.
         kbest=round(N*kbest/100);
     else
         kbest=N; %eq.9.
     end
        [Ms ds]=sort(M,'descend');
     for i=1:N
         E(i,:)=zeros(1,dim);
         for ii=1:kbest
             j=ds(ii);
             if j~=i
                R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
             for k=1:dim 
                 E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
             end
             end
         end
     end
    a=E.*G;  %note that Mp(i)/Mi(i)=1
end

function [M]=massCalculation(fit,min_flag) 
    Fmax=max(fit); Fmin=min(fit); Fmean=mean(fit); 
    [i N]=size(fit);
    if Fmax==Fmin
       M=ones(N,1);
    else
       if min_flag==1 %for minimization
          best=Fmin;worst=Fmax; %eq.17-18.
       else %for maximization
          best=Fmax;worst=Fmin; %eq.19-20.
       end
       M=(fit-worst)./(best-worst); %eq.15,
    end
    M=M./sum(M); %eq. 16.
end

function [X,V]=move(X,a,V)
    [N,dim]=size(X);
    V=rand(N,dim).*V+a; %eq. 11.
    X=X+V; %eq. 12.
end