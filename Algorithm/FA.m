function Data_o = FA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                            
X = Data_i.X;                        
Ffun=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
CostFunction=Data_i.fobj;   nVar=Data_i.dim;        VarSize=[1 nVar];
VarMin=Data_i.lb(1);    VarMax=Data_i.ub(1);        MaxIt=Data_i.maxIter;       nPop=Data_i.pop;
gamma=1;    beta0=2;    alpha=0.2;      alpha_damp=0.98;    delta=0.05*(VarMax-VarMin);   m=2;

if isscalar(VarMin) && isscalar(VarMax)
    dmax = (VarMax-VarMin)*sqrt(nVar);
else
    dmax = norm(VarMax-VarMin);
end
firefly.Position=[];
firefly.Cost=[];
pop=repmat(firefly,nPop,1);

BestSol.Cost=inf;

for i=1:nPop
   pop(i).Position=X(i,:);
   pop(i).Cost=Ffun(i);
   
   if pop(i).Cost<=BestSol.Cost
       BestSol=pop(i);
   end
end

for it=1:MaxIt
    newpop=repmat(firefly,nPop,1);
    for i=1:nPop
         newpop(i).Cost = inf;
         for j=1:nPop
            if pop(j).Cost < pop(i).Cost
                rij=norm(pop(i).Position-pop(j).Position)/dmax;
                beta=beta0*exp(-gamma*rij^m);
                e=delta*unifrnd(-1,+1,VarSize);
                newsol.Position = pop(i).Position ...
                                + beta*rand(VarSize).*(pop(j).Position-pop(i).Position) ...
                                + alpha*e;
                newsol.Position=max(newsol.Position,VarMin);
                newsol.Position=min(newsol.Position,VarMax);
                
                newsol.Cost=CostFunction(newsol.Position);
                 if newsol.Cost <= newpop(i).Cost
                    newpop(i) = newsol;
                    if newpop(i).Cost<=BestSol.Cost
                        BestSol=newpop(i);
                        if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>BestSol.Cost
                            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestSol.Cost;
                            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                        end
                    end
                end
            end    
         end    
    end
    pop=[pop
         newpop];  %#ok
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    % Truncate
    pop=pop(1:nPop);
    % Store Best Cost Ever Found
    IterCurve(it)=BestSol.Cost;
    %Damp Mutation Coefficient
    alpha = alpha*alpha_damp;
end    
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end