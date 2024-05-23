function Data_o = GPC(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
CostFunction=Data_i.fobj;       nVar=Data_i.dim;        VarSize=[1 nVar];  nPop=Data_i.pop;
VarMin=Data_i.lb(1);     VarMax=Data_i.ub(1);      MaxIteration=Data_i.maxIter;

pSS= 0.5;    Tetha = 14;
G = 9.8; 
MuMin = 1;           % Minimum Friction 
MuMax = 10;          % Maximum Friction
stone.Position=[]; 
stone.Cost=[];
pop=repmat(stone,nPop,1);
best_worker.Cost=inf;
for i=1:nPop
   pop(i).Position=X(i,:);
   pop(i).Cost=Cost(i);
   if pop(i).Cost<=best_worker.Cost
       best_worker=pop(i);          % as Pharaoh's special agent
   end
end

for it=1:MaxIteration
    newpop=repmat(stone,nPop,1);
    for i=1:nPop
        newpop(i).Cost = inf;
       
        V0= rand(1,1);                          % Initial Velocity                                      
        Mu= MuMin+(MuMax-MuMin)*rand(1,1);      % Friction

        d = (V0^2)/((2*G)*(sind(Tetha)+(Mu*cosd(Tetha))));                  % Stone Destination
        x = (V0^2)/((2*G)*(sind(Tetha)));                                   % Worker Movement
        epsilon=unifrnd(-((VarMax-VarMin)/2),((VarMax-VarMin)/2),VarSize);  % Epsilon
        newsol.Position = (pop(i).Position+d).*(x*epsilon);                 % Position of Stone and Worker
        newsol.Position = (pop(i).Position+d)+(x*epsilon);                  % Note: In some cases or some problems use this instead of the previous line to get better results

        newsol.Position=max(newsol.Position,VarMin);
        newsol.Position=min(newsol.Position,VarMax);
        
        % Substitution
        z=zeros(size(pop(i).Position));
        k0=randi([1 numel(pop(i).Position)]);
        for k=1:numel(pop(i).Position)
            if k==k0 || rand<=pSS
                z(k)=newsol.Position(k);
            else
                z(k)=pop(i).Position(k);
            end
        end
        
        newsol.Position=z;
        newsol.Cost=CostFunction(newsol.Position);
        
        if newsol.Cost <= newpop(i).Cost
           newpop(i) = newsol;
           if newpop(i).Cost<=best_worker.Cost
               best_worker=newpop(i);
           end
        end
    end
    pop=[pop 
         newpop];
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    pop=pop(1:nPop);
    IterCurve(it)=pop(1).Cost;
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>pop(1).Cost
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=pop(1).Cost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=SortOrder(1);
    end
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end