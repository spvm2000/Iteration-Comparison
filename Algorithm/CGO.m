function Data_o = CGO(Data_i,Data_o)
rand_num=[];                           
ti=clock;                              
pop.Position = Data_i.X;                        
pop.Cost=Data_i.F_value;                   
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ObjFuncName =Data_i.fobj;   Var_Number = Data_i.dim;    LB =Data_i.lb;  UB =Data_i.ub;
MaxIter =Data_i.maxIter;   Seed_Number =Data_i.pop;  
Seed=Data_i.X;
Fun_eval=Data_i.F_value';
%% Search Process of the CGO
for Iter=1:MaxIter
    for i=1:Seed_Number
        % Update the best Seed
        [~,idbest]=min(Fun_eval);
        BestSeed=Seed(idbest,:);
        %% Generate New Solutions
        % Random Numbers
        I=randi([1,2],1,12); % Beta and Gamma
        Ir=randi([0,1],1,5);
        % Random Groups
        RandGroupNumber=randperm(Seed_Number,1);
        RandGroup=randperm(Seed_Number,RandGroupNumber);
        % Mean of Random Group
        MeanGroup=mean(Seed(RandGroup,:)).*(length(RandGroup)~=1)...
            +Seed(RandGroup(1,1),:)*(length(RandGroup)==1);   
        % New Seeds
        Alfa(1,:)=rand(1,Var_Number);
        Alfa(2,:)= 2*rand(1,Var_Number)-1;
        Alfa(3,:)= (Ir(1)*rand(1,Var_Number)+1);
        Alfa(4,:)= (Ir(2)*rand(1,Var_Number)+(~Ir(2)));   
        ii=randi([1,4],1,3);
        SelectedAlfa=Alfa(ii,:);

        NewSeed(1,:)=Seed(i,:)+SelectedAlfa(1,:).*(I(1)*BestSeed-I(2)*MeanGroup);
        NewSeed(2,:)=BestSeed+SelectedAlfa(2,:).*(I(3)*MeanGroup-I(4)*Seed(i,:));
        NewSeed(3,:)=MeanGroup+SelectedAlfa(3,:).*(I(5)*BestSeed-I(6)*Seed(i,:));
        NewSeed(4,:)=unifrnd(LB,UB);
        for j=1:4
            % Checking/Updating the boundary limits for Seeds
            NewSeed(j,:)=BoundaryCheck(NewSeed(j,:),UB,LB,Data_i.dim);
            % Evaluating New Solutions
            Fun_evalNew(j,:)=feval(ObjFuncName, NewSeed(j,:));
        end
        Seed=[Seed; NewSeed];       
        Fun_eval=[Fun_eval; Fun_evalNew];
    end
    [Fun_eval, SortOrder]=sort(Fun_eval);
    Seed=Seed(SortOrder,:);
    [BestFitness,idbest]=min(Fun_eval);   
    BestSeed=Seed(idbest,:);
    Seed=Seed(1:Seed_Number,:);
    Fun_eval=Fun_eval(1:Seed_Number,:);
    if BestFitness<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestFitness;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=idbest(1);
    end    
    IterCurve(Iter)=BestFitness;
%     disp(['Iteration ' num2str(Iter) ': Best Cost = ' num2str(IterCurve(Iter))]);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end