function Data_o = CRYSTAL(Data_i,Data_o)
rand_num=[];                           
ti=clock;                               
Pop = Data_i.X;                        
Cost=Data_i.F_value;                   
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ObjFuncName =Data_i.fobj;  Var_Number =Data_i.dim; 
LB =Data_i.lb;  UB =Data_i.ub;
MaxIteation =Data_i.maxIter;  Cr_Number =Data_i.pop;  Crystal=Data_i.X;
Fun_eval=Data_i.F_value;
% The best Crystal
[BestFitness,idbest]=min(Fun_eval);

Crb=Crystal(idbest,:);
Eval_Number=Cr_Number;
Iter=0;
while Iter<MaxIteation
    for i=1:Cr_Number
        %% Generate New Crystals
        % Main Crystal
        Crmain=Crystal(randperm(Cr_Number,1),:);
        % Random-selected Crystals
        RandNumber=randperm(Cr_Number,1);
        RandSelectCrystal=randperm(Cr_Number,RandNumber);
        % Mean of randomly-selected Crystals
        Fc=mean(Crystal(RandSelectCrystal,:)).*(length(RandSelectCrystal)~=1)...
            +Crystal(RandSelectCrystal(1,1),:)*(length(RandSelectCrystal)==1);   
        % Random numbers (-1,1)
       r=2*rand-1;       r1=2*rand-1;
       r2=2*rand-1;     r3=2*rand-1;        
        % New Crystals
       NewCrystal(1,:)=Crystal(i,:)+r*Crmain;
       NewCrystal(2,:)=Crystal(i,:)+r1*Crmain+r2*Crb;
       NewCrystal(3,:)=Crystal(i,:)+r1*Crmain+r2*Fc;
       NewCrystal(4,:)=Crystal(i,:)+r1*Crmain+r2*Crb+r3*Fc;

        for i2=1:4
            % Checking/Updating the boundary limits for Crystals
            NewCrystal(i2,:)=BoundaryCheck(NewCrystal(i2,:),UB,LB,Data_i.dim);
            % Evaluating New Crystals
            Fun_evalNew(i2)=feval(ObjFuncName, NewCrystal(i2,:));
            % Updating the Crystals
            if Fun_evalNew(i2)<Fun_eval(i)
                Fun_eval(i)=Fun_evalNew(i2);
                Crystal(i,:)=NewCrystal(i2,:);
            end
            % Updation the Number of Function Evalutions
            Eval_Number=Eval_Number+1;
        end
    end % End of One Iteration
    Iter=Iter+1;
   % The best Crystal
    [BestFitness,idbest]=min(Fun_eval);
    Crb=Crystal(idbest,:);
    BestCr=Crb;
    if BestFitness<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestFitness;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=idbest(1);
    end    
    IterCurve(Iter)=BestFitness;
end 
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end