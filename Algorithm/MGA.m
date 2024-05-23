function Data_o = MGA(Data_i,Data_o)
rand_num=[];                        
ti=clock;                              
Compan.Position = Data_i.X;                       
Compan.Fun_Eval=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ObjFuncName = Data_i.fobj;      Var_Number = Data_i.dim;    VarMin =Data_i.lb;
VarMax =Data_i.ub;          MaxIteration = Data_i.maxIter;  NCompan = Data_i.pop;
Globalbest=0;
[value, Index1]=min(Compan.Fun_Eval);
BestSoFar.Fun_Eval(1) =value;
BestSoFar.Position (1,:)= Compan.Position(Index1,:);

for Iter=2:Data_i.maxIter
    Index=NCompan.*(0:Var_Number-1)+randperm(NCompan,Var_Number);
    CompnNew.Position(1,:)= Compan.Position(Index)+unifrnd(-1,1).*randn(1,Var_Number);
    
    Index=randperm(NCompan,1);
    Index2= randperm(NCompan,Index);
    CMs=randn(Index,1);
    CMs=CMs/sum(CMs);
    CompnNew.Position(2,:)=sum(CMs.*Compan.Position(Index2,:));
    CompnNew.Position = max(CompnNew.Position,VarMin);
    CompnNew.Position = min(CompnNew.Position,VarMax);

    for ind3=1:size(CompnNew.Position,1)
        CompnNew.Fun_Eval(ind3)=feval(ObjFuncName,CompnNew.Position(ind3,:));
    end
    AllCompn.Position=[Compan.Position;CompnNew.Position];
    AllCompn.Fun_Eval=[Compan.Fun_Eval,CompnNew.Fun_Eval];
    [~, Index1]=sort(AllCompn.Fun_Eval);
    Compan.Position=AllCompn.Position(Index1(1:NCompan),:);
    Compan.Fun_Eval=AllCompn.Fun_Eval(Index1(1:NCompan));

    To_UpdateG= Compan.Fun_Eval(1) < BestSoFar.Fun_Eval(Iter-1);
    best = (1-To_UpdateG).*BestSoFar.Position(Iter-1,:)+To_UpdateG.*Compan.Position(1,:);
    BestSoFar.Fun_Eval (Iter)= (1-To_UpdateG).*BestSoFar.Fun_Eval(Iter-1)+To_UpdateG.*Compan.Fun_Eval(1);
    BestSoFar.Position(Iter,:)=best;
    
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>min(CompnNew.Fun_Eval)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=min(CompnNew.Fun_Eval);
        [a,b]=min(CompnNew.Fun_Eval);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=b(1);
    end    
    IterCurve(Iter)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
    Showindex=1;
    Iter=Iter+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end