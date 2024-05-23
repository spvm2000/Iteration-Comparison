function Data_o = ICA(Data_i,Data_o)
rand_num=[];                        
ti=clock;                             
X = Data_i.X;                       
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
CostFunction=Data_i.fobj;       nVar=Data_i.dim;        VarSize=[1 nVar];   VarMin=Data_i.lb(1);    
VarMax=Data_i.ub(1);            MaxIt=Data_i.maxIter;       nPop=Data_i.pop;        

nEmp=10;        alpha=1;        beta=1.5;   pRevolution=0.05;   mu=0.1;   zeta=0.2;
global ProblemSettings;
ProblemSettings.CostFunction=CostFunction;
ProblemSettings.nVar=nVar;
ProblemSettings.VarSize=VarSize;
ProblemSettings.VarMin=VarMin;
ProblemSettings.VarMax=VarMax;

global ICASettings;
ICASettings.MaxIt=MaxIt;
ICASettings.nPop=nPop;
ICASettings.nEmp=nEmp;
ICASettings.alpha=alpha;
ICASettings.beta=beta;
ICASettings.pRevolution=pRevolution;
ICASettings.mu=mu;
ICASettings.zeta=zeta;
emp=CreateInitialEmpires(Data_i);
for it=1:Data_i.maxIter
    emp=AssimilateColonies(emp);
    emp=DoRevolution(emp);
    emp=IntraEmpireCompetition(emp);
    emp=UpdateTotalCost(emp);
    emp=InterEmpireCompetition(emp);
    imp=[emp.Imp];
    [~, BestImpIndex]=min([imp.Cost]);
    BestSol=imp(BestImpIndex);
    IterCurve(it)=BestSol.Cost;
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>BestSol.Cost
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestSol.Cost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=BestImpIndex;
    end
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function emp=CreateInitialEmpires(Data_i)
    global ProblemSettings;
    global ICASettings;
    CostFunction=ProblemSettings.CostFunction;
    nVar=ProblemSettings.nVar;
    VarSize=ProblemSettings.VarSize;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    nPop=ICASettings.nPop;
    nEmp=ICASettings.nEmp;
    nCol=nPop-nEmp;
    alpha=ICASettings.alpha;
    empty_country.Position=[];
    empty_country.Cost=[];
    country=repmat(empty_country,nPop,1);
    for i=1:nPop
        country(i).Position=Data_i.X(i,:);
        country(i).Cost=Data_i.F_value(i);
    end
    costs=[country.Cost];
    [~, SortOrder]=sort(costs);
    country=country(SortOrder);
    imp=country(1:nEmp);
    col=country(nEmp+1:end);
    empty_empire.Imp=[];
    empty_empire.Col=repmat(empty_country,0,1);
    empty_empire.nCol=0;
    empty_empire.TotalCost=[]; 
    emp=repmat(empty_empire,nEmp,1);
    % Assign Imperialists
    for k=1:nEmp
        emp(k).Imp=imp(k);
    end
    % Assign Colonies
    P=exp(-alpha*[imp.Cost]/max([imp.Cost]));
    P=P/sum(P);
    for j=1:nCol
        k=RouletteWheelSelection(P);
        emp(k).Col=[emp(k).Col
                    col(j)];
        emp(k).nCol=emp(k).nCol+1;
    end
    emp=UpdateTotalCost(emp);
end

function i=RouletteWheelSelection(P)
    r=rand;
    C=cumsum(P);
    i=find(r<=C,1,'first');
end

function emp=IntraEmpireCompetition(emp)
    nEmp=numel(emp);
    for k=1:nEmp
        for i=1:emp(k).nCol
            if emp(k).Col(i).Cost<emp(k).Imp.Cost
                imp=emp(k).Imp;
                col=emp(k).Col(i);
                emp(k).Imp=col;
                emp(k).Col(i)=imp;
            end
        end
    end
end

function emp=UpdateTotalCost(emp)
    global ICASettings;
    zeta=ICASettings.zeta;
    nEmp=numel(emp);
    for k=1:nEmp
        if emp(k).nCol>0
            emp(k).TotalCost=emp(k).Imp.Cost+zeta*mean([emp(k).Col.Cost]);
        else
            emp(k).TotalCost=emp(k).Imp.Cost;
        end
    end
end

function emp=InterEmpireCompetition(emp)
    if numel(emp)==1
        return;
    end

    global ICASettings;
    alpha=ICASettings.alpha;
    TotalCost=[emp.TotalCost];
    [~, WeakestEmpIndex]=max(TotalCost);
    WeakestEmp=emp(WeakestEmpIndex);
    P=exp(-alpha*TotalCost/max(TotalCost));
    P(WeakestEmpIndex)=0;
    P=P/sum(P);
    if any(isnan(P))
        P(isnan(P))=0;
        if all(P==0)
            P(:)=1;
        end
        P=P/sum(P);
    end
        
    if WeakestEmp.nCol>0
        [~, WeakestColIndex]=max([WeakestEmp.Col.Cost]);
        WeakestCol=WeakestEmp.Col(WeakestColIndex);
        WinnerEmpIndex=RouletteWheelSelection(P);
        WinnerEmp=emp(WinnerEmpIndex);
        WinnerEmp.Col(end+1)=WeakestCol;
        WinnerEmp.nCol=WinnerEmp.nCol+1;
        emp(WinnerEmpIndex)=WinnerEmp;
        WeakestEmp.Col(WeakestColIndex)=[];
        WeakestEmp.nCol=WeakestEmp.nCol-1;
        emp(WeakestEmpIndex)=WeakestEmp;
    end
    
    if WeakestEmp.nCol==0
        WinnerEmpIndex2=RouletteWheelSelection(P);
        WinnerEmp2=emp(WinnerEmpIndex2);
        
        WinnerEmp2.Col(end+1)=WeakestEmp.Imp;
        WinnerEmp2.nCol=WinnerEmp2.nCol+1;
        emp(WinnerEmpIndex2)=WinnerEmp2;
        
        emp(WeakestEmpIndex)=[];
    end
end

function emp=AssimilateColonies(emp)
    global ProblemSettings;
    CostFunction=ProblemSettings.CostFunction;
    VarSize=ProblemSettings.VarSize;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    global ICASettings;
    beta=ICASettings.beta;
    nEmp=numel(emp);
    for k=1:nEmp
        for i=1:emp(k).nCol
            emp(k).Col(i).Position = emp(k).Col(i).Position ...
                + beta*rand(VarSize).*(emp(k).Imp.Position-emp(k).Col(i).Position);
            emp(k).Col(i).Position = max(emp(k).Col(i).Position,VarMin);
            emp(k).Col(i).Position = min(emp(k).Col(i).Position,VarMax);
            emp(k).Col(i).Cost = CostFunction(emp(k).Col(i).Position);
        end
    end
end

function emp=DoRevolution(emp)
    global ProblemSettings;
    CostFunction=ProblemSettings.CostFunction;
    nVar=ProblemSettings.nVar;
    VarSize=ProblemSettings.VarSize;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    global ICASettings;
    pRevolution=ICASettings.pRevolution;
    mu=ICASettings.mu;
    nmu=ceil(mu*nVar);
    sigma=0.1*(VarMax-VarMin);
    nEmp=numel(emp);
    for k=1:nEmp
        NewPos = emp(k).Imp.Position + sigma*randn(VarSize);
        jj=randsample(nVar,nmu)';
        NewImp=emp(k).Imp;
        NewImp.Position(jj)=NewPos(jj);
        NewImp.Cost=CostFunction(NewImp.Position);
        if NewImp.Cost<emp(k).Imp.Cost
            emp(k).Imp = NewImp;
        end
        
        for i=1:emp(k).nCol
            if rand<=pRevolution
                NewPos = emp(k).Col(i).Position + sigma*randn(VarSize);
                jj=randsample(nVar,nmu)';
                emp(k).Col(i).Position(jj) = NewPos(jj);
                emp(k).Col(i).Position = max(emp(k).Col(i).Position,VarMin);
                emp(k).Col(i).Position = min(emp(k).Col(i).Position,VarMax);
                emp(k).Col(i).Cost = CostFunction(emp(k).Col(i).Position);
            end
        end
    end
end
