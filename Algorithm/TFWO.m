function Data_o = TFWO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                       
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
nWh=5;      nPop=Data_i.pop;    nObW=(nPop-nWh)/nWh;     VarMin=Data_i.lb;      VarMax=Data_i.ub;  
CostFunction=Data_i.fobj;     nOb=nPop-nWh;   nVar=Data_i.dim;   MaxDecades=Data_i.maxIter;
global ProblemSettings;
ProblemSettings.CostFunction=CostFunction;
ProblemSettings.nVar=nVar;
ProblemSettings.VarMin=VarMin;
ProblemSettings.VarMax=VarMax;

global TFWOSettings;
TFWOSettings.nPop=nPop;
TFWOSettings.nWh=nWh;
TFWOSettings.nOb=nOb;
TFWOSettings.nObW=nObW;
TFWOSettings.MaxDecades=MaxDecades;

Whirlpool=Initialize(Data_i,Data_o);
BestSol.Position=[];
BestSol.Cost=[];
IterCurve=zeros(TFWOSettings.MaxDecades,1);
MeanCost=zeros(TFWOSettings.MaxDecades,1);

for Decade=1:MaxDecades
    Whirlpool=Effectsofwhirlpools(Whirlpool, Decade);%Pseudocodes From 1 To 5    
    Whirlpool=Pseudocode6(Whirlpool, Decade); %Pseudocode 6th   
    WhirlpoolCost=[Whirlpool.Cost];
    [BestWhirlpoolCost BestWhirlpoolIndex]=min(WhirlpoolCost);
    BestWhirlpool=Whirlpool(BestWhirlpoolIndex);
    BestSol.Position=BestWhirlpool.Position;
    BestSol.Cost=BestWhirlpool.Cost;
    S(Decade,:)=BestSol.Position;
    IterCurve(Decade)=BestWhirlpoolCost;
    MeanCost(Decade)=mean(WhirlpoolCost);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>BestWhirlpoolCost
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestWhirlpoolCost;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=BestWhirlpoolIndex;
    end    
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Decade;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve';
end

function Whirlpool=Initialize(Data_i,Data_o)
    global ProblemSettings;
    global TFWOSettings;
    
    CostFunction=ProblemSettings.CostFunction;
    nVar=ProblemSettings.nVar;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    
    nPop=TFWOSettings.nPop;
    nWh=TFWOSettings.nWh;
    nOb=TFWOSettings.nOb;
    nObW=TFWOSettings.nObW;


    EmptyObject.Position=[];
    EmptyObject.Cost=[];
    EmptyObject.delta=[];
    Objects=repmat(EmptyObject,nPop,1);
    for k=1:nPop
        Objects(k).Position=Data_i.X(k,:);
        Objects(k).Cost=Data_i.F_value(k);
        Objects(k).delta=0;
    end
    [SortedCosts CostsSortOrder]=sort([Objects.Cost]);
    Objects=Objects(CostsSortOrder);
    
    EmptyWhirlpool.Position=[];
    EmptyWhirlpool.Cost=[];
    EmptyWhirlpool.TotalCost=[];
    EmptyWhirlpool.nObW=[];
    EmptyWhirlpool.delta=[];
    EmptyWhirlpool.Objects=[];
    Whirlpool=repmat(EmptyWhirlpool,nWh,1);
    for i=1:nWh
        Whirlpool(i).Position=Objects(i).Position;
        Whirlpool(i).Cost=Objects(i).Cost;
        Whirlpool(i).delta=Objects(i).delta;
    end
    
    Objects=Objects(nWh+1:end);
    if isempty(Objects)
        return;
    end
    
    WhirlpoolCosts=[Whirlpool.Cost];

    Objects=Objects(randperm(nOb));
    
    for i=1:nWh
         Whirlpool(i).nObW=nObW;
         Whirlpool(i).Objects=Objects(1:nObW);
         Objects=Objects(nObW+1:end);
    end
end

function Whirlpool=Effectsofwhirlpools(Whirlpool, Decade)
    global ProblemSettings;
    global TFWOSettings;
    CostFunction=ProblemSettings.CostFunction;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    nVar=ProblemSettings.nVar;
    
    for i=1:numel(Whirlpool)
        for j=1:Whirlpool(i).nObW
            if numel(Whirlpool)~=1
                J=[];
                S=[];
                E=[];
                D=[];
                AA = 1:numel(Whirlpool);
                AA(i)=[];
                for t=1:AA
                    J(t)=(abs(Whirlpool(t).Cost)^1)*((abs(sum(Whirlpool(t).Position))-(sum(Whirlpool(i).Objects(j).Position)))^0.5);
                end
                  %%%%%%%%%%%%% min
                S=min(J);
                [ E D]=find(S==J);
                
                d=rand(1, nVar).*(Whirlpool(D(1)).Position-Whirlpool(i).Objects(j).Position);
                
                %%%%%%%%%%%%% max
                S2=max(J);
                [ E2 D2]=find(S2==J);
                d2=rand(1, nVar).*(Whirlpool(D2(1)).Position-Whirlpool(i).Objects(j).Position);  
            end
            
            if numel(Whirlpool)==1
                d=rand(1, nVar).*(Whirlpool(i).Position-Whirlpool(i).Objects(j).Position);
                d2=0;
                D(1)=i; 
            end
            Whirlpool(i).Objects(j).delta=Whirlpool(i).Objects(j).delta+ (rand)*rand*pi;
            eee= Whirlpool(i).Objects(j).delta;
            fr0=(cos(eee));
            fr10=(-sin(eee));
            
            x=((fr0.*(d))+(fr10.*(d2)))*(1+abs(fr0*fr10*1));
            RR=(Whirlpool(i).Position-x);
            RR =min(max(RR,VarMin),VarMax);
            Cost=CostFunction(RR ) ;
            if  Cost<= Whirlpool(i).Objects(j).Cost
                Whirlpool(i).Objects(j).Cost=Cost;
                Whirlpool(i).Objects(j).Position=RR;
            end
            %%%%%%%%%%Pseudo-code 3:
            FE_i=(abs(cos(Whirlpool(i).Objects(j).delta)^2*sin(Whirlpool(i).Objects(j).delta)^2))^2;
            
            %              Q=Q^(2);
            if rand<(FE_i)
                k=randi([1 nVar]);
                Whirlpool(i).Objects(j).Position(k)=unifrnd(VarMin(k),VarMax(k));
                Whirlpool(i).Objects(j).Cost=CostFunction(Whirlpool(i).Objects(j).Position);
            end
        end
    end
    %%%%%%%%%% Pseudo-code 4:
    J2=[];
    for t=1:numel(Whirlpool)
        J2(t)=(Whirlpool(t).Cost);
    end
    S2=min(J2);
    [E2 D2]=find(S2==J2);
    d2=Whirlpool(D2(1)).Position;
    for i=1:numel(Whirlpool)
        J=[];
        E=[];
        D=[];
        for t=1:numel(Whirlpool)
            J(t)=Whirlpool(t).Cost*(abs((sum(Whirlpool(t).Position))-(sum(Whirlpool(i).Position))));
            if t==i
                J(t)=inf;
            end
        end
        S=min(J);
        [E D]=find(S==J);
        %%%%%%%%%%
        Whirlpool(i).delta=Whirlpool(i).delta+ (rand)*rand*pi;
        d=Whirlpool(D(1)).Position-Whirlpool(i).Position;
        fr=abs(cos(Whirlpool(i).delta)+sin(Whirlpool(i).delta));
        x= fr*rand(1, nVar).*(d);
        Whirlpool1(i).Position=Whirlpool(D(1)).Position-x;   
        Whirlpool1(i).Position=min(max(Whirlpool1(i).Position,VarMin),VarMax);
        Whirlpool1(i).Cost=CostFunction(Whirlpool1(i).Position);
        %%%%%%Pseudo-code 5:%%selection Whirlpool
        if Whirlpool1(i).Cost<=Whirlpool(i).Cost
            Whirlpool(i).Position= Whirlpool1(i).Position;
            Whirlpool(i).Cost= Whirlpool1(i).Cost; 
        end
    end
    if S2<Whirlpool(D2(1)).Cost
        Whirlpool(i).Position=  d2;
        Whirlpool(i).Cost= S2; 
    end
end

function Whirlpool=Pseudocode6(Whirlpool, Decade)
    global ProblemSettings;
    global TFWOSettings;
   
    CostFunction=ProblemSettings.CostFunction;
    nVar=ProblemSettings.nVar;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    
    nPop=TFWOSettings.nPop;
    nWh=TFWOSettings.nWh;
    nOb=TFWOSettings.nOb;
    for i=1:numel(Whirlpool)
        cc=[Whirlpool(i).Objects.Cost];
        [min_cc min_cc_index]=min(cc);
        if min_cc<=Whirlpool(i).Cost
            BestObject=Whirlpool(i).Objects(min_cc_index);
            Whirlpool(i).Objects(min_cc_index).Position=Whirlpool(i).Position;
            Whirlpool(i).Objects(min_cc_index).Cost=Whirlpool(i).Cost;
            Whirlpool(i).Position=BestObject.Position;
            Whirlpool(i).Cost=BestObject.Cost; 
        end
    end
end