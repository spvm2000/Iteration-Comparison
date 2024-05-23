function Data_o = POAA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  

%% Problem Definition
MaxIterations = Data_i.maxIter;     
NumAgents = Data_i.pop;            
Dim = Data_i.dim;              
Pos = Data_i.X;               
Fit = Data_i.F_value;         
LowerBound = Data_i.lb;
UpperBound = Data_i.ub;
    
    ConvergenceCurve=zeros(1,MaxIterations);
        
    NumPeacock=5; % the number of leader (Peacock)
    NumPeahen=round((NumAgents-NumPeacock)*0.3); % the number of peahen
    NumPeacockCub=NumAgents-NumPeacock-NumPeahen; % the number of peacock cub

    SearchRadius0=(UpperBound-LowerBound)*0.2; % initial dance radius of peacock

    % initialization
    empty_peacock.Position=[];
    empty_peacock.Fitness=[];

    PeacockPopulation0=repmat(empty_peacock,[NumAgents,1]);
    Peahen=repmat(empty_peacock,[NumPeahen,1]);
    PeacockCub=repmat(empty_peacock,[NumPeacockCub,1]);

    for k=1:NumAgents
        PeacockPopulation0(k).Position=Pos(k,:);
        PeacockPopulation0(k).Fitness=Fit(k);
    end

    PeacockPopulation=PeacockPopulation0;
    [~,index]=sort([PeacockPopulation.Fitness]);
    PeacockPopulation=PeacockPopulation(index);

    ConvergenceCurve(1)=PeacockPopulation(1).Fitness;
    [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Fit);

    % main loop
    for it=2:MaxIterations

        SearchRadius=SearchRadius0-(SearchRadius0-0)*(it/MaxIterations)^0.01;
        alpha=0.9-(0.9-0.4)*(it/MaxIterations)^2;
        delta=0.1+(1-0.1)*(it/MaxIterations)^0.5;
        step=0.1+(1-0.1)*(it/MaxIterations);

        Peacock=PeacockPopulation(1:NumPeacock);
        if rand<1
            X_random=2*rand(1,Dim)-1;
            Peacock(1).Position=Peacock(1).Position+1*SearchRadius.*X_random/(eps+norm(X_random));
        end
        if rand<0.9
            X_random=2*rand(1,Dim)-1;
            Peacock(2).Position=Peacock(2).Position+1.5*SearchRadius.*X_random/(eps+norm(X_random));
        end
        if rand<0.8
            X_random=2*rand(1,Dim)-1;
            Peacock(3).Position=Peacock(3).Position+2*SearchRadius.*X_random/(eps+norm(X_random));
        end
        if rand<0.6
            X_random=2*rand(1,Dim)-1;
            Peacock(4).Position=Peacock(4).Position+3*SearchRadius.*X_random/(eps+norm(X_random));
        end
        if rand<0.3
            X_random=2*rand(1,Dim)-1;
            Peacock(5).Position=Peacock(5).Position+5*SearchRadius.*X_random/(eps+norm(X_random));
        end
        for k=1:NumPeacock
            flag4ub=Peacock(k).Position>UpperBound;
            flag4lb=Peacock(k).Position<LowerBound;
            Peacock(k).Position=~(flag4ub+flag4lb).*Peacock(k).Position+flag4ub.*UpperBound+flag4lb.*LowerBound;
            Peacock(k).Fitness=Data_i.fobj(Peacock(k).Position);
            if Peacock(k).Fitness < PeacockPopulation(k).Fitness
                PeacockPopulation(k)=Peacock(k);
            end
        end

        for k=1:NumPeahen
            r1=rand();
            if r1 <= 1 && r1 >=0.6
                Peahen(k).Position=PeacockPopulation(NumPeacock+k).Position+3*step*(PeacockPopulation(1).Position-PeacockPopulation(NumPeacock+k).Position);
            end
            if r1 < 0.6 && r1 >=0.4
                Peahen(k).Position=PeacockPopulation(NumPeacock+k).Position+3*step*(PeacockPopulation(2).Position-PeacockPopulation(NumPeacock+k).Position);
            end
            if r1 < 0.4 && r1 >=0.2
                Peahen(k).Position=PeacockPopulation(NumPeacock+k).Position+3*step*(PeacockPopulation(3).Position-PeacockPopulation(NumPeacock+k).Position);
            end
            if r1 < 0.2 && r1 >=0.1
                Peahen(k).Position=PeacockPopulation(NumPeacock+k).Position+3*step*(PeacockPopulation(4).Position-PeacockPopulation(NumPeacock+k).Position);
            end
            if r1 < 0.1 && r1 >=0
                Peahen(k).Position=PeacockPopulation(NumPeacock+k).Position+3*step*(PeacockPopulation(5).Position-PeacockPopulation(NumPeacock+k).Position);
            end
            flag4ub=Peahen(k).Position>UpperBound;
            flag4lb=Peahen(k).Position<LowerBound;
            Peahen(k).Position=~(flag4ub+flag4lb).*Peahen(k).Position+flag4ub.*UpperBound+flag4lb.*LowerBound;
            Peahen(k).Fitness=Data_i.fobj(Peahen(k).Position);
            if Peahen(k).Fitness < PeacockPopulation(NumPeacock+k).Fitness
                PeacockPopulation(NumPeacock+k)=Peahen(k);
            end
        end

        for k=1:NumPeacockCub
            PeacockCub(k)=PeacockPopulation(NumPeacock+NumPeahen+k);
            
            r2=rand;
            if r2>0.8 && r2<=1
                SelectedPeacock=PeacockPopulation(1);
            elseif r2>0.6 && r2<=0.8
                SelectedPeacock=PeacockPopulation(2);
            elseif r2>0.4 && r2<=0.6
                SelectedPeacock=PeacockPopulation(3);
            elseif r2>0.2 && r2<=0.4
                SelectedPeacock=PeacockPopulation(4);
            else 
                SelectedPeacock=PeacockPopulation(5);
            end
            
            PeacockCub(k).Position=PeacockCub(k).Position+alpha*Levy(Dim).*( PeacockPopulation(1).Position - PeacockCub(k).Position )+delta*( SelectedPeacock.Position - PeacockCub(k).Position );
            
            flag4ub=PeacockCub(k).Position>UpperBound;
            flag4lb=PeacockCub(k).Position<LowerBound;
            PeacockCub(k).Position=~(flag4ub+flag4lb).*PeacockCub(k,:).Position+flag4ub.*UpperBound+flag4lb.*LowerBound;
            PeacockCub(k).Fitness=Data_i.fobj(PeacockCub(k).Position);
            if PeacockCub(k).Fitness < PeacockPopulation(NumPeacock+NumPeahen+k).Fitness
                PeacockPopulation(NumPeacock+NumPeahen+k)=PeacockCub(k,:);
            end
        end

        Peacock=PeacockPopulation(1:NumPeacock);
        
        Xrandom=2*rand(1,Dim)-1;
        Direction1=Peacock(1,:).Position-Peacock(2,:).Position;
        Direction2=Xrandom-(Xrandom*Direction1')/(Direction1*Direction1'+eps)*Direction1;
        Direction2=Direction2/norm(Direction2+eps)*norm(Direction1);
        Peacock(2,:).Position=Peacock(2,:).Position+step*Direction1+rand*Direction2;
        
        Xrandom=2*rand(1,Dim)-1;
        Direction1=Peacock(1,:).Position-Peacock(3,:).Position;
        Direction2=Xrandom-(Xrandom*Direction1')/(Direction1*Direction1'+eps)*Direction1;
        Direction2=Direction2/norm(Direction2+eps)*norm(Direction1);
        Peacock(3,:).Position=Peacock(3,:).Position+step*Direction1+rand*Direction2;
        
        Xrandom=2*rand(1,Dim)-1;
        Direction1=Peacock(1,:).Position-Peacock(4,:).Position;
        Direction2=Xrandom-(Xrandom*Direction1')/(Direction1*Direction1'+eps)*Direction1;
        Direction2=Direction2/norm(Direction2+eps)*norm(Direction1);
        Peacock(4,:).Position=Peacock(4,:).Position+step*Direction1+rand*Direction2;
        
        Xrandom=2*rand(1,Dim)-1;
        Direction1=Peacock(1,:).Position-Peacock(5,:).Position;
        Direction2=Xrandom-(Xrandom*Direction1')/(Direction1*Direction1'+eps)*Direction1;
        Direction2=Direction2/norm(Direction2+eps)*norm(Direction1);
        Peacock(5,:).Position=Peacock(5,:).Position+step*Direction1+rand*Direction2;
        
        for k=1:NumPeacock
            flag4ub=Peacock(k).Position>UpperBound;
            flag4lb=Peacock(k).Position<LowerBound;
            Peacock(k).Position=~(flag4ub+flag4lb).*Peacock(k).Position+flag4ub.*UpperBound+flag4lb.*LowerBound;
            Peacock(k).Fitness=Data_i.fobj(Peacock(k).Position);
            if Peacock(k).Fitness < PeacockPopulation(k).Fitness
                PeacockPopulation(k)=Peacock(k);
            end
        end

        [~,index]=sort([PeacockPopulation.Fitness]);
        PeacockPopulation=PeacockPopulation(index);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
        ConvergenceCurve(1,it)=PeacockPopulation(1).Fitness;
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = PeacockPopulation(1).Fitness;
    end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=ConvergenceCurve;  
end

function o=Levy(d)
beta=3/2;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=rands(1,d)*sigma;
v=rand(1,d);
step=u./abs(v).^(1/beta);
% Eq. (3.9)
o=step;
end

