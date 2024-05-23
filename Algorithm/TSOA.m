function Data_o = TSOA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
X = Data_i.X;                      
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
ns=Data_i.pop;      SN=Data_i.dim;  CostFunction=Data_i.fobj;   nV=Data_i.dim;  ConstraintFunction=Data_i.fobj;

Empty.Location = [];
Empty.Cost = inf;
Galaxy_Center = repmat (Empty, 1, 1);
region = repmat (Empty, ns*SN, 1);
selested_regions = repmat (Empty, ns, 1);
Stars = repmat (Empty, ns, 1);
Stars_sorted = zeros(ns,1);
Ranks = 1:1:ns;
Stars_Ranks = zeros(ns,1);
Luminosity = zeros(ns,1);
Star_RanksNormal = zeros(ns,1);
Distance = zeros(ns,1);
Transit0 = zeros(ns,1);
SN_P = repmat (Empty, SN, 1);
Bests=region;


Vmin=Data_i.lb;
Vmax=Data_i.ub;

        
%% Galaxy Phase

% Initial Location of The Center of the Galaxy
for i=1:Data_i.pop
    Galaxy_Center.Location = X(i,:);
    Galaxy_Center.Cost = Cost(i);
end
% Galactic Habitate Zone of the Galaxy
for l = 1:(ns*SN)
    zone = randi(2);
    if zone ==1
        difference = rand().*(Galaxy_Center.Location)-(unifrnd(Vmin,Vmax,1,nV));
    else
        difference = rand().*(Galaxy_Center.Location)+(unifrnd(Vmin,Vmax,1,nV));
    end
    Noise = ((rand(1,nV)).^3).*(unifrnd(Vmin,Vmax,1,nV));
    region(l).Location = Galaxy_Center.Location + difference - Noise;
    region(l).Location = max(region(l).Location, Vmin);
    region(l).Location = min(region(l).Location, Vmax);
    region(l).Cost = Penalty(region(l),ConstraintFunction,CostFunction);
end

% Selection of Stars from the Galactic Habitate Zone of the Galaxy
[Sort,index]=sort([region.Cost]);
for i = 1:ns
    selested_regions(i) = region(index(1,i));
    for k = 1:SN
        zone = randi(2);
        if zone ==1
            difference = rand().*(selested_regions(i).Location)-rand().*(unifrnd(Vmin,Vmax,1,nV));
        else
            difference = rand().*(selested_regions(i).Location)+rand().*(unifrnd(Vmin,Vmax,1,nV));
        end
        Noise = ((rand(1,nV)).^3).*(unifrnd(Vmin,Vmax,1,nV));
        new.Location = selested_regions(i).Location + difference - Noise;
        new.Location = max(new.Location, Vmin);
        new.Location = min(new.Location, Vmax);
        new.Cost = CostFunction(new.Location);
        new.Cost = Penalty(new,ConstraintFunction,CostFunction);
        if new.Cost < Stars(i).Cost
            Stars(i) = new;
        end
    end
end

% Initial Location of the Best Planets (Start Point: Its Star)
Best_Planets = Stars;

% Specification of the Best Planet
[Sort,index]=sort([Best_Planets.Cost]);
Best_Planet = Best_Planets(index(1,1));

% Telescope Location
Telescope.Location = unifrnd(Vmin,Vmax,1,nV);

% Determination of the Luminosity of the Stars
for i = 1:ns
    Stars_sorted(i,1) = Stars(i).Cost;
end
Stars_sorted = sort (Stars_sorted);
for i = 1:ns
    for ii = 1:ns
        if Stars(i).Cost == Stars_sorted(ii,1)
            Stars_Ranks(i,1) = Ranks(1,ii);
            Star_RanksNormal(i,1) = (Stars_Ranks(i,1))./ns;
        end
    end
    Distance(i,1) = sum((Stars(i).Location-Telescope.Location).^2).^0.5;
    Luminosity(i,1) = Star_RanksNormal(i,1)/((Distance(i,1))^2);
end
Luminosity_new = Luminosity;
Stars2 = Stars;

% Loops of the TS Algorithm
for it = 1:Data_i.maxIter
    %% Transit Phase
    Transit = Transit0;
    Luminosity = Luminosity_new;
    
    for i = 1:ns
        difference = (2*rand()-1).*(Stars(i).Location);
        Noise = ((rand(1,nV)).^3).*(unifrnd(Vmin,Vmax,1,nV));
        Stars2(i).Location = Stars(i).Location + difference - Noise;
        Stars2(i).Location = max(Stars2(i).Location, Vmin);
        Stars2(i).Location = min(Stars2(i).Location, Vmax);
        Stars2(i).Cost = Penalty(Stars2(i),ConstraintFunction,CostFunction);
    end
    
    for i = 1:ns
        Stars_sorted(i,1) = Stars2(i).Cost;
    end
    Stars_sorted = sort (Stars_sorted);
    for i = 1:ns
        for ii = 1:ns
            if Stars2(i).Cost == Stars_sorted(ii,1)
                Stars_Ranks(i,1) = Ranks(1,ii);
                Star_RanksNormal(i,1) = (Stars_Ranks(i,1))./ns;
            end
        end
        Distance(i,1) = sum((Stars2(i).Location-Telescope.Location).^2).^0.5;
        Luminosity_new(i,1) = Star_RanksNormal(i,1)/((Distance(i,1))^2);
        if Luminosity_new(i,1) < Luminosity(i,1)
            Transit (i,1) = 1;      % Has transit been observed?  0 = No; 1 = Yes
        end
    end
    Stars = Stars2;
    
    %% Location Phase (Exploration)
    for i = 1:ns
        if Transit (i,1) == 1
            
            % Determination of the Location of the Planet
            Luminosity_Ratio = Luminosity_new(i,1)/Luminosity(i,1);
            Planet.Location = (rand().*Telescope.Location + Luminosity_Ratio.*Stars(i).Location)./2;
            
            for k = 1:SN
                zone = randi(3);
                if zone ==1
                    new.Location = Planet.Location - (2*rand()-1).*(unifrnd(Vmin,Vmax,1,nV));
                elseif zone ==2
                    new.Location = Planet.Location + (2*rand()-1).*(unifrnd(Vmin,Vmax,1,nV));
                else
                    new.Location = Planet.Location + (2.*rand(1,nV)-1).*(unifrnd(Vmin,Vmax,1,nV));
                end
                new.Location = max(new.Location, Vmin);
                new.Location = min(new.Location, Vmax);
                new.Cost = CostFunction(new.Location);
                SN_P(k) = new;
            end
            SUM = 0;
            for k = 1:SN
                SUM = SUM+SN_P(k).Location;
            end
            new.Location = SUM./SN;
            new.Cost = Penalty(new,ConstraintFunction,CostFunction);
            
            if new.Cost < Best_Planets(i).Cost
                Best_Planets(i) = new;
            end
            
        else  % No Transit observed: Neighbouring planets
            
            Neighbor.Location = (rand().*Stars(i).Location + rand().*(unifrnd(Vmin,Vmax,1,nV)))./2;
            
            for k = 1:SN
                zone = randi(3);
                if zone ==1
                    Neighbor.Location = Neighbor.Location - (2*rand()-1).*(unifrnd(Vmin,Vmax,1,nV));
                elseif zone ==2
                    Neighbor.Location = Neighbor.Location + (2*rand()-1).*(unifrnd(Vmin,Vmax,1,nV));
                else
                    Neighbor.Location = Neighbor.Location + (2.*rand(1,nV)-1).*(unifrnd(Vmin,Vmax,1,nV));
                end
                Neighbor.Location = max(Neighbor.Location, Vmin);
                Neighbor.Location = min(Neighbor.Location, Vmax);
                Neighbor.Cost = CostFunction (Neighbor.Location);
                SN_P(k) = Neighbor;
            end
            SUM = 0;
            for k = 1:SN
                SUM = SUM+SN_P(k).Location;
            end
            Neighbor.Location = SUM./SN;
            Neighbor.Cost = Penalty(Neighbor,ConstraintFunction,CostFunction);
            
            if Neighbor.Cost < Best_Planets(i).Cost
                Best_Planets(i) = Neighbor;
            end
        end
        
        if Best_Planets(i).Cost < Best_Planet.Cost
            Best_Planet = Best_Planets(i);
        end
    end
    
    %% Signal Amplification of the Best Planets (Exploitation)
    for i = 1:ns
        for k = 1:SN
            RAND = randi(2);
            if RAND ==1
                Power = randi(SN*ns);
                Coefficient = 2*rand();
                Noise = ((rand(1,nV)).^Power).*(unifrnd(Vmin,Vmax,1,nV));
            else
                Power = randi(SN*ns);
                Coefficient = 2*rand();
                Noise = -((rand(1,nV)).^Power).*(unifrnd(Vmin,Vmax,1,nV));
            end
            chance = randi(2);
            if chance ==1
                new.Location = Best_Planets(i).Location - Coefficient.*Noise;
            else
                new.Location = (rand().*Best_Planets(i).Location) - Coefficient.*Noise;
            end
            new.Location = max(new.Location, Vmin);
            new.Location = min(new.Location, Vmax);
            new.Cost = CostFunction(new.Location);
            new.Cost = Penalty(new,ConstraintFunction,CostFunction);
            
            if new.Cost < Best_Planets(i).Cost
                Best_Planets(i) = new;
            end
        end
        if Best_Planets(i).Cost < Best_Planet.Cost
            Best_Planet = Best_Planets(i);
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Best_Planet.Cost
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Planet.Cost;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end    
        end
    end
    
    %% Results
    Best_Planet.Cost = CostFunction(Best_Planet.Location);
    IterCurve(it) = Best_Planet.Cost;  
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function[Fitness]=Penalty(P0,ConstraintFunction,CostFunction)

%% Step 1: Evaluate solutions infeasibility
fname = fieldnames(P0);
P0.(fname{2}) = CostFunction(P0.(fname{1}));  % Determination of the objective value for P0
Const_P = (ConstraintFunction(P0.(fname{1})))';
n_const = length(Const_P);    % Number of constraints
Constraints_Values = zeros(1,n_const);
Constraints_Values(1,:) = (ConstraintFunction(P0.(fname{1})))';
for j = 1:n_const
    if Constraints_Values(1,j) < 0
        Constraints_Values(1,j) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%% PENALTY %%%%%%%%%%%%%%%%%
SUM = sum(Constraints_Values);
%% Step 2: Calculate the Fitness
S=0;
for j=1:n_const
    S = S+ (Constraints_Values(1,j))^2;
end
if SUM ~=0
    Fitness = P0.(fname{2}) + (10^20)*S;
else
    Fitness = P0.(fname{2});
end
end