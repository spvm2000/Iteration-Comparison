function Data_o = MSA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
G = Data_i.maxIter;        
SearchAgents_no = Data_i.pop;         
d = Data_i.dim;               
Moth_pos = Data_i.X;                   
Moth_fitness = Data_i.F_value;     
ub=Data_i.ub;
lb=Data_i.lb;
Nc=6;% Number of Pathfinders: 4 <=  Nc  <= 20% of SearchAgents_no

%Initialize the fitnesses and positions of moths
Moth_fitness=Moth_fitness';
[Moth_fitness ,location]=sort(Moth_fitness);
Moth_pos=Moth_pos(location,:);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Moth_fitness(1);
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Moth_pos(1);

Convergence_curve=inf*zeros(1,Data_i.maxIter);
g=1;% Loop counter

% Main loop
while g<(G+1) 
pMoth_pos=Moth_pos;
pMoth_fitness=Moth_fitness; 
%____________________ 1. Reconnaissance phase : Pathfinder moths ____________________________________________________    

Lights=Moth_pos(1:Nc,:);   Light_fitness=Moth_fitness(1:Nc); % best moths considered to be pathfinders

[Lights,Light_fitness] = Pathfinders(Lights,Light_fitness,Nc,Data_i.fobj,ub,lb); % improve pathfinders

Moth_fitness(1:Nc)=Light_fitness;      Moth_pos(1:Nc,:)=Lights; % insertion in the swarm

% Sharing of luminescence intensities using Eqs.(28) and (29). 
for i=1:Nc
    if Light_fitness(i)>=0
        Light_fitness(i)=1/(1+Light_fitness(i));
    else
        Light_fitness(i)=1+abs(Light_fitness(i));
    end
end

for i=1:Nc
    probability(i)=Light_fitness(i)/sum(Light_fitness);
end %#ok<*AGROW>

[~,R] = histc(rand(1,SearchAgents_no),cumsum([0;probability(:)]));a1 = 1:Nc;new_Light = Lights(a1(R),:);

%__________________2. Transverse orientation phase : prospector moths_________________________________    
% inspired from Moth-flame Optimization (MFO) [14]: with each variable as an integrated unit

Predictor_no=round((SearchAgents_no-Nc)*(1-g/(G))); %No. of prospectors using Eq.(30).
a=-1-g/G;

for i=Nc+1:Predictor_no+Nc %#ok<*NBRAK>
       t=(a-1)*rand()+1; spiral=exp(t).*cos(t.*2*pi);    
       D_to_Light=abs(new_Light(i,:)-Moth_pos(i,:));
       Moth_pos(i,:)=D_to_Light.*spiral+new_Light(i,:);% updating prospectors using Eq.(31)          
end

% Return back prospectors that go beyond the boundaries
Moth_pos(Nc+1:Predictor_no+Nc,:) =Bound_Checking(Moth_pos(Nc+1:Predictor_no+Nc,:),ub,lb);

%Fitness evaluation for prospectors
for i2 = Nc+1:Predictor_no+Nc
 Moth_fitness(i2,:) = Data_i.fobj(Moth_pos(i2,:));
end

[~, location] = min(Moth_fitness);    gx = Moth_pos(location,:);

%_______________________ 3. Celestial navigation phase: onlooker moths____________________
% 

% 3.1 Gaussian walks
 x3=round((SearchAgents_no-Predictor_no-Nc)*1/2);j2=SearchAgents_no;

 for j1=1:x3%#ok<*ALIGN> 
    j2=j1+Predictor_no+Nc;
    % updating prospectors using Eqs.(32-33)
    Moth_pos(j2,:)=normrnd(gx(1,:),(log(g)/g)*(abs((Moth_pos(j2,:)-gx(1,:)))), [1 size(Moth_pos,2)])+(randn*gx(1,:)-randn*Moth_pos(j2,:));%gaussian(large step away)
 
    out=Moth_pos(j2,:)<lb | Moth_pos(j2,:) > ub;
    Moth_pos(j2,out)=normrnd(pMoth_pos(j2,out),(abs((rand*pMoth_pos(j2,out)-rand*gx(1,out)))), [1 size(Moth_pos(j2,out),2)]);%(step around);    
 end
   
% 3.2 Associative learning mechanism with immediate memory (ALIM)
x3=SearchAgents_no-j2; % no. of onlookers moves with ALIM

for j1=1:x3%#ok<*ALIGN> 
    j3=j1+j2;
    % updating prospectors using Eq.(34)
    Moth_pos(j3,:)=Moth_pos(j3,:)+0.001*unifrnd(lb-Moth_pos(j3,:),ub-Moth_pos(j3,:))+(2*g/G)*rand(size(gx)).*(gx(1,:)-Moth_pos(j3,:))+(1-g/G)*rand(size(gx)).*(new_Light(j3,:)-Moth_pos(j3,:));%step away like PSO
end

% Apply Position Limits for onlookers
Moth_pos(Predictor_no+Nc+1:SearchAgents_no,:) = Bound_Checking(Moth_pos(Predictor_no+Nc+1:SearchAgents_no,:),ub,lb);

%Fitness evaluation for onlookers 
for i2 = Predictor_no+Nc+1:SearchAgents_no
Moth_fitness(i2,:) = Data_i.fobj(Moth_pos(i2,:));
end
%===================================================================== 
% select the best moths
[Moth_fitness,I]=unique([Moth_fitness, pMoth_fitness],'first');  Moth_fitness=Moth_fitness(1:SearchAgents_no);
dMoth_pos=[Moth_pos; pMoth_pos];Moth_pos=dMoth_pos(I(1:SearchAgents_no),:);
% Find the global best solution
[Best_score, location] = min(Moth_fitness);    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=location;
Convergence_curve(g)=Best_score;
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Best_score;
g=g+1; 
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=g;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end


%____________________ 1. Reconnaissance phase : Pathfinder moths ____________________________________________________    
function [Lights,Light_fitness] = Pathfinders(Lights,Light_fitness,Nc,fobj,ub,lb)

    Lights1=Lights;%trail vector
    
    %Proposed adaptive crossover operation based on population diversity
    
    C_V=std(Lights,0,1)./abs(mean(Lights));% sensing_distance using Eqs.(18) and (19)
    nmu=mean(C_V);mcol=find(C_V <= nmu);% Select crossover points using condition of Eq.(20)
    
    %levy mutations for crossover using Eq.(24)
    beta=3/2;sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(Nc,length(mcol),2)*sigma;v=randn(Nc,length(mcol),2);step=u./abs(v).^(1/beta);L=abs(0.01*step);% Draw Levy flight sample
    
    % Proposed difference vectors L関y-mutation 
    for i=1:Nc
        if Nc<6; I=randperm(Nc);  II=find(I ~= i);I=I(II(1:3));
        Lights1(i,mcol)=Lights(I(1),mcol)+L(i,:,1).*(Lights(I(2),mcol)-Lights(I(3),mcol));%揇E/rand/1�
        else
        I=randperm(Nc);  II=find(I ~= i);I=I(II(1:5));%mating vectors (donners)
        Lights1(i,mcol)=Lights(I(1),mcol)+(L(i,:,1).*(Lights(I(2),mcol)-Lights(I(3),mcol))+L(i,:,2).*(Lights(I(4),mcol)-Lights(I(5),mcol)));%揇E/rand/2�:using Eq.(25) 
        end
    end
    
    % Return back the Pathfinders that go beyond the boundaries
    Lights1 =Bound_Checking(Lights1,ub,lb);

%selection for updating Pathfinder moths using Eq.(27)
    for i=1:size(Lights1,1)
    Light_fitness1(i,:) = fobj(Lights1(i,:)); 
        if Light_fitness1(i,1)< Light_fitness(i,1)
            Light_fitness(i,:)=Light_fitness1(i,:);
            Lights(i,:)=Lights1(i,:);
        end
    end
end

% return back populations
function Moth_pos = Bound_Checking(Moth_pos,ub,lb)
    for i=1:size(Moth_pos,1)
        Flag4ub=Moth_pos(i,:)>ub; 
        Flag4lb=Moth_pos(i,:)<lb;
        Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end 
end