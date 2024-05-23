function Data_o = BWOA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
Positions = Data_i.X;                       
Fitness=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
SearchAgents_no=Data_i.pop;     Max_iter=Data_i.maxIter;   
lb=Data_i.lb;       ub=Data_i.ub; dim=Data_i.dim;  fobj=Data_i.fobj;
[vMin minIdx]= min(Fitness);  % the min fitness value vMin and the position minIdx
theBestVct= Positions(minIdx,:);  % the best vector 
[vMax maxIdx]= max(Fitness);  % the min fitness value vMin and the position minIdx
Convergence_curve=zeros(1,Max_iter);
Convergence_curve(1)= vMin;
pheromone= getPheromone( Fitness, vMin, vMax);
IterCurve=zeros(1,Data_i.maxIter);
for t=1:Max_iter
    beta= -1 + 2* rand();  % -1 < beta2 < 1     section 3.2.1 in the paper, page 4
    m= 0.4 + 0.5 *rand();
    for r=1:SearchAgents_no
        P= rand();
        r1= round(1+ (SearchAgents_no-1)* rand());
        if P >= 0.3  % spiral search   Eq. 11 in the paper, page 4
             v(r,:)= theBestVct -cos(2*pi*beta)*Positions(r,:);
        else         % direct search Eq. 11
             v(r,:)=theBestVct -m *Positions(r1,:);
        end
        if pheromone(r) <= 0.3
            band=1;
            while band 
                r1= round(1+ (SearchAgents_no-1)* rand());
                r2= round(1+ (SearchAgents_no-1)* rand());
                if r1 ~= r2 
                   band=0;
               end
            end
            v(r,:)=theBestVct + (Positions(r1,:)-((-1)^getBinary)*Positions(r2,:))/2;
        end
        Flag4ub=v(r,:)>ub;
        Flag4lb=v(r,:)<lb;
        v(r,:)=(v(r,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        Fnew= fobj(v(r,:));
        if Fnew <= Fitness(r);
            Positions(r,:)= v(r,:);
            Fitness(r)= Fnew;
        end
        if Fnew <= vMin
            theBestVct=v(r,:);
            vMin= Fnew;
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>vMin
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=vMin;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=r;
            end    
        end
    end
    [vMax maxIdx]= max(Fitness);
    pheromone= getPheromone( Fitness, vMin, vMax);
    IterCurve(t)= vMin;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [ val] = getBinary( )  
    if rand() < 0.5
         val= 0;
    else
         val=1;
    end
end

function [ o ] = getPheromone(fit, min, max ) % Eq.12 in the paper
    for i=1:size(fit,2)
         o(i)= (max-fit(i))/(max-min);
    end
end