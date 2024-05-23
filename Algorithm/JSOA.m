function Data_o = JSOA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
parameters.maxIteration = Data_i.maxIter;      
parameters.SearchAgents = Data_i.pop;            
parameters.dim=Data_i.dim;
Positions = Data_i.X;               
Fitness = Data_i.F_value;         
parameters.lb = Data_i.lb;
parameters.ub = Data_i.ub;

[vMin ,minIdx]= min(Fitness);           
[vMax ,maxIdx]= max(Fitness); 
theBestVct=  Positions(minIdx,:);       % the best vector 
theWorstVct= Positions(maxIdx,:);       % the worst vector
Convergence_curve=zeros(1,parameters.maxIteration);
Convergence_curve(1)= vMin;

Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=vMin;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=minIdx;

pheromone= getPheromone( Fitness, vMin, vMax);  % Used in Jumping spider pheromone rates Eq.9
gravity= 9.80665; %m/seg^2  
vo= 100; % 100 mm/seg
for t=1:parameters.maxIteration       
   for r=1:parameters.SearchAgents
       if rand() < 0.5  % Attack, persecution and jumping on the prey 
               if rand()  < 0.5
                     %****************************************************************
                     % Jumping on the prey represented by equation of
                     % projectile motion. Eq. 6 on paper 
                     radians=  90* rand()*pi/180; 
                     v(r,:)=  (Positions(r,:) .* tan(radians))-((gravity .*(Positions(r,:).^2)) ./ (2*(vo^2)*(cos(radians)^2)));
                     %***************************************
               else
                    % Persecution represented by the uniformly accelerated
                    % rectilinear motion. Eq.2 on paper
                    ban=1;
                     while(ban)
                            r1= round(1+ (parameters.SearchAgents-1)* rand());
                            if (r ~= r1)
                                 ban=0;
                            end 
                     end
                     v(r,:)= 0.5 * (Positions(r,:) -Positions(r1,:)); 
                end
       else   % Searching for prey 
             if rand < 0.5    % Global search. Eq. 8 on paper
                  e1=  CauchyRand(0,1);
                  v(r,:)=  theBestVct +( theBestVct-theWorstVct)  * e1;                
             else             % Local search. Eq. 7 on paper 
                 walk= -2 + 4 * rand(); % -2 < d < 2   Uniformly distributed pseudorandom numbers
                 e=  randn(); % Normally distributed pseudorandom numbers.
                  v(r,:)= theBestVct + walk*(0.5-e);
             end
       end       
       if pheromone(r) <= 0.3  % Jumping spider pheromone rates. Eq.10, Algorithm 1 on paper.
             band=1;
             while band
               r1= round(1+ (parameters.SearchAgents-1)* rand());
               r2= round(1+ (parameters.SearchAgents-1)* rand());
               if r1 ~= r2
                   band=0;
               end
             end  % Eq.10 on paper
                  v(r,:)=   theBestVct + (Positions(r1,:)-((-1)^getBinary)*Positions(r2,:))/2;
       end    
     %****************************************************************
        Flag4ub=v(r,:)>parameters.ub;
        Flag4lb=v(r,:)<parameters.lb;
        v(r,:)=(v(r,:).*(~(Flag4ub+Flag4lb)))+parameters.ub.*Flag4ub+parameters.lb.*Flag4lb;
    % Evaluate new solutions
     Fnew= Data_i.fobj(v(r,:));
     % Update if the solution improves
     if Fnew <= Fitness(r)
        Positions(r,:)= v(r,:);
        Fitness(r)= Fnew;
     end
     if Fnew <= vMin
         theBestVct= v(r,:);   
         vMin= Fnew;
         if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>vMin
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=vMin;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=r;
         end    
     end 
   end
   [vMax ,maxIdx]= max(Fitness); 
   theWorstVct= Positions(maxIdx,:);
   % update max and pheromons
   pheromone= getPheromone( Fitness, vMin, vMax);
   Convergence_curve(t)= vMin; 
end

Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end
%***********************************[End JSOA Algorithm]
%*******************************************************
function [ o ] =  getPheromone(  fit, min, max )
    for i=1:size(fit,2)
         o(i)= (max-fit(i))/(max-min);
    end
end

function [ val] = getBinary( )
    if rand() < 0.5
         val= 0;
    else
         val=1;
    end
end

function [cauchy] = CauchyRand(m,c)
    cauchy = c*tan(pi*(rand()-0.5)) + m;
end