function Data_o = MVO(Data_i,Data_o)
rand_num=[];                        
ti=clock;                              
LightRays = Data_i.X;               
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
Best_universe_Inflation_rate=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
Best_universe=Data_i.X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
Universes=Data_i.X;

WEP_Max=1;
WEP_Min=0.2;

Time=1;
Max_time=Data_i.maxIter;
while Time<Max_time+1
    %Eq. (3.3) in the paper
    WEP=WEP_Min+Time*((WEP_Max-WEP_Min)/Max_time);
    % Eq. (3.4) in the paper
    TDR=1-((Time)^(1/6)/(Max_time)^(1/6));
    Inflation_rates=zeros(1,size(Universes,1));

    for i=1:size(Universes,1)
        %Boundary checking (to bring back the universes inside search
        % space if they go beyoud the boundaries
        Flag4ub=Universes(i,:)>Data_i.ub;
        Flag4lb=Universes(i,:)<Data_i.lb;
        Universes(i,:)=(Universes(i,:).*(~(Flag4ub+Flag4lb)))+Data_i.ub.*Flag4ub+Data_i.lb.*Flag4lb;
        Inflation_rates(1,i)=Data_i.fobj(Universes(i,:));
        %Elitism
    end

    [sorted_Inflation_rates,sorted_indexes]=sort(Inflation_rates);
    for newindex=1:Data_i.pop
        Sorted_universes(newindex,:)=Universes(sorted_indexes(newindex),:);
    end
    % Eq. (3.1) 
    normalized_sorted_Inflation_rates=normr(sorted_Inflation_rates);
    Universes(1,:)= Sorted_universes(1,:);
    
    for i=2:size(Universes,1)%Starting from 2 since the firt one is the elite
        Back_hole_index=i;
        for j=1:size(Universes,2)
            r1=rand();
            if r1<normalized_sorted_Inflation_rates(i)
                White_hole_index=RouletteWheelSelection(-sorted_Inflation_rates);% for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates
                if White_hole_index==-1
                    White_hole_index=1;
                end
                %Eq. (3.1) in the paper
                Universes(Back_hole_index,j)=Sorted_universes(White_hole_index,j);
            end
            
            if (size(Data_i.lb,2)==1)
                %Eq. (3.2) in the paper if the boundaries are all the same
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((Data_i.ub-Data_i.lb)*rand+Data_i.lb);
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((Data_i.ub-Data_i.lb)*rand+Data_i.lb);
                    end
                end
            end
            
            if (size(Data_i.lb,2)~=1)
                %Eq. (3.2) in the paper if the upper and lower bounds are
                %different for each variables
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((Data_i.ub(j)-Data_i.lb(j))*rand+Data_i.lb(j));
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((Data_i.ub(j)-Data_i.lb(j))*rand+Data_i.lb(j));
                    end
                end
            end
        end
    end
    
    for i=1:Data_i.pop
        if Inflation_rates(1,i)<Best_universe_Inflation_rate
            Best_universe_Inflation_rate=Inflation_rates(1,i);
            Best_universe=Universes(i,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_universe_Inflation_rate;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end    

    IterCurve(Time)=Best_universe_Inflation_rate;
    Time=Time+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Time;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function choice = RouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand() * accumulation(end);
  chosen_index = -1;
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end