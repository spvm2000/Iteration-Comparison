function Data_o = AJSO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
nVar=Data_i.dim;                           % Number of Decision Variables
VarSize=[1 nVar];                  % Size of Decision Variables Matrix
IterCurve=zeros(1,Data_i.maxIter);
%% AJS Parameters
MaxIt = Data_i.maxIter;     
nPop = Data_i.pop;           
popi = Data_i.X;             
popCost = Data_i.F_value;    

%% AJS Main Loop
for it=1:MaxIt
    Meanvl=mean(popi,1);
    [value,index]=sort(popCost);
    BestSol=popi(index(1),:);
    BestCost=popCost(index(1));
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
    for i=1:nPop
        % Calculate time control c(t) using Eq. (17);
        Ar=(1-it*((1)/MaxIt))*(2*rand-1);
        rand_num(it,i)=Ar;
        if abs(Ar)>=0.5
            %% Folowing to ocean current using Eq. (11)
            newsol = popi(i,:)+ rand(VarSize).*(BestSol - 3*rand*Meanvl);
            % Check the boundary using Eq. (19)
            newsol = BoundaryCheck(newsol,Data_i.ub,Data_i.lb,Data_i.dim);
            % Evaluation
            newsolCost = Data_i.fobj(newsol);
            % Comparison
            if newsolCost<popCost(i)
                popi(i,:) = newsol;
                popCost(i)=newsolCost;
                if popCost(i) < BestCost
                    BestCost=popCost(i);
                    BestSol = popi(i,:);
                    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                end
            end
        else
            %% Moving inside swarm
            if rand<=(1-Ar)
                % Determine direction of jellyfish by Eq. (15)
                j=i;
                while j==i
                    j=randperm(nPop,1);
                end
                Step = popi(i,:) - popi(j,:);
                if popCost(j) < popCost(i)
                    Step = -Step;
                end
                % Active motions (Type B) using Eq. (16)
                newsol = popi(i,:) + rand(VarSize).*Step;
            else
                % Passive motions (Type A) using Eq. (12)
                newsol = popi(i,:) + 0.1*(Data_i.ub-Data_i.lb)*rand;
            end
            % Check the boundary using Eq. (19)
            newsol = BoundaryCheck(newsol,Data_i.ub,Data_i.lb,Data_i.dim);
            % Evaluation
            newsolCost = Data_i.fobj(newsol);
            % Comparison
            if newsolCost<popCost(i)
                popi(i,:) = newsol;
                popCost(i)=newsolCost;
                if popCost(i) < BestCost
                    BestCost=popCost(i);
                    BestSol = popi(i,:);
                    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                end
            end
        end
    end
   IterCurve(it)=min(popCost);    
   Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = IterCurve(it); 
end

Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);          
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                           
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve; 
end

