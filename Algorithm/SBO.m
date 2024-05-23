function Data_o = SBO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 
%% Problem Definition
MaxIt = Data_i.maxIter;        
nPop = Data_i.pop;         
numbervar = Data_i.dim;               
pop.Position = Data_i.X;                   
pop.Cost = Data_i.F_value;     
upperbound=Data_i.ub;
lowerbound=Data_i.lb;

% VarSize=[1 numbervar];  

%% SBO Parameters
%Greatest step size
alpha=0.94;
%Mutation probability
pMutation=0.05;
%Percent of the difference between the upper and lower limit
Z=0.02;
%proportion of space width
sigma=Z*(max(upperbound)-min(lowerbound));

%% Initialization
[bestfit,index]=min([pop.Cost]);
elite=pop.Position(index,:);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=bestfit;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
BestCost=zeros(1,Data_i.maxIter);
%% SBO Main Loop
for it=1:MaxIt
    newpop=pop;
    
    %Calculating the Fitness of each bower
    F=zeros(nPop,1);
    for i=1:nPop
        if pop.Cost(i)>=0
            F(i)=1/(1+pop.Cost(i));
        else
            F(i)=1+abs(pop.Cost(i));
        end
    end
    
    %Calculating the probability of each bower
    P=F/sum(F);
    
    %changes at any bower
    for i=1:nPop
        for k=1:numbervar
                
                % Select target bower                 
                j=RouletteWheelSelection(P);
                
                % Calculating Step Size
                lambda=alpha/(1+P(j));
                
                newpop.Position(i,k)=pop.Position(i,k) ...
                    +lambda*(((pop.Position(j,k)+elite(k))/2)-pop.Position(i,k));
                
                % Mutation
            if rand<=pMutation
                newpop.Position(i,k)=newpop.Position(i,k)+(sigma*randn);
            end
                
        end  
        % Evaluation
        newpop.Cost(i)=Data_i.fobj(newpop.Position(i,:));
    end 
          
     maxpop=[pop
         newpop]; %#ok
  
    % Sort Population
    [tempop, SortOrder]=sort([maxpop.Cost]);
    pop.Cost=tempop(1:nPop);
    tempop_Position=[pop.Position;newpop.Position];
    for pos=1:nPop
        sortpos=SortOrder(pos);
        Position(pos,:)=tempop_Position(sortpos,:);
    end
    pop.Position=Position;
    % Update Best Solution Ever Found
    elite=pop.Position(1,:);
    BestCost(it)=pop.Cost(1);
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=BestCost;
end

function f=RouletteWheelSelection(P)

    r=rand;
    C=cumsum(P);
    f=find(r<=C,1,'first');

end