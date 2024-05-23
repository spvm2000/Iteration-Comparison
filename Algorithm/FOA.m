function Data_o = FOA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Ffun=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
dim=Data_i.dim;       Iterations=Data_i.maxIter;      Eval=Data_i.fobj;
 area_limit=Data_i.pop;  Life_time=15;   Transfer_rate=15;

Forest.P.Index_of_best=Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter); 
Forest.P.area_limit=area_limit;         % The limitation of the forest
Forest.P.Life_time=Life_time;           % The maximum allowed age of a tree 
Forest.P.Transfer_rate=Transfer_rate;   % The percentage of candidate population 
Forest.P.Dimension=dim;                 % The dimension of the problem domain
Forest.P.Llimit=Data_i.lb(1);            % The lower limit of the variables
Forest.P.Ulimit=Data_i.ub(1);            % The upper limit of the variables
Forest.P.MaxIterations=Iterations;      % Maximum number of iterations
Forest.P.dx=(abs(Forest.P.Ulimit)/5);   % dx is a small value used in local seeding. This value is not used in binary problems and in discrete problem, this value should be rounded.
if dim<5
    Forest.P.LSC=1; % Local seeding changes (1/5 of the dimension)
    Forest.P.GSC=1; % Global seeding changes
else
    Forest.P.LSC=floor((2*Forest.P.Dimension)/10); % 20 percent (not optimal) of the dimension used in local seeding
    Forest.P.GSC=floor((1*Forest.P.Dimension)/10); % 10 percent (not optimal) of the dimension used in global seeding   
end
% 初始化个体
Forest.T(:,:)=random('unif',Forest.P.Llimit,Forest.P.Ulimit,Forest.P.area_limit,Forest.P.Dimension+2);
for q=1:size(Forest.T,1)
     Forest.T(q,1:Forest.P.Dimension)=X(q,:);
     Forest.T(q,Forest.P.Dimension+1)=Eval(Forest.T(q,1:Forest.P.Dimension));
     Forest.T(q,Forest.P.Dimension+2)=0;
end
bestTree=Forest.T(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
for i=1:Forest.P.MaxIterations
    Forest=LocalSeeding(Forest,Eval);
    [Forest,candidate]=PopulationLimiting(Forest);
    Forest=GlobalSeeding(candidate,Forest,Eval);
    [Forest,bestTree]=UpdateBestTree(Forest,bestTree);
    IterCurve(i)=bestTree(Forest.P.Dimension+1);
    if bestTree(Forest.P.Dimension+1)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=bestTree(Forest.P.Dimension+1);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Forest.P.Index_of_best;
    end    
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=i;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function Forest=GlobalSeeding(candidate, Forest, Eval)
index=[];
New_Trees=[];
% choosing GSC variables randomly
for i=1:Forest.P.GSC
    index=[index, randi([1,Forest.P.Dimension])];
end 
% choose a percentage of candidate population
RevoSize=floor(((size(candidate,1))*Forest.P.Transfer_rate)/100);
SelectedTrees=candidate(1:RevoSize, :);
for w=1:size(SelectedTrees,1)    
    tempTree=SelectedTrees(w,:);
    for j=index
        tempTree(1,j)=random('unif',Forest.P.Llimit,Forest.P.Ulimit,1,1);
    end
    tempTree(1,Forest.P.Dimension+1)=Eval(tempTree(1,1:Forest.P.Dimension));
    tempTree(1,Forest.P.Dimension+2)=0;
    New_Trees=[New_Trees; tempTree(1,:)];
end
Forest.T=[Forest.T; New_Trees];
end

function Forest=LocalSeeding(Forest, Eval)
newtrees=[];
for q=1:size(Forest.T,1)
    if Forest.T(q,Forest.P.Dimension+2)==0 %Dimension+2 shows the Age of each tree
        current=Forest.T(q,:);
        % choosing LSC variables randomly
        move=[];
        for i=1:Forest.P.LSC
            move=[move, randi([1 Forest.P.Dimension])]; 
        end
        childs=[];
        for ww=move
            temp=current;
            add=random('unif',-Forest.P.dx ,Forest.P.dx ,1,1);
            temp(1,ww)=temp(1,ww)+add;
            if temp(1,ww)<Forest.P.Llimit
                temp(1,ww)=Forest.P.Llimit;
            elseif temp(1,ww)>Forest.P.Ulimit
                temp(1,ww)=Forest.P.Ulimit;
            end
            temp(1,Forest.P.Dimension+1)=Eval(temp(1,1:Forest.P.Dimension));
            temp(1,Forest.P.Dimension+2)=0;
            childs=[childs;temp];
        end
        
        newtrees=[newtrees;childs];
        childs=[];   
    end   
end
% Increasing the Age of all trees except for new trees
for u=1:size(Forest.T,1)
    Forest.T(u,Forest.P.Dimension+2)=Forest.T(u,Forest.P.Dimension+2)+1;
end
% adding the new trees with Age 0 to the Forest
Forest.T=[Forest.T; newtrees];
end

function [Forest,candidate]=PopulationLimiting(Forest)
temp=[];
candidate=[];
for i=1:size(Forest.T,1)
    if  Forest.T(i,Forest.P.Dimension+2)< Forest.P.Life_time
        temp=[temp; Forest.T(i,:)]; % trees with "Age" less than Life_time parameter
    else
        candidate=[candidate; Forest.T(i,:)]; % trees with "Age" bigger than Life_time parameter
    end
end
Forest.T=[];
Forest.T=temp;
[p,q]=sort(Forest.T(:,Forest.P.Dimension+1));% Dimension+1 shows the fitness
Forest.T=Forest.T(q,:);
if size(Forest.T,1)>Forest.P.area_limit % removing extra trees
    candidate=[candidate;Forest.T(Forest.P.area_limit+1:size(Forest.T,1),:)];
    Forest.T(Forest.P.area_limit+1:size(Forest.T,1), :)=[];
end
end

function [Forest, bestTree]=UpdateBestTree(Forest, bestTree)
Index_of_best=1;
for u=1:size(Forest.T,1)
    if Forest.T(u,Forest.P.Dimension+1)<=bestTree(1,Forest.P.Dimension+1)% Dimension+1 is the fitness value
        bestTree=Forest.T(u,:);
        Index_of_best=u;
    end
end
 if  Forest.T(Index_of_best,Forest.P.Dimension+1)<bestTree(1,Forest.P.Dimension+1) % update the best tree
     bestTree=Forest.T(Index_of_best,:);
     bestTree(1,Forest.P.Dimension+2)=0; % Dimension+2 shows the Age of each tree
     Forest.T(Index_of_best,:)=[];
     Forest.T=[bestTree;Forest.T];
 else % update the age of the best tree
     updated_best=Forest.T(Index_of_best,:);
     updated_best(1,Forest.P.Dimension+2)=0;
     Forest.T(Index_of_best,:)=[];
     Forest.T=[updated_best;Forest.T];
     bestTree=updated_best;
 end

end