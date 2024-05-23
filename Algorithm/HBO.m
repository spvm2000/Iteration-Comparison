function Data_o = HBO(Data_i,Data_o)
rand_num=[];                        
ti=clock;                            
Solutions = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
searchAgents=Data_i.pop;        
cycles = floor(Data_i.maxIter/25);
degree = 3;
treeHeight = ceil((log10(searchAgents * degree - searchAgents + 1)/log10(degree))); %Starting from 1
fevals = 0;
Leader_pos=zeros(1,Data_i.dim);
Leader_score=inf;
fitnessHeap = zeros(searchAgents, 2) + inf;
for c = 1:searchAgents
    fevals = fevals +1;
    fitnessHeap(c, 1) = fitness(c);
    fitnessHeap(c, 2) = c;
    %Heapifying
    t = c;
    while t > 1
        parentInd = floor((t+1)/degree);
        if fitnessHeap(t, 1) >= fitnessHeap(parentInd,1)
            break;
        else
            tempFitness = fitnessHeap(t,:);
            fitnessHeap(t,:) = fitnessHeap(parentInd,:);
            fitnessHeap(parentInd,:) = tempFitness;
        end
        t = parentInd;
    end
    
    if fitness(c) <= Leader_score
        Leader_score = fitness(c);
        Leader_pos = Solutions(c,:);
    end
end

colleaguesLimits = colleaguesLimitsGenerator(degree,searchAgents);
itPerCycle =Data_i.maxIter/cycles;
qtrCycle = itPerCycle / 4;
for it=1:Data_i.maxIter
    gamma = (mod(it, itPerCycle)) / qtrCycle;
    gamma = abs(2-gamma);
    for c = Data_i.pop:-1:2
        if c == 1 %Dealing with root
            continue;
        else
            parentInd = floor((c+1)/degree);
            curSol = Solutions(fitnessHeap(c,2), :); %Sol to be updated
            parentSol = Solutions(fitnessHeap(parentInd,2), :); %Sol to be updated with reference to
            
            if colleaguesLimits(c,2) > searchAgents
                colleaguesLimits(c,2) = searchAgents;
            end
            colleagueInd = c;
            while colleagueInd == c
                colleagueInd = randi([colleaguesLimits(c,1) colleaguesLimits(c,2)]);
            end
            colleagueSol = Solutions(fitnessHeap(colleagueInd,2), :); %Sol to be updated with reference to
            
            %Position Updating
            for j = 1:Data_i.dim
                p1 =  (1 - it/(Data_i.maxIter));
                p2 = p1+(1- p1)/2;
                r = rand();
                rn = (2*rand()-1);
                
                if r < p1      %To skip any dim to update
                    continue;
                elseif  r < p2
                    D = abs(parentSol(j) - curSol(j));
                    curSol(1, j) = parentSol(j) + rn * gamma * D;
                else
                    if fitnessHeap(colleagueInd,1) < fitnessHeap(c,1)
                        D = abs(colleagueSol(j) - curSol(j));
                        curSol(1, j) = colleagueSol(j) + rn * gamma * D;
                    else
                        D = abs(colleagueSol(j) - curSol(j));
                        curSol(1, j) = curSol(j) + rn * gamma * D;
                    end
                end
            end
        end
        curSol(1, :)=BoundaryCheck(curSol(1, :),Data_i.ub,Data_i.lb,Data_i.dim);
        newFitness = Data_i.fobj(curSol);
        fevals = fevals +1;
        if newFitness < fitnessHeap(c,1)
            fitnessHeap(c,1) = newFitness;
            Solutions(fitnessHeap(c,2), :) = curSol;
        end
        if newFitness < Leader_score
            Leader_score = newFitness;
            Leader_pos = curSol;
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Leader_score
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Leader_score;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=c;
            end
        end
        t = c;
        while t > 1
            parentInd = floor((t+1)/degree);
            if fitnessHeap(t, 1) >= fitnessHeap(parentInd,1)
                break;
            else
                tempFitness = fitnessHeap(t,:);
                fitnessHeap(t,:) = fitnessHeap(parentInd,:);
                fitnessHeap(parentInd,:) = tempFitness;
            end
            t = parentInd;
        end
    end
    IterCurve(it)=Leader_score;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [colleaguesLimits]= colleaguesLimitsGenerator(degree,searchAgents)
    colleaguesLimits = zeros(searchAgents,2);
    for c = searchAgents: -1 : 1
        hi = ceil((log10(c * degree - c + 1)/log10(degree))) - 1;
        lowerLim = ((degree * degree^(hi-1) - 1)/(degree-1) + 1);
        upperLim = (degree * degree^hi - 1)/(degree-1);
        colleaguesLimits(c,1) = lowerLim;
        colleaguesLimits(c,2) = upperLim;
    end
end