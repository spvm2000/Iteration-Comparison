function Data_o = HGS(Data_i,Data_o)
rand_num=[];                           
ti=clock;                               
X = Data_i.X;                        
AllFitness=Data_i.F_value;                   
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
N=Data_i.pop;   Function_name=Data_i.fobj;     dimSize=Data_i.dim;
bestPositions=zeros(1,Data_i.dim);    Max_iter=Data_i.maxIter;
tempPosition=zeros(N,Data_i.dim);
Destination_fitness=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);  
Worstest_fitness=-inf;
VC1 = ones(N,1);
weight3 = ones(N,Data_i.dim);
weight4 = ones(N,Data_i.dim);

IterCurve=zeros(1,Max_iter);
it=1; 
hungry = zeros(1,size(X,1));
count=0;

while  it <= Data_i.maxIter
    VC2 = 0.03; %The variable of variation control 
    sumHungry = 0;%record the sum of each hungry 
    %sort the fitness
    for i=1:size(X,1)
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>Data_i.ub;
        Flag4lb=X(i,:)<Data_i.lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+Data_i.ub.*Flag4ub+Data_i.lb.*Flag4lb;
        AllFitness(i) = Data_i.fobj(X(i,:));
    end
    [AllFitnessSorted,IndexSorted] = sort(AllFitness);
    bestFitness = AllFitnessSorted(1);
    worstFitness = AllFitnessSorted(size(X,1));
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(IndexSorted(1),:);
        Destination_fitness = bestFitness;
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Destination_fitness;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=IndexSorted(1);
        count=0;
    end
    if worstFitness > Worstest_fitness
        Worstest_fitness = worstFitness;
    end
    
    for i = 1:size(X,1)
         %calculate the variation control of all positions
         VC1(i) = sech(abs(AllFitness(i)-Destination_fitness));    
         %calculate the hungry of each position
        if Destination_fitness == AllFitness(i)
            hungry(1,i) = 0;
            count = count+1;
            tempPosition(count,:)=X(i,:);
        else
            temprand = rand();
            c = (AllFitness(i)-Destination_fitness)/(Worstest_fitness-Destination_fitness)*temprand*2*(Data_i.ub-Data_i.lb);
            if c<100
                b=100*(1+temprand);
            else
                b=c;
            end   
            hungry(1,i) = hungry(1,i)+ max(b); 
            sumHungry = sumHungry + hungry(1,i);
        end
    end 
    
    %calculate the hungry weight of each position
    for i=1:size(X,1)
        for j=2:size(X,2)
            weight3(i,j) = (1-exp(-abs(hungry(1,i)-sumHungry)))*rand()*2;
            if rand()<VC2
                weight4(i,j) = hungry(1,i)*size(X,1)/sumHungry*rand();
            else
                weight4(i,j) = 1;
            end
        end
    end
    % Update the Position of search agents
    shrink=2*(1-it/Max_iter); % a decreases linearly fron 2 to 0
    for i=1:size(X,1)
        if rand<VC2
            X(i,:) = X(i,j)*(1+randn(1));
        else
            if count==0  count=1;  end
            A = randi([1,count]);
            for j=1:size(X,2)
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    X(i,j) = weight4(i,j)*tempPosition(A,j)+vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));
                else
                    X(i,j) = weight4(i,j)*tempPosition(A,j)-vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));
                end
            end
        end
    end
    IterCurve(it)=Destination_fitness;
    it=it+1;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end