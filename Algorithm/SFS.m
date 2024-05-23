function Data_o = SFS(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
point = Data_i.X;                        
FirstFit=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
S.Start_Point = Data_i.pop;            
S.Maximum_Generation  = Data_i.maxIter;     
S.Maximum_Diffusion = 2;
S.Function_Name=Data_i.fobj;
S.Walk = 1; % *Important
S.plot = 1;     S.ShowResult = 0; 
%S.Walk = 1 ----> SFS uses the first Gaussian walk(usually SIMPLE Problems)
%S.Walk = 0 ----> SFS uses the second Gaussian walk(usually HARD Problems)
S.Ndim = Data_i.dim;    S.Lband = Data_i.lb;    S.Uband = Data_i.ub;
[Sorted_FitVector, Indecis] = sort(FirstFit);       point = point(Indecis,:);
BestPoint = point(1, :);        F = Sorted_FitVector(1);
for G = 1:S.Maximum_Generation
    New_Point = [];
    FitVector = [];
    for i = 1 : size(point,1)
        %creating new points based on diffusion process
        [NP, fit] = Diffusion_Process(point(i,:),S,G,BestPoint);
        New_Point = [New_Point;NP];
        FitVector = [FitVector,fit];
    end

    BestFit = min(FitVector);
    BestPoint = New_Point(find(FitVector == BestFit),:);
    S.Start_Point = size(New_Point,1);
    fit = FitVector';
    [sortVal, sortIndex] = sort(fit);

    for i=1:1:S.Start_Point     
        Pa(sortIndex(i)) = (S.Start_Point - i + 1) / S.Start_Point;
    end

    RandVec1 = randperm(S.Start_Point);
    RandVec2 = randperm(S.Start_Point);

    for i = 1 : S.Start_Point
        for j = 1 : size(New_Point,2)
            if rand > Pa(i)
                P(i,j) = New_Point(RandVec1(i),j) - rand*(New_Point(RandVec2(i),j) - ...
                    New_Point(i,j));
            else
                P(i,j)= New_Point(i,j);
            end
        end
    end
    P = BoundaryCheck(P,S.Uband,S.Lband,S.Ndim);%for checking bounds
    Fit_FirstProcess = [];

    for i = 1 : S.Start_Point
        Fit_FirstProcess = [Fit_FirstProcess;feval(S.Function_Name,P(i,:))];
    end
    for i=1:S.Start_Point
        if Fit_FirstProcess(i,:)<=fit(i,:)
            New_Point(i,:)=P(i,:);
            fit(i,:)=Fit_FirstProcess(i,:);
        end
    end
    FitVector = fit;
    [SortedFit,SortedIndex] = sort(FitVector);
    New_Point = New_Point(SortedIndex,:);
    BestPoint = New_Point(1,:);%first point is the best
    F = [F;FitVector(1,1)];
    F = sort(F);
    if S.ShowResult == 1
        fprintf('Iteration: %i', G);
        fprintf(',    Best result: %e \n', F(1,1));
    end
    fbest = FitVector(1,:);
    if fbest <= F(1,:)
        pbest = New_Point(1,:);
        fbest = F(1,:);
    end
    point = New_Point;
    Pa = sort(SortedIndex/S.Start_Point, 'descend');
    for i = 1 : S.Start_Point
       if rand > Pa(i)
           %selecting two different points in the group
           R1 = ceil(rand*size(point,1));
           R2 = ceil(rand*size(point,1));
            while R1 == R2
                R2 = ceil(rand*size(point,1));
            end
            
            if rand < .5
                ReplacePoint = point(i,:) - ...
                    rand * (point(R2,:) - BestPoint);
                ReplacePoint = BoundaryCheck(ReplacePoint,S.Uband,S.Lband,S.Ndim);
            else
                ReplacePoint = point(i,:) + ...
                    rand * (point(R2,:) - point(R1,:));
                ReplacePoint = BoundaryCheck(ReplacePoint,S.Uband,S.Lband,S.Ndim);
            end
            
            if feval(S.Function_Name, ReplacePoint) < ...
                    feval(S.Function_Name, point(i,:))
                point(i,:) = ReplacePoint;
            end
       end
    end
    fbest = min(F);
    [value,index]=sort(F);
    if value(1)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=value(1);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
    end
    IterCurve(G)=value(1);
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=G;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [createPoint, fitness] = Diffusion_Process(Point,S,g,BestPoint)
    %calculating the maximum diffusion for each point
    NumDiffiusion = S.Maximum_Diffusion;
    New_Point = Point;
    
    %Diffiusing Part*******************************************************
    for i = 1 : NumDiffiusion
        %consider which walks should be selected.
        if rand < S.Walk 
            GeneratePoint = normrnd(BestPoint, (log(g)/g)*(abs((Point - BestPoint))), [1 size(Point,2)]) + ...
                (randn*BestPoint - randn*Point);
        else
            GeneratePoint = normrnd(Point, (log(g)/g)*(abs((Point - BestPoint))),...
                [1 size(Point,2)]);
        end

        New_Point = [New_Point;GeneratePoint];

    end
    %check bounds of New Point
    New_Point = BoundaryCheck(New_Point,S.Uband,S.Lband,S.Ndim);
    %sorting fitness
    fitness = [];
    for i = 1 : size(New_Point,1)
        fitness = [fitness;feval(S.Function_Name,New_Point(i,:))];
    end

    [fit_value,fit_index] = sort(fitness);
    fitness = fit_value(1,1);
    New_Point = New_Point(fit_index,:);
    createPoint = New_Point(1,:);
    %======================================================================
end