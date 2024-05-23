function Data_o = YSGA(Data_i,Data_o)
    ti=clock;                             
    [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
% Initialization
    n = Data_i.pop;  
    maxGeneration = Data_i.maxIter;
    dim = Data_i.dim;
    groupN = 4;
    range = [Data_i.lb(1),Data_i.ub(1)];
    f = Data_i.fobj;
    globalBest.Sol = Data_i.X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);         
    globalBest.Fitness = Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);                      
    globalBest.Group = Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter);                          
    FitHist=zeros(1,Data_i.maxIter);

    fish.Pos = cell(1,groupN);
    fish.Fitness = cell(1,groupN);
    
    fish.groupBest = Data_i.X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
    FitHist = zeros(1,maxGeneration);               
    flags = zeros(1,groupN);
    
    [fish] = ysga_init(f,fish,n,range,groupN,dim,Data_i.X);
% Main iterative process
    for i = 1:maxGeneration
        % For each group of fish
        for j=1:groupN
           zn = cell2mat(fish.Fitness(1,j));
           alpha = -1+i*((-1)/maxGeneration);
           % Chaser fish movement
           [fish] = ysga_hunt(fish,j,range,f,dim,i,maxGeneration,globalBest);
           % Blocker fish movement
           [fish] = ysga_block(fish,j,range,alpha,dim);  
           xn2 = cell2mat(fish.Pos(1,j));
           z2 = zeros(size(xn2(:,1)));
           ni = size(z2);
           for k = 1:ni(1) 
               z2(k) = f(xn2(k,:)); 
           end
           fish.Fitness{1,j} = z2;
           [zn2,Ind2] = sort(cell2mat(fish.Fitness(1,j)), 'ascend');
           for k=1:dim
                xn2(:,k) = xn2(Ind2,k);
           end
           fish.Pos{1,j} = xn2;
           fish.Fitness{1,j} = zn2;
           %Update best found global solution
           if zn2(1) < globalBest.Fitness
               globalBest.Fitness = zn2(1);
               globalBest.Sol = xn2(1,:);
               globalBest.Group = j;
           %Update flag for change of zone    
           elseif zn2(1) > globalBest.Fitness
               flags(j) = flags(j) + 1; 
           end
           %Update best found solution of the current group
           if zn2(1) < zn(1)
               fish.groupBest = zn2(1);
           elseif zn2(1)== zn(1)
               flags(j) = flags(j) + 1; 
           end
        end
        
        for j=1:groupN
            if flags(j) >= 10 && j ~= globalBest.Group
                %Start change zone behaviour
                ysga_changezone(fish,j,range,dim,globalBest);
                flags(j) = 0;
            end
        end
       FitHist(i) = globalBest.Fitness;             
    end  
    bestFit = globalBest.Fitness;
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=bestFit;             
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=globalBest.Group;        

Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=i;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=FitHist;
end

function [fish] = ysga_init(f,fish,n,range,groupN,dim,X)
    xo = X';
%     for i=1:dim
%         xo(i,:) = rand(1,n)*(range(2)-range(1))+range(1);
%     end
    [idx,~] = kmeans(transpose(xo),groupN,'Distance','sqeuclidean','Replicates',5);
    for j=1:groupN
        fish.Pos{1,j} = transpose(xo(:,idx==j));
        xo2 = cell2mat(fish.Pos(1,j));
        zo2 = zeros(size(xo2(:,1)));
        ni = size(zo2);
        for k = 1:ni(1) 
            zo2(k) = f(xo2(k,:)); 
        end
        fish.Fitness{1,j} = zo2;
        [zn,Ind]=sort(cell2mat(fish.Fitness(1,j)), 'ascend');
        xn = zeros(size(xo2));
        for k=1:dim
             xn(:,k) = xo2(Ind,k);
        end
        fish.Pos{1,j} = xn;
        fish.Fitness{1,j} = zn;
    end
end

function [fish] = ysga_hunt(fish,j,range,f,dim,iter,maxGeneration,globalBest)
    xo = transpose(cell2mat(fish.Pos(1,j)));
    zo = transpose(cell2mat(fish.Fitness(1,j)));
    xn = zeros(dim,10);
    zn = zeros(1,10);
    beta = 1.99 + (.001*iter/(maxGeneration/10));
    if beta > 2
        beta = 2;
    end
    levySteps = levyStep(10,dim,beta);
    if (transpose(xo(:,1))-globalBest.Sol) ~= 0
        for i = 1:10
             levySteps(i,:) = 1*levySteps(i,:).*(transpose(xo(:,1))-globalBest.Sol);
             xn(:,i) = xo(:,1) + transpose(levySteps(i,:));  
        end
    else
        for i = 1:dim
            xn(i,:) = xo(i,1) + transpose(levySteps(:,i));  
        end
    end
         
    for i = 1:10
       zn(i) = f(transpose(xn(:,i)));
    end
    [Fitness,Ind]=sort(zn, 'ascend');
    xn = xn(:,Ind);
    zn = Fitness;
    if zn(1) < zo(1)
        for i = 1:dim
            xo(i,1) = xn(i,1);
        end
        zo(1) = zn(1);
    end
    [xo] = findrange(xo,range,dim);
    fish.Pos{1,j} = transpose(xo);
    fish.Fitness{1,j} = zo;
end

function [fish] = ysga_block(fish,j,range,alpha,dim)
    xo = transpose(cell2mat(fish.Pos(1,j)));
    xn = xo;
    ni = size(xn(1,:),2);
    b = 1;
    for i = 2:ni
        t = (alpha-1)*rand+1;
        rx = zeros(1,dim);
        for k = 1:dim
           r = (rand*2-1);
           rx(k) = abs(xn(k,1) * r - xn(k,i)); 
           xn(k,i) = rx(k) * exp(b.*t) .* cos(t.*2*pi) + xn(k,1);
        end
    end
    [xn] = findrange(xn,range,dim);
    fish.Pos{1,j} = transpose(xn);
end

function [fish] = ysga_changezone(fish,j,range,dim,globalBest)
    xo = transpose(cell2mat(fish.Pos(1,j)));
    for i=1:dim
         xo(i,:) =  (globalBest.Sol(i) + xo(i,:))/2;
    end
    [xo] = findrange(xo,range,dim);
    fish.Pos{1,j} = transpose(xo);
end

function [xo] = findrange(xo,range,dim)
    for i=1:dim
        for j=1:length(xo(1,:))
            if xo(i,j)<=range(1)
                xo(i,j)=range(1);
            end
            if xo(i,j)>=range(2)
                xo(i,j)=range(2);
            end
        end
    end
end

function [z] = levyStep(n,m,beta)
% This function implements Levy's flight. 

% For more information see 
%'Multiobjective cuckoo search for design optimization Xin-She Yang, Suash Deb'. 

% Coded by Hemanth Manjunatha on Nov 13 2015.

% Input parameters
% n     -> Number of steps 
% m     -> Number of Dimensions 
% beta  -> Power law index  % Note: 1 < beta < 2

% Output 
% z     -> 'n' levy steps in 'm' dimension

    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator 
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator

    sigma_u = (num/den)^(1/beta);% Standard deviation

    u = random('Normal',0,sigma_u^2,n,m); 
    
    v = random('Normal',0,1,n,m);

    z = u./(abs(v).^(1/beta));
end