function Data_o = SALPSA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                

%% Problem Definition
Max_iter = Data_i.maxIter;            
N = Data_i.pop;                       
dim = Data_i.dim;                     
SalpPositions = Data_i.X;             
SalpFitness = Data_i.F_value;         
% if size(ub,1)==1
%     ub=ones(dim,1)*ub;
%     lb=ones(dim,1)*lb;
% end

Convergence_curve = zeros(1,Max_iter);

[sorted_salps_fitness,sorted_indexes]=min(SalpFitness);
FoodPosition=SalpPositions(sorted_indexes,:);
FoodFitness=sorted_salps_fitness;
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(SalpFitness);

%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while l<Max_iter+1
    
    c1 = 2*exp(-(4*l/Max_iter)^2); % Eq. (3.2) in the paper
    
    for i=1:size(SalpPositions,1)
        
        SalpPositions= SalpPositions';
        
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    SalpPositions(j,i)=FoodPosition(j)+c1*((Data_i.ub(j)-Data_i.lb(j))*c2+Data_i.lb(j));
                else
                    SalpPositions(j,i)=FoodPosition(j)-c1*((Data_i.ub(j)-Data_i.lb(j))*c2+Data_i.lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>N/2 && i<N+1
            point1=SalpPositions(:,i-1);
            point2=SalpPositions(:,i);
            
            SalpPositions(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        SalpPositions= SalpPositions';
    end
    
    for i=1:size(SalpPositions,1)
        
        SalpPositions=BoundaryCheck(SalpPositions,Data_i.ub,Data_i.lb,Data_i.dim);
        
        SalpFitness(1,i)=Data_i.fobj(SalpPositions(i,:));
        
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=SalpPositions(i,:);
            FoodFitness=SalpFitness(1,i);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    Convergence_curve(l)=FoodFitness;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = FoodFitness;
    l = l + 1;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end


