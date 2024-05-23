function Data_o = AVOA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
Curve=zeros(1,Data_i.maxIter);
%% Problem Definition
max_iter = Data_i.maxIter;      
pop_size = Data_i.pop;            
variables_no = Data_i.dim;              
lower_bound = Data_i.lb;
upper_bound = Data_i.ub;
X = Data_i.X;               
fit = Data_i.F_value;         
    % initialize Best_vulture1, Best_vulture2
    Best_vulture1_X=zeros(1,variables_no);
    Best_vulture1_F=inf;
    Best_vulture2_X=zeros(1,variables_no);
    Best_vulture2_F=inf;

    [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);
    %%  Controlling parameter
    p1=0.6;
    p2=0.4;
    p3=0.6;
    alpha=0.8;
    betha=0.2;
    gamma=2.5;

    %%Main loop
    current_iter=0; % Loop counter

    while current_iter < max_iter
        

        for i=1:pop_size

            if current_iter == 0
                current_vulture_F=fit(i);
                current_vulture_X=X(i,:);
            end
            % Update the first best two vultures if needed
            if current_vulture_F<Best_vulture1_F
                Best_vulture1_F=current_vulture_F; % Update the first best bulture
                Best_vulture1_X=current_vulture_X;
            end

            if current_vulture_F>Best_vulture1_F && current_vulture_F<Best_vulture2_F
                Best_vulture2_F=current_vulture_F; % Update the second best bulture
                Best_vulture2_X=current_vulture_X;
            end
        end

        a=unifrnd(-2,2,1,1)*((sin((pi/2)*(current_iter/max_iter))^gamma)+cos((pi/2)*(current_iter/max_iter))-1);
        P1=(2*rand+1)*(1-(current_iter/max_iter))+a;

        % Update the location
        for i=1:size(X,1)
            current_vulture_X = X(i,:);  % pick the current vulture back to the population
            F=P1*(2*rand()-1);  

            random_vulture_X=random_select(Best_vulture1_X,Best_vulture2_X,alpha,betha);
            
            if abs(F) >= 1 % Exploration:
                current_vulture_X = exploration1(current_vulture_X, random_vulture_X, F, p1, upper_bound, lower_bound);
            elseif abs(F) < 1 % Exploitation:
                current_vulture_X = exploitation2(current_vulture_X, Best_vulture1_X, Best_vulture2_X, random_vulture_X, F, p2, p3, variables_no);
            end

            X(i,:) = current_vulture_X; % place the current vulture back into the population
        end
        
        for i=1:pop_size
            X = BoundaryCheck(X,Data_i.ub,Data_i.lb,Data_i.dim);
            % Calculate the fitness of the population
            current_vulture_F=Data_i.fobj(X(i,:));
            % Update the last bestfitness
            if current_vulture_F<Best_vulture1_F
                Best_vulture1_F=current_vulture_F; 
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end 

        current_iter=current_iter+1;
        Curve(current_iter)=Best_vulture1_F;
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Best_vulture1_F;
   end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=current_iter;                               
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Curve;
end

function [random_vulture_X]=random_select(Best_vulture1_X,Best_vulture2_X,alpha,betha)

    probabilities=[alpha, betha ];
    
    if (rouletteWheelSelection( probabilities ) == 1)
            random_vulture_X=Best_vulture1_X;
    else
            random_vulture_X=Best_vulture2_X;
    end

end

function [index] = rouletteWheelSelection(x)

    index=find(rand() <= cumsum(x) ,1,'first');

end

function [current_vulture_X] = exploration1(current_vulture_X, random_vulture_X, F, p1, upper_bound, lower_bound)

    if rand<p1
        current_vulture_X=random_vulture_X-(abs((2*rand)*random_vulture_X-current_vulture_X))*F;
    else
        current_vulture_X=(random_vulture_X-(F)+rand()*((upper_bound-lower_bound)*rand+lower_bound));
    end
    
end

function [current_vulture_X] = exploitation2(current_vulture_X, Best_vulture1_X, Best_vulture2_X,random_vulture_X, F, p2, p3, variables_no)
                                                                      
% phase 1
    if  abs(F)<0.5
        if rand<p2
            A=Best_vulture1_X-((Best_vulture1_X.*current_vulture_X)./(Best_vulture1_X-current_vulture_X.^2))*F;
            B=Best_vulture2_X-((Best_vulture2_X.*current_vulture_X)./(Best_vulture2_X-current_vulture_X.^2))*F;
            current_vulture_X=(A+B)/2;
        else
            current_vulture_X=random_vulture_X-abs(random_vulture_X-current_vulture_X)*F.*levyFlight(variables_no);
        end
    end
    % phase 2
    if  abs(F)>=0.5
        if rand<p3
            current_vulture_X=(abs((2*rand)*random_vulture_X-current_vulture_X))*(F+rand)-(random_vulture_X-current_vulture_X);
        else
            s1=random_vulture_X.* (rand()*current_vulture_X/(2*pi)).*cos(current_vulture_X);
            s2=random_vulture_X.* (rand()*current_vulture_X/(2*pi)).*sin(current_vulture_X);
            current_vulture_X=random_vulture_X-(s1+s2);
        end
    end
end

function [ o ]=levyFlight(d)
  
    beta=3/2;

    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;
    v=randn(1,d);
    step=u./abs(v).^(1/beta);

    o=step;

end


