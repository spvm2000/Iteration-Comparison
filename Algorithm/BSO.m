function Data_o = BSO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
popu = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve = zeros(1,Data_i.maxIter);
n_p = Data_i.pop;   % population size 100
n_d = Data_i.dim;   % dimension 10
funName=Data_i.fobj;    rang_l=Data_i.lb(1);    rang_r=Data_i.ub(1);    max_iteration = Data_i.maxIter;
n_c = 5;            % number of clusters 5      
prob_one_cluster = 0.8;     stepSize = ones(1,n_d);
popu_sorted  = rang_l + (rang_r - rang_l) * rand(n_p,n_d);      n_iteration = 1;
IterCurve(1)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
prob = zeros(n_c,1);        best = zeros(n_c,1); 
centers = rang_l + (rang_r - rang_l) * rand(n_c,n_d);  % initialize best individual in each cluster
centers_copy = rang_l + (rang_r - rang_l) * rand(n_c,n_d);
fitness_popu = 1000000*ones(n_p,1);  % store fitness value for each individual
fitness_popu_sorted = 1000000*ones(n_p,1);
indi_temp = zeros(1,n_d);
while n_iteration <= max_iteration
    cluster = kmeans(popu, n_c,'Distance','cityblock','Start',centers,'EmptyAction','singleton');
    fit_values = 100000000000000000000000000.0*ones(n_c,1); 
    number_in_cluster = zeros(n_c,1);
    for idx = 1:n_p
        number_in_cluster(cluster(idx,1),1)= number_in_cluster(cluster(idx,1),1) + 1;
        % find the best individual in each cluster
        if fit_values(cluster(idx,1),1) > fitness_popu(idx,1)  % minimization
            fit_values(cluster(idx,1),1) = fitness_popu(idx,1);
            best(cluster(idx,1),1) = idx;
        end
    end
    counter_cluster = zeros(n_c,1);
    acculate_num_cluster = zeros(n_c,1);
    for idx =2:n_c
        acculate_num_cluster(idx,1) = acculate_num_cluster((idx-1),1) + number_in_cluster((idx-1),1);
    end
    for idx = 1:n_p
        counter_cluster(cluster(idx,1),1) = counter_cluster(cluster(idx,1),1) + 1 ;
        temIdx = acculate_num_cluster(cluster(idx,1),1) +  counter_cluster(cluster(idx,1),1);
        popu_sorted(temIdx,:) = popu(idx,:);
        fitness_popu_sorted(temIdx,1) = fitness_popu(idx,1);
    end
    for idx = 1:n_c
        centers(idx,:) = popu(best(idx,1),:);        
    end
    
    centers_copy = centers;
    if (rand() < 0.2) %  select one cluster center to be replaced by a randomly generated center
        cenIdx = ceil(rand()*n_c);
        centers(cenIdx,:) = rang_l + (rang_r - rang_l) * rand(1,n_d);
    end 

    for idx = 1:n_c
        prob(idx,1) = number_in_cluster(idx,1)/n_p;
        if idx > 1
            prob(idx,1) = prob(idx,1) + prob(idx-1,1);
        end
    end

    for idx = 1:n_p
        r_1 = rand();  % probability for select one cluster to form new individual
        if r_1 < prob_one_cluster % select one cluster
            r = rand();
            for idj = 1:n_c
                if r < prob(idj,1)                      
                    if rand() < 0.4  % use the center
                       indi_temp(1,:) = centers(idj,:); 
                    else % use one randomly selected  cluster
                        indi_1 = acculate_num_cluster(idj,1) + ceil(rand() * number_in_cluster(idj,1));
                        indi_temp(1,:) = popu_sorted(indi_1,:);  
                    end
                    break
                end
            end
        else % select two clusters
            % pick two clusters 
            cluster_1 = ceil(rand() * n_c);
            indi_1 = acculate_num_cluster(cluster_1,1) + ceil(rand() * number_in_cluster(cluster_1,1));
            
            cluster_2 = ceil(rand() * n_c);
            indi_2 = acculate_num_cluster(cluster_2,1) + ceil(rand() * number_in_cluster(cluster_2,1));
            
            tem = rand();
            if rand() < 0.5 %use center
                indi_temp(1,:) = tem * centers(cluster_1,:) + (1-tem) * centers(cluster_2,:); 
            else   % use randomly selected individuals from each cluster            
                indi_temp(1,:) = tem * popu_sorted(indi_1,:) + (1-tem) * popu_sorted(indi_2,:); 
            end
        end         
        
        stepSize = logsig(((0.5*max_iteration - n_iteration)/20)) * rand(1,n_d);
        indi_temp(1,:) = indi_temp(1,:) + stepSize .* normrnd(0,1,1,n_d);
        % if better than the previous one, replace it
        fv = Data_i.fobj(indi_temp);
        if fv < fitness_popu(idx,1)  % better than the previous one, replace
            fitness_popu(idx,1) = fv;
            popu(idx,:) = indi_temp(1,:);
        end 
        if fitness_popu(n_p,1)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fitness_popu(n_p,1);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=n_p;
        end
    end

    for idx = 1:n_c
        popu(best(idx,1),:) = centers_copy(idx,:);  
        fitness_popu(best(idx,1),1) = fit_values(idx,1);
    end   
    % record the best fitness in each iteration
    IterCurve(n_iteration) = min(fit_values);
    n_iteration=n_iteration+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=n_iteration;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end