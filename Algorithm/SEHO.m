function Data_o = SEHO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 
%% Problem Definition
itern = Data_i.maxIter;      
N = Data_i.pop;            
dims=Data_i.dim;               
xd = min(Data_i.lb);
xu = max(Data_i.ub);
    %% STEP 1: Initialize the animal population "A"
            A = Data_i.X;                                                       % Initialize "N" random positions within the bounds [xd,xu]
    %% STEP 2: Divide "A" in two groups: herd members "H" and predators "P"
        % 2.1. Define the number of herd members and predators (N_h and N_p)
            preys_rate = [0.7,0.9];                                             % Lower and Upper herd members' Rate (recommended values are 0.7 and 0.9)
            rate = preys_rate(1)+(preys_rate(2) - preys_rate(1))*rand;          % Randomly define herd members' between lower and upper values
            N_h = round(N*rate);                                                % Herd members' population size
            N_p = N - N_h;                                                      % Predators' population size   
        % 2.2. Assign N_h individuals from "A" as herd members and 
        % calculate their fintess values
            H = A(1:N_h,:);                                                     % Define the first "N_h" members in "A" as herd members
            f_h = zeros(N_h,1);                                                 % Initialize herd members fitness values as zeros
            for ind = 1:N_h                                                     % FOR each member in "H"
                f_h(ind) = Data_i.fobj(H(ind,:));                                         %   Evaluate Fitness Function "f(.)"
            end                                                                 % END
        % 2.3. Assign N_p individuals from "A" as predators and 
        % calculate their fitness values
            P = A(N_h+1:N,:);                                                   % Define the remaining "N_p" members in "A" as herd members
            f_p = zeros(N_p,1);                                                 % Initialize herd members fitness values as zeros
            for ind = 1:N_p                                                     % FOR each member in "P"
                f_p(ind) = Data_i.fobj(P(ind,:));                                             % Evaluate Fitness Function "f(.)"
            end                                                                 % END
        % 2.4. Initialize Global Memory
            fit = [f_h;f_p];                                                    % Group fitness values from "H" and "P" as a single vector
            [~,b_idx] = min(fit);                                               % Find the index "b_idx" of the best solution
            xbest = A(b_idx,:);                                                 % Initial Best Position 
            bfit = fit(b_idx);                                                  % Initial Best Fitness Value
            [~,w_idx] = max(fit);                                               % Find the index "w_idx" of the worst solution
            wfit = fit(w_idx);
        % 2.5. Initialize Historical Global Memory
            xbesth = zeros(itern,dims);                                         % Initialize HISTORICAL GLOBAL BEST vector
            bfith = zeros(1,itern);                                             % Initialize HISTORICAL BEST FITNESS vector        
    %% STEP 3: Calculate Survival Values for each member in "H" and "P"
            [SV_h, SV_p] = SV(bfit, wfit, f_h, f_p);                            % Calculate Survival Values "SV" for each member in "H" and "P"
    %% STEP EX1: Build Function Contour (For Visualization Purposes)
%             [X,Y,Z] = contour_plot(Data_i.fobj,xd,xu);                                    % Create Function "f" Contour Data (visualization purposes)
%% B. INITIALIZE ITERATIONS
    for it=1:itern                                                              % FOR each ITERATION "it"
    %% STEP 4: Move all members in "H" by appying herd movement operators
        H = HerdMove(H,SV_h,P,SV_p,xbest,dims,xd,xu);                               % Herd Movement Operators
    %% STEP 5: Move all members in "P" by applying predators movement operators
        P = PredMove(P,H,SV_h,dims,xd,xu);                                          % Predator Movement Operators
    %% STEP 6: Re-calculate Survival Values for each member in "H" and "P"
        for i=1:N_h                                                                 % FOR each PREY "i"
            f_h(i)=Data_i.fobj(H(i,:));                                                           % Evaluate Cost Function "f" for PREY "i"
        end                                                                         % END
        for i=1:N_p                                                                 % FOR each PREDATOR "i"
            f_p(i)=Data_i.fobj(P(i,:));                                                           % Evaluate Cost Function "f" for PREDATOR "i"
        end                                                                         % END
        [SV_h, SV_p] = SV(bfit, wfit, f_h, f_p);                                    % Calculate Survival Values "SV" for each member in "H" and "P"
    %% STEP 7: Perform Predation Phase
        [L_idx, K_idx] = Predation(P,H,SV_p,SV_h,dims);                             % Determine ALIVE and KILLED herd members indexes
    %% STEP 8: Perform Herd Restoration Phase
        if any(K_idx)                                                               % IF there is any killed herd member
            [H, f_h] = Restoration(H,f_h,SV_h,L_idx,K_idx,dims,Data_i.fobj);                      % Replace KILLED herd members with new Solutions
            [SV_h, SV_p] = SV(bfit, wfit, f_h, f_p);                                    % Calculate Survival Values "SV" for each member in "H" and "P"
        end                                                                         % END
    %% STEP EX2: Update Global Memory
        A = [H; P];                                                                 % Combine PREYS and PREDATORS into a single set of solutions
        fit = [f_h; f_p];                                                           % Combine PREYS and PREDATOS Fitness Values into a single set 
        [~,b_idx] = min(fit);                                                       % Find the index "Ibe" of the BEST SOLUTION for Iteration "It"
        gbest2 = A(b_idx,:);                                                        % Get BEST POSITION "gbest2" for Iteration "It"
        bfit2 = fit(b_idx);                                                         % Get BEST FITNESS "bfit2" for Iteration "It"       
        if bfit2 <= bfit                                                            % IF the new BEST Fitness is better than that on the Global Memory
            bfit = bfit2;                                                           % Update current BEST FITNESS on the Global Memory with "bfit2"
            xbest = gbest2;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=b_idx;       % Update current BEST POSITION on the Global Memory with "gbest2"
        end                                                                         % END
        [~,w_idx] = max(fit);                                                       % Find the INDEX "Ibe" of the Initial WORST SOLUTION
        wfit2 = fit(w_idx);                                                         % Get BEST FITNESS "wfit2" for Iteration "It"
        if wfit2 >= wfit                                                            % IF the new WORST Fitness is better than that on the Global Memory
            wfit = wfit2;                                                           % Update current WORST FITNESS on the Global Memory with "wfit2"
        end                                                                         % END
        xbesth(it,:) = xbest;                                                       % Update Historical BEST POSITION for Iteration "It"
        bfith(it) = bfit;                                                           % Update Historical BEST FITNESS for Iteration "It"
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = bfit;
    %% STEP EX3: Plot Positions in 2D contour (Visualization Purposes)
%         contour(X,Y,Z);                                                             % Generate 2D Function Contour using previously generated data
%         hold on                                                                     % Hold ON graphics
%         plot(H(:,1),H(:,2),'r.',P(:,1),P(:,2),...                                   % Plot positions from "H" (red dots) and "P" (blue crosses) and Global Best Position (Green circle) for the current iteration
%             'bx',xbest(:,1),xbest(:,2),'go');
%         axis([xd xu xd xu])                                                         % Set Plot's Axis using the user specified function bounds [xd,xu] 
%         grid on                                                                     % Enable Plot's Grid
%         drawnow                                                                     % Draw Pending Graphics
%         hold off                                                                    % Hold OFF graphics
    end                                                                         %END
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=bfith;
end                                                                         % END of Function (SHO)

function [SV_h, SV_p] = SV(bfit, wfit, f_h, f_p)
% Calculate Survival Values for all members in "H" and "P"
    SV_h = (f_h - wfit)./(bfit - wfit);
    SV_p = (f_p - wfit)./(bfit - wfit);
% NaN and Inf values are set to 1 and 0 respectively (Sigular cases)
    SV_h(isnan(SV_h)) = 1;  SV_h(isinf(SV_h)) = 0;
    SV_p(isnan(SV_p)) = 1;  SV_p(isinf(SV_p)) = 0;
end

function [H] = HerdMove(H,SV_h,P,SV_p,xbest,dims,xd,xu)
    scale = abs(xu-xd);
    hL_idx = find(SV_h == max(SV_h));
    if size(hL_idx,1) > 1
        hL_idx = hL_idx(randi([1,size(hL_idx,1)]));
    end
    HF_idx = find(1:size(H,1)~=hL_idx);
    SV_mean = mean(SV_h);
    SV_h_wgt = repmat(SV_h,[1,dims]);
    h_M = sum(SV_h_wgt.*H)./(sum(SV_h_wgt)+.00001);
    SV_p_wgth = repmat(SV_p,[1,dims]);
    p_M = sum(SV_p_wgth.*P)./(sum(SV_p_wgth)+.00001);
%% Herd Leader Movement Operators
    h_L = H(hL_idx,:);
    if SV_h(hL_idx) == 1
    % Seemingly Cooperative Leadership Movement
        v_bi = p_M - h_L;
        r_bi = sqrt(sum(v_bi.^2,2))./scale;
        W_L = exp(-r_bi^2);
        move_L = -2*W_L*(v_bi);
    else
    % Openly Selfish Leadership Movement
        v_bi = xbest - h_L;
        r_bi = sqrt(sum(v_bi.^2,2))./scale;
        W_L = exp(-r_bi^2);
        move_L = 2*W_L*(v_bi);
    end
    alpha = rand(1,dims);
    H(hL_idx,:) = H(hL_idx,:) + alpha.*move_L;
    
%% Herd Members Movement Operators                                                  
    for i = HF_idx
        if SV_h(i)>=rand
        % Apply Herd Following Members Movement Operators
            if SV_h(i)>SV_mean
            % Apply Nearest Neighbor Movement Rule
                HF_idx = find(SV_h(HF_idx) > SV_h(i));
                if isempty(HF_idx)
                    move_iCi = zeros(1,dims);
                else
                    v_bi = H(HF_idx,:) - repmat(H(i,:),[size(HF_idx,1),1]);
                    r_bi = sqrt(sum(v_bi.^2,2))./scale;
                    Ci_idx = find(r_bi == min(r_bi));
                    if size(Ci_idx,1) > 1
                        Ci_idx = Ci_idx(randi([1,size(Ci_idx,1)]));
                    end
                    h_Ci = H(HF_idx(Ci_idx),:);
                    W_iCi = SV_h(HF_idx(Ci_idx))*exp(-(r_bi(Ci_idx).^2));
                    move_iCi = W_iCi*(h_Ci-H(i,:));
                end
                r_iL = norm(h_L - H(i,:))./scale;
                W_iL = SV_h(hL_idx)*exp(-(r_iL.^2));
                move_iL = W_iL*(h_L-H(i,:));
                beta = rand(1,dims);
                gamma = rand(1,dims);
                H(i,:) = H(i,:) + 2*(beta.*move_iL + gamma.*move_iCi);
            else
            %Apply Local Crowded Horizon Movement Rule
                delta = rand(1,dims);
                move_iM = (h_M - H(i,:));
                H(i,:) = H(i,:) + 2*delta.*move_iM;
            end
        else
        % Apply Herd Desertion Movement Operator
            V = randn(1,dims);
            v_i = (V./sqrt(V*V'))/scale;
            W_rand = (1 - SV_h(i));
            rand_move = W_rand*v_i;
            v_bi = xbest - H(i,:);
            r_bi = sqrt(sum(v_bi.^2,2))./scale;
            W_dp = exp(-r_bi^2);
            move_bi = W_dp*(v_bi);
            beta = rand(1,dims);
            gamma = rand(1,dims);
            H(i,:) = H(i,:) + 2*(beta.*move_bi + gamma.*rand_move);
        end
    end
    H(H<xd) = xd; H(H>xu) = xu;
end

function [P] = PredMove(P,H,SV_h,dims,xd,xu)
    for i = 1:size(P,1)
        PV_h = 1 - SV_h;
        h_r = H(roulette_Selection(PV_h'),:);
        move_ir = (h_r - P(i,:));
        omega = rand(1,dims);
        P(i,:) = P(i,:)+ 2*omega.*move_ir;
    end
    P(P<xd) = xd; P(P>xu) = xu;
end

function [L_idx,K_idx] = Predation(P,H,SV_p,SV_h,dims)
    R = sum(max(H) - min(H))/(2*dims);
    L_idx = (1:size(H,1))';
    K_idx = [];
    K_count = 0;
    for i = 1:size(P,1)
        Weak_idx = L_idx(SV_h(L_idx) < SV_p(i));
        if isempty(Weak_idx) ~= 1
            v_ij = H(Weak_idx,:) - repmat(P(i,:),[size(Weak_idx,1),1]);
            r_ij = sqrt(sum(v_ij.^2,2));
            Tpi_idx = find(r_ij < R);
            if isempty(Tpi_idx) ~= 1
                H_pihi = (1 - SV_h(Weak_idx(Tpi_idx))).*exp(-(r_ij(Tpi_idx).^2));
                target = roulette_Selection(H_pihi');            
                K_count = K_count + 1;
                K_idx(K_count,1) = Weak_idx(target);
                L_idx = L_idx(ismember(L_idx,K_idx(K_count,1))==0);
            else
                % PREDATOR "i" doesn't "KILL" any PREY
            end
        else
            % PREDATOR "i" doesn't "KILL" any PREY
        end
    end
end                                                                        

function [H, f_h] = Restoration(H,f_h,SV_h,L_idx,K_idx,dims,f)
    alive_H = H(L_idx,:);
    alive_SVh = SV_h(L_idx); 
    for i = 1:size(K_idx,1)
        for gen = 1:dims
            pose_i = roulette_Selection(alive_SVh');
            new_h_aux(1,gen) = alive_H(pose_i,gen);
        end
        new_h(i,:) = new_h_aux(randperm(length(new_h_aux)));
        new_fh(i) = f(new_h(i,:));           
    end
    H(K_idx,:) = new_h;
    f_h(K_idx) = new_fh;
end

function selection = roulette_Selection(data)
    roulette = cumsum(data/sum(data));
    trial = rand;
    selection = size(find(roulette<trial),2) + 1;
end

function [X,Y,Z] = contour_plot(f, xd, xu)
    n_samples = 50;
    x = linspace(xd,xu,n_samples);
    y = linspace(xd,xu,n_samples);
    [X,Y] = meshgrid(x,y);
    Z = zeros(n_samples);
    for i = 1:n_samples
        for j = 1:n_samples
            Z(i,j) = f([x(i),y(j)]);
        end
    end
end