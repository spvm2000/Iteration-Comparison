function Data_o = BWO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 

%% Problem Definition
Max_it = Data_i.maxIter;      
Npop = Data_i.pop;            
nD = Data_i.dim;              
pos = Data_i.X;               
fit = Data_i.F_value;         

Curve = inf*ones(1,Max_it);
[fvalbest,index]=min(fit);
xposbest = pos(index,:);
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);

T = 1;
while T <= Max_it
    newpos = pos;
    WF = 0.1-0.05*(T/Max_it);  % The probability of whale fall
    kk = (1-0.5*T/Max_it)*rand(Npop,1); % The probability in exploration or exploitation
    rand_num(T,:) = kk;

    for i = 1:Npop
        if kk(i) > 0.5 % exploration phase
            r1 = rand(); r2 = rand();
            RJ = ceil(Npop*rand);   % Roulette Wheel Selection
            while RJ == i
                RJ = ceil(Npop*rand);
            end
            if nD <= Npop/5
                params = randperm(nD,2);
                newpos(i,params(1)) = pos(i,params(1))+(pos(RJ,params(1))-pos(i,params(2)))*(r1+1)*sin(r2*360);
                newpos(i,params(2)) = pos(i,params(2))+(pos(RJ,params(1))-pos(i,params(2)))*(r1+1)*cos(r2*360);
            else
                params=randperm(nD);
                for j = 1:floor(nD/2)
                    newpos(i,2*j-1) = pos(i,params(2*j-1))+(pos(RJ,params(1))-pos(i,params(2*j-1)))*(r1+1)*sin(r2*360);
                    newpos(i,2*j) = pos(i,params(2*j))+(pos(RJ,params(1))-pos(i,params(2*j)))*(r1+1)*cos(r2*360);
                end
            end
        else  % exploitation phase
            r3 = rand(); r4 = rand(); C1 = 2*r4*(1-T/Max_it);
            RJ = ceil(Npop*rand);   % Roulette Wheel Selection
            while RJ == i
                RJ = ceil(Npop*rand);
            end
            alpha=3/2;
            sigma=(gamma(1+alpha)*sin(pi*alpha/2)/(gamma((1+alpha)/2)*alpha*2^((alpha-1)/2)))^(1/alpha); % Levy flight
            u=randn(1,nD).*sigma;
            v=randn(1,nD);
            S=u./abs(v).^(1/alpha);
            KD = 0.05;
            LevyFlight=KD.*S;
            newpos(i,:) = r3*xposbest - r4*pos(i,:) + C1*LevyFlight.*(pos(RJ,:)-pos(i,:));
        end
        % boundarycheck
        newpos = BoundaryCheck(newpos,Data_i.ub,Data_i.lb,Data_i.dim);
        newfit = Data_i.fobj(newpos(i,:));    % fitness calculation
        if newfit < fit(i)
            pos(i,:) = newpos(i,:);
            fit(i) = newfit;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    
    for i = 1:Npop
        % whale falls
        if kk(i) <= WF
            RJ = ceil(Npop*rand); r5 = rand(); r6 = rand(); r7 = rand();
            C2 = 2*Npop*WF;
            stepsize2 = r7*(Data_i.ub-Data_i.lb)*exp(-C2*T/Max_it);
            newpos(i,:) = (r5*pos(i,:) - r6*pos(RJ,:)) + stepsize2;
            % boundarycheck
            newpos = BoundaryCheck(newpos,Data_i.ub,Data_i.lb,Data_i.dim);
            newfit = Data_i.fobj(newpos(i,:));    % fitness calculation
            if newfit < fit(i)
                pos(i,:) = newpos(i,:);
                fit(i) = newfit;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
    end

    [fval,index]=min(fit);

    if fval<fvalbest
        fvalbest = fval;
        xposbest = pos(index,:);
    end
    
    T = T+1;
    if T<=Data_i.maxIter
        Curve(T) = fvalbest;
    end
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fvalbest;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=T;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Curve;
end
