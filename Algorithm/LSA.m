function Data_o = LSA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
Dpoint = Data_i.X;                        
fitness=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
LB=Data_i.lb;      UB=Data_i.ub;       dim=Data_i.dim;     D=dim;  T=Data_i.maxIter;
ch_time = 0; % reset        
max_ch_time = 10;   N=Data_i.pop;
direct = sign(unifrnd(-1,1,1,dim));
for t = 1:T
    for i=1:Data_i.pop
        fit(i) = Data_i.fobj(Dpoint(i,:));
    end
    Ec = fit;

    ch_time = ch_time+1;
    if ch_time >= max_ch_time
        [Ms ds]=sort(Ec,'ascend');
        Dpoint(ds(N),:) = Dpoint(ds(1),:); % Eliminate the worst channel
        Ec(ds(N)) = Ec(ds(1)); % Update  
        ch_time = 0; % reset
    end

    [Ms ds]=sort(Ec,'ascend');
    best = Ec(ds(1));
    worst = Ec(ds(N));

    Energy = 2.05 - 2*exp(-5*(T-t)/T);

    for d = 1:D
        Dpoint_test = Dpoint(ds(1),:);
        Dpoint_test(d) = Dpoint_test(d)+direct(d)*0.005*(UB(d)-LB(d));
        fv_test = Data_i.fobj(Dpoint_test);
        if fv_test < best % If better, +ve direction
            direct(d) = direct(d);
        else
            direct(d) = -1*direct(d);
        end
    end

    for i = 1:N
        dist=Dpoint(i,:)-Dpoint(ds(1),:);
        for d = 1:D
            if Dpoint(i,:)==Dpoint(ds(1),:)
                Dpoint_temp(d) = Dpoint(i,d)+direct(d)*abs(normrnd(0,Energy));
            else
                if dist(d)<0
                    Dpoint_temp(d) = Dpoint(i,d)+exprnd(abs(dist(d)));
                else
                    Dpoint_temp(d) = Dpoint(i,d)-exprnd(dist(d));
                end
            end
            if (Dpoint_temp(d)>UB(d))||(Dpoint_temp(d)<LB(d))
                Dpoint_temp(d) = rand(1)*(UB(d)-LB(d))+LB(d); % Re-initialized
            end
        end
        
        fv = Data_i.fobj(Dpoint_temp);
        if fv < Ec(i)
            Dpoint(i,:) = Dpoint_temp;
            Ec(i) = fv;
            % Focking procedure
            if rand < 0.01
                for d = 1:D 
                    Dpoint_fock(d) = UB(d)+LB(d)-Dpoint_temp(d);% Focking
                end
                fock_fit = Data_i.fobj(Dpoint_fock); % Evaluate
                if fock_fit < Ec(i) 
                    Dpoint(i,:) = Dpoint_fock; % Replace the channel
                    Ec(i) = fock_fit;
                end
            end
        end   
    end
    IterCurve(t) = min(fit);
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end