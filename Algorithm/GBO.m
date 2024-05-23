function Data_o = GBO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Cost);
IterCurve=zeros(1,Data_i.maxIter);
nV = Data_i.dim;                        % Number f Variables
pr = 0.5;                               % Probability Parameter
[~,Ind] = sort(Cost);     
Best_Cost = Cost(Ind(1));        % Determine the vale of Best Fitness
Best_X = X(Ind(1),:);
Worst_Cost=Cost(Ind(end));       % Determine the vale of Worst Fitness
Worst_X=X(Ind(end),:);
it=0;
while it<Data_i.maxIter
    beta = 0.2+(1.2-0.2)*(1-(it/Data_i.maxIter)^3)^2;                        % Eq.(14.2)
    alpha = abs(beta.*sin((3*pi/2+sin(3*pi/2*beta))));              % Eq.(14.1)
    
    for i=1:Data_i.pop
        A1 = fix(rand(1,Data_i.pop)*Data_i.pop)+1;                                  % Four positions randomly selected from population        
        r1 = A1(1);r2 = A1(2);   
        r3 = A1(3);r4 = A1(4); 

        Xm = (X(r1,:)+X(r2,:)+X(r3,:)+X(r4,:))/4;                   % Average of Four positions randomly selected from population        
        ro = alpha.*(2*rand-1);ro1 = alpha.*(2*rand-1);        
        eps = 5e-3*rand;

        DM = rand.*ro.*(Best_X-X(r1,:));Flag = 1;                   % Direction of Movement Eq.(18)
        GSR=GradientSearchRule(ro1,Best_X,Worst_X,X(i,:),X(r1,:),DM,eps,Xm,Flag);      
        DM = rand.*ro.*(Best_X-X(r1,:));
        X1 = X(i,:) - GSR + DM;

        DM = rand.*ro.*(X(r1,:)-X(r2,:));Flag = 2;
        GSR=GradientSearchRule(ro1,Best_X,Worst_X,X(i,:),X(r1,:),DM,eps,Xm,Flag); 
        DM = rand.*ro.*(X(r1,:)-X(r2,:));
        X2 = Best_X - GSR + DM;

        Xnew=zeros(1,Data_i.dim);
        for j=1:Data_i.dim                                                  
            ro=alpha.*(2*rand-1);                       
            X3=X(i,j)-ro.*(X2(j)-X1(j));           
            ra=rand;rb=rand;
            Xnew(j) = ra.*(rb.*X1(j)+(1-rb).*X2(j))+(1-ra).*X3;     % Eq.(27)          
        end

        if rand<pr           
            k = fix(rand*Data_i.pop)+1;
            f1 = -1+(1-(-1)).*rand();f2 = -1+(1-(-1)).*rand();         
            ro = alpha.*(2*rand-1);
            Xk = unifrnd(Data_i.lb,Data_i.ub,1,Data_i.dim);    % Eq.(28.8)

            L1=rand<0.5;u1 = L1.*2*rand+(1-L1).*1;u2 = L1.*rand+(1-L1).*1;
            u3 = L1.*rand+(1-L1).*1;                                    
            L2=rand<0.5;            
            Xp = (1-L2).*X(k,:)+(L2).*Xk;                           % Eq.(28.7)
                                                 
            if u1<0.5
                Xnew = Xnew + f1.*(u1.*Best_X-u2.*Xp)+f2.*ro.*(u3.*(X2-X1)+u2.*(X(r1,:)-X(r2,:)))/2;     
            else
                Xnew = Best_X + f1.*(u1.*Best_X-u2.*Xp)+f2.*ro.*(u3.*(X2-X1)+u2.*(X(r1,:)-X(r2,:)))/2;   
            end

            % Check if solutions go outside the search space and bring them back
            Flag4ub=Xnew>Data_i.ub;
            Flag4lb=Xnew<Data_i.lb;
            Xnew=(Xnew.*(~(Flag4ub+Flag4lb)))+Data_i.ub.*Flag4ub+Data_i.lb.*Flag4lb;  
            Xnew_Cost=Data_i.fobj(Xnew);
            
            % Update the Best Position        
            if Xnew_Cost<Cost(i)
                X(i,:)=Xnew;
                Cost(i)=Xnew_Cost;
                if Cost(i)<Best_Cost
                    Best_X=X(i,:);
                    Best_Cost=Cost(i);
                    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Cost;
                    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                end            
            end
           % Update the Worst Position 
            if Cost(i)>Worst_Cost
                Worst_X= X(i,:);
                Worst_Cost= Cost(i);
            end
        end
    end
    it=it+1;
    if it>Data_i.maxIter
         break;
    else
        IterCurve(it)=Best_Cost;
    end
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end


function GSR=GradientSearchRule(ro1,Best_X,Worst_X,X,Xr1,DM,eps,Xm,Flag)
    nV = size(X,2);
    Delta = 2.*rand.*abs(Xm-X);                            % Eq.(16.2)
    Step = ((Best_X-Xr1)+Delta)/2;                         % Eq.(16.1)
    DelX = rand(1,nV).*(abs(Step));                        % Eq.(16)
    
    GSR = randn.*ro1.*(2*DelX.*X)./(Best_X-Worst_X+eps);   % Gradient search rule  Eq.(15)
    if Flag == 1
      Xs = X - GSR + DM;                                   % Eq.(21)
    else
      Xs = Best_X - GSR + DM;
    end    
    yp = rand.*(0.5*(Xs+X)+rand.*DelX);                    % Eq.(22.6)
    yq = rand.*(0.5*(Xs+X)-rand.*DelX);                    % Eq.(22.7)
    GSR = randn.*ro1.*(2*DelX.*X)./(yp-yq+eps);            % Eq.(23)   
end