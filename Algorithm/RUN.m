function Data_o = RUN(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
Cost=fitness;
Xnew2=zeros(1,Data_i.dim);
Best_X = X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);
t=0;
while t<Data_i.maxIter
    f=20.*exp(-(12.*(t/Data_i.maxIter)));  % (Eq.17.6) 
    Xavg = mean(X);                
    SF=2.*(0.5-rand(1,Data_i.pop)).*f;     %  (Eq.17.5)

    for i=1:Data_i.pop
        [~,ind_l] = min(Cost);
        lBest = X(ind_l,:);

        [A,B,C]=RndX(Data_i.pop,i);         % Determine Three Random Indices of Solutions
        [~,ind1] = min(Cost([A B C]));

        gama = rand.*(X(i,:)-rand(1,Data_i.dim).*(Data_i.ub-Data_i.lb)).*exp(-4*t/Data_i.maxIter);  
        Stp=rand(1,Data_i.dim).*((Best_X-rand.*Xavg)+gama);
        DelX = 2*rand(1,Data_i.dim).*(abs(Stp));

        if Cost(i)<Cost(ind1)                
            Xb = X(i,:);
            Xw = X(ind1,:);
        else
            Xb = X(ind1,:);
            Xw = X(i,:);
        end

        SM = RungeKutta(Xb,Xw,DelX);
        L=rand(1,Data_i.dim)<0.5;
        Xc = L.*X(i,:)+(1-L).*X(A,:);  % (Eq. 17.3)
        Xm = L.*Best_X+(1-L).*lBest;   % (Eq. 17.4)
          
        vec=[1,-1];
        flag = floor(2*rand(1,Data_i.dim)+1);
        r=vec(flag);                   % An Interger number 
        
        g = 2*rand;
        mu = 0.5+.1*randn(1,Data_i.dim);

        if rand<0.5
            Xnew = (Xc+r.*SF(i).*g.*Xc) + SF(i).*(SM) + mu.*(Xm-Xc);
        else
            Xnew = (Xm+r.*SF(i).*g.*Xm) + SF(i).*(SM)+ mu.*(X(A,:)-X(B,:));
        end 

        FU=Xnew>Data_i.ub;FL=Xnew<Data_i.lb;Xnew=(Xnew.*(~(FU+FL)))+Data_i.ub.*FU+Data_i.lb.*FL; 
        CostNew=Data_i.fobj(Xnew);
        
        
        % update locate value
        if CostNew<Cost(i)
            X(i,:)=Xnew;
            Cost(i)=CostNew;
        end
        t=t+1;
        if t<=Data_i.maxIter
            IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
        else
            break;
        end

        if rand<0.5
            EXP=exp(-5*rand*t/Data_i.maxIter);
            r = floor(Unifrnd(-1,2,1,1));

            u=2*rand(1,Data_i.dim); 
            w=Unifrnd(0,2,1,Data_i.dim).*EXP;               %(Eq.19-1)
            
            [A,B,C]=RndX(Data_i.pop,i);
            Xavg=(X(A,:)+X(B,:)+X(C,:))/3;           %(Eq.19-2)         
            
            beta=rand(1,Data_i.dim);
            Xnew1 = beta.*(Best_X)+(1-beta).*(Xavg); %(Eq.19-3)

            for j=1:Data_i.dim
                if w(j)<1 
                    Xnew2(j) = Xnew1(j)+r*w(j)*abs((Xnew1(j)-Xavg(j))+randn);
                else
                    Xnew2(j) = (Xnew1(j)-Xavg(j))+r*w(j)*abs((u(j).*Xnew1(j)-Xavg(j))+randn);
                end
            end
            
            FU=Xnew2>Data_i.ub;FL=Xnew2<Data_i.lb;Xnew2=(Xnew2.*(~(FU+FL)))+Data_i.ub.*FU+Data_i.lb.*FL;
            CostNew=Data_i.fobj(Xnew2);
            t=t+1;
            if t<=Data_i.maxIter
                IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
            else
                break;
            end
            % update global and locate value
            if CostNew<Cost(i)
               X(i,:)=Xnew2;
               Cost(i)=CostNew;
               if Cost(i)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
                   Best_X=X(i,:);
                   Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Cost(i);
                   Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
               end
               IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
            else
                if rand<w(randi(Data_i.dim)) 
                    SM = RungeKutta(X(i,:),Xnew2,DelX);
                    Xnew = (Xnew2-rand.*Xnew2)+ SF(i)*(SM+(2*rand(1,Data_i.dim).*Best_X-Xnew2));  % (Eq. 20)
                    
                    FU=Xnew>Data_i.ub;FL=Xnew<Data_i.lb;Xnew=(Xnew.*(~(FU+FL)))+Data_i.ub.*FU+Data_i.lb.*FL;
                    CostNew=Data_i.fobj(Xnew);
                    
                    % update locate value
                    if CostNew<Cost(i)
                        X(i,:)=Xnew;
                        Cost(i)=CostNew;
                    end
                    if Cost(i)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
                       Best_X=X(i,:);
                       Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Cost(i);
                       Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
                   end
                   IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);

                    t=t+1;
                    if t<=Data_i.maxIter
                        IterCurve(t)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
                    else
                        break;
                    end
                end
            end
        end    
    end
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [A,B,C]=RndX(nP,i)
    Qi=randperm(nP);Qi(Qi==i)=[];
    A=Qi(1);B=Qi(2);C=Qi(3);
end

function z=Unifrnd(a,b,c,dim)
    a2 = a/2;
    b2 = b/2;
    mu = a2+b2;
    sig = b2-a2;
    z = mu + sig .* (2*rand(c,dim)-1);
end

function SM=RungeKutta(XB,XW,DelX)
    dim=size(XB,2);
    C=randi([1 2])*(1-rand);
    r1=rand(1,dim);
    r2=rand(1,dim);
    
    K1 = 0.5*(rand*XW-C.*XB);
    K2 = 0.5*(rand*(XW+r2.*K1.*DelX/2)-(C*XB+r1.*K1.*DelX/2));
    K3 = 0.5*(rand*(XW+r2.*K2.*DelX/2)-(C*XB+r1.*K2.*DelX/2));
    K4 = 0.5*(rand*(XW+r2.*K3.*DelX)-(C*XB+r1.*K3.*DelX));
    
    XRK = (K1+2.*K2+2.*K3+K4);
    SM=1/6*XRK;
end
