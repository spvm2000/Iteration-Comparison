function Data_o = RSA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
Ffun=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Ffun);
IterCurve=zeros(1,Data_i.maxIter);
t=0;                         % starting iteration
Alpha=0.1;                   % the best value 0.1
Beta=0.1;                    % the best value 0.005
Xnew=X;
Ffun_new=Ffun;
Best_F=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
Best_P=X(Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter),:);

while t<Data_i.maxIter
    ES=2*randi([-1 1])*(1-(t/Data_i.maxIter));  % Probability Ratio
    for i=2:size(X,1) 
        for j=1:size(X,2)  
            R=Best_P(1,j)-X(randi([1 size(X,1)]),j)/((Best_P(1,j))+eps);
            P=Alpha+(X(i,j)-mean(X(i,:)))/(Best_P(1,j)*(Data_i.ub(j)-Data_i.lb(j))+eps);
            Eta=Best_P(1,j)*P;
            if (t<Data_i.maxIter/4)
                Xnew(i,j)=Best_P(1,j)-Eta*Beta-R*rand;    
            elseif (t<2*Data_i.maxIter/4 && t>=Data_i.maxIter/4)
                Xnew(i,j)=Best_P(1,j)*X(randi([1 size(X,1)]),j)*ES*rand;
            elseif (t<3*Data_i.maxIter/4 && t>=2*Data_i.maxIter/4)
                Xnew(i,j)=Best_P(1,j)*P*rand;
            else
                Xnew(i,j)=Best_P(1,j)-Eta*eps-R*rand;
            end
            Flag_UB=Xnew(i,:)>Data_i.ub; % check if they exceed (up) the boundaries
            Flag_LB=Xnew(i,:)<Data_i.lb; % check if they exceed (down) the boundaries
            Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+Data_i.ub.*Flag_UB+Data_i.lb.*Flag_LB;
            Ffun_new(1,i)=Data_i.fobj(Xnew(i,:));
            if Ffun_new(1,i)<Ffun(1,i)
                X(i,:)=Xnew(i,:);
                Ffun(1,i)=Ffun_new(1,i);
            end
            if Ffun(1,i)<Best_F                 %global best value
                Best_F=Ffun(1,i);
                Best_P=X(i,:);
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_F;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
            t=t+1;
            if t>Data_i.maxIter
               break;
            else
               IterCurve(t)=Best_F;  %Update the convergence curve
            end
            
        end
    end    
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end