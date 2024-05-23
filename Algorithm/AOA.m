function Data_o = AOA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
Materials_no=Data_i.pop;
C1=2;C2=6;
C3=1;C4=2;  %standard Optimization functions

u=.9;l=.1;   %paramters in Eq. (12)
X=Data_i.X;%initial positions Eq. (4)
den=rand(Materials_no,Data_i.dim);           % Eq. (5)
vol=rand(Materials_no,Data_i.dim);
acc=Data_i.lb+rand(Materials_no,Data_i.dim).*(Data_i.ub-Data_i.lb);% Eq. (6)
Y=Data_i.F_value;

[Scorebest, Score_index] = min(Y);
Xbest = X(Score_index,:);
den_best=den(Score_index,:);
vol_best=vol(Score_index,:);
acc_best=acc(Score_index,:);
acc_norm=acc;
Max_iter=Data_i.maxIter;

for t = 1:Max_iter
    TF=exp(((t-Max_iter)/(Max_iter)));   % Eq. (8)
    if TF>1
        TF=1;
    end
    d=exp((Max_iter-t)/Max_iter)-(t/Max_iter); % Eq. (9)
    acc=acc_norm;
    r=rand();
    for i=1:Data_i.pop
        den(i,:)=den(i,:)+r*(den_best-den(i,:));   % Eq. (7)
        vol(i,:)=vol(i,:)+r*(vol_best-vol(i,:));
        if TF<.45%collision
            mr=randi(Data_i.pop);
            acc_temp(i,:)=(den(mr,:)+(vol(mr,:).*acc(mr,:)))./(rand*den(i,:).*vol(i,:));   % Eq. (10)
        else
            acc_temp(i,:)=(den_best+(vol_best.*acc_best))./(rand*den(i,:).*vol(i,:));   % Eq. (11)
        end
    end
    
    acc_norm=((u*(acc_temp-min(acc_temp(:))))./(max(acc_temp(:))-min(acc_temp(:))))+l;   % Eq. (12)

    for i=1:Data_i.pop
        if TF<.4
            for j=1:size(X,2)
                mrand=randi(Materials_no);
                Xnew(i,j)=X(i,j)+C1*rand*acc_norm(i,j).*(X(mrand,j)-X(i,j))*d;  % Eq. (13)
            end
        else
            for j=1:size(X,2)
                p=2*rand-C4;  % Eq. (15)
                T=C3*TF;
                if T>1
                    T=1;
                end
                if p<.5
                    Xnew(i,j)=Xbest(j)+C2*rand*acc_norm(i,j).*(T*Xbest(j)-X(i,j))*d;  % Eq. (14)
                else
                    Xnew(i,j)=Xbest(j)-C2*rand*acc_norm(i,j).*(T*Xbest(j)-X(i,j))*d;
                end
            end
        end 
        Xnew(i,:)=BoundaryCheck(Xnew(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
    end

    for i=1:Data_i.pop
        v=Data_i.fobj(Xnew(i,:));
        if v<Y(i)
            X(i,:)=Xnew(i,:);
            Y(i)=v;
        end
    end
    [var_Ybest,var_index] = min(Y);
    IterCurve(t)=var_Ybest;
    if var_Ybest<Scorebest              
        Scorebest=var_Ybest;
        Score_index=var_index;
        Xbest = X(var_index,:);
        den_best=den(Score_index,:);
        vol_best=vol(Score_index,:);
        acc_best=acc_norm(Score_index,:);

        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=var_Ybest;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=var_index;
    end
end    
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end