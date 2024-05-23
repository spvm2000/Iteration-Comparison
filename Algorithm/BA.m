function Data_o = BA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 
Sol=Data_i.X;        Fitness=Data_i.F_value;
[fmin,I]=min(Fitness);
best=Sol(I,:);
A=1;            % Initial loudness (constant or decreasing)
r0=1;           % The initial pulse rate (constant or decreasing)
alpha=0.97;     % Parameter alpha
gamma=0.1;      % Parameter gamma
% Frequency range
Freq_min=0;     % Frequency minimum
Freq_max=2;     % Frequency maximum
% 初始化全局最优适应度值
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fmin;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=I;

d=Data_i.dim;       n=Data_i.pop;
% Initialization of all the arrays
Freq=zeros(n,1);   % Frequency-tuning array
v=zeros(n,d);      % Equivalnet velocities or increments
Lb=Data_i.lb;      % Lower bounds
Ub=Data_i.ub;      % Upper bounds
IterCurve=zeros(1,Data_i.maxIter);

% Start the iterations -- Bat Algorithm (essential part)  %
for t=1:Data_i.maxIter 
    r=r0*(1-exp(-gamma*t));
    A=alpha*A;
    for i=1:n
        Freq(i)=Freq_min+(Freq_max-Freq_min)*rand;
        v(i,:)=v(i,:)+(Sol(i,:)-best)*Freq(i);
        S(i,:)=Sol(i,:)+v(i,:);
        % Check a switching condition
        if rand<r
            S(i,:)=best+0.1*randn(1,d)*A;
        end
        
        % Check if the new solution is within the simple bounds
        S(i,:)=BoundaryCheck(S(i,:),Ub,Lb,d);
        % Evaluate new solutions
        Fnew=Data_i.fobj(S(i,:));
        % If the solution improves or not too loudness
        if ((Fnew<=Fitness(i)) && (rand>A))
            Sol(i,:)=S(i,:);
            Fitness(i)=Fnew;
        end
        % Update the current best solution
        if Fnew<=fmin
            best=S(i,:);
            fmin=Fnew;
        end
    end
   Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fmin;
   IterCurve(t)=fmin;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound vector
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bound vector 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
end
%%%%% ============ end ====================================


