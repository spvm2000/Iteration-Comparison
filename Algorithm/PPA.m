function Data_o = PPA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
nest = Data_i.X;                       
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
n =Data_i.pop;      lb=Data_i.lb;   ub=Data_i.ub; d=Data_i.dim;     fobj=Data_i.fobj;
maxiter=Data_i.maxIter; 

[fitness,I]=sort(fitness);nest=nest(I,:);
fmin=fitness(1);  bestnest=nest(1,:);
wMax = 0.9;
wMin = 0.3;
w = linspace(wMax, wMin, maxiter);
Cats_v = 0.25 * nest;
[GrowthRateCrows,GrowthRateCats,GrowthRateCuckoos] =Growth_rate(n,maxiter);
for t=1:maxiter
    nCats=GrowthRateCats(t);%current no. of Cats
    nCrows=GrowthRateCrows(t);%current no. of Crows
    nCuckoos=GrowthRateCuckoos(t);

    new_Crows=get_Crows(nest(1:nCrows,:),lb,ub);

    Rolette_index = RankingSelection(n, nCuckoos);
    parasitized_nests = nest(Rolette_index(1:nCuckoos),:);
    new_Cuckoos = get_Cuckoos(parasitized_nests,lb,ub,t,maxiter);

    non_parasitized_nests = randperm( n);
    non_parasitized_nests(Rolette_index(1:nCuckoos)) = [];
    Cats_nests = nest(non_parasitized_nests( 1:nCats),:); % Choosing
    new_Cats_v = Cats_v(non_parasitized_nests( 1:nCats),:);
    [new_Cats, new_Cats_v] = get_Cats(new_Cats_v, Cats_nests, bestnest,t,maxiter, w(t),lb,ub);
    
    Cats_v(non_parasitized_nests( 1:nCats),:) = new_Cats_v;
    newnest = [new_Crows;new_Cuckoos;new_Cats];

    [nest,fitness] = get_best_nest(fobj,nest,newnest,fitness);
    [fitness,I] = sort(fitness);
    nest = nest(I,:);
    fmin = fitness(1);
    bestnest = nest(1,:);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fmin
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fmin;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=I(1);
    end
    IterCurve(t) = fmin;
end    
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [Cats, Cats_v] = ...
    get_Cats(Cats_v, Cats, gx,t,maxiter, w,lb,ub)
perCnt = 1-t/maxiter/4; 
vlb = perCnt*lb;
vUp = perCnt*ub;
c=2-t/maxiter; % c (social  learning) is modified to be decreased linearly from 2 to 1

for iSz2 = 1:size(Cats,1)
    V_vec = w  * Cats_v(iSz2,:) + c*rand*(gx - Cats(iSz2,:)); % Eq. (13)
    % Apply simple bounds/limits to velocity of cats
    V_vec = max(V_vec, vlb);
    V_vec = min(V_vec, vUp);
    Cats_v(iSz2,:) = V_vec;
end
Cats = Cats + Cats_v; % Eq. (14)
Cats = simplebounds( Cats, lb, ub);
end

function [nest,fitness]=get_best_nest(fobj,nest,newnest,fitness)
for j=1:size(nest,1)
     fnew = fobj(newnest(j,:)); 
    if fnew<=fitness(j)
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end
end

function [GrowthRateCrows,GrowthRateCats, GrowthRateCuckoos] =Growth_rate(n,maxiter)
    GrowthRateCrows= round(n*linspace(2/3,1/2, maxiter));%Growth rate of Crows
    GrowthRateCats= round(n*linspace(0.01,1/3, maxiter));%Growth rate of Cats
    GrowthRateCuckoos=n-GrowthRateCrows-GrowthRateCats;
end


function nest=get_Crows(nest,lb,ub)
n=size(nest,1);
beta=3/2; sigma=0.6965745;
num=ceil(n*rand(1,n));
    for i=1:n
        s=nest(i,:);
        u=randn(size(s))*sigma;
        v=randn(size(s));
        step=u./abs(v).^(1/beta);           % Eq.(7)
        stepsize=0.1*step.*randn(size(s));  % Eq.(8)
        s= s+stepsize.*(nest(num(i),:)-s); % Eq.(6)
        out=s<lb | s > ub;
        s(:,out)=lb(out)+(ub(out)-lb(out)).*rand(1,nnz(out)); % Eq.(9)
        nest(i,:)=s;
    end
end

function choice = RankingSelection(n, N)
    ranking = (1:n);% for sorted fitness
    weights=ranking/sum(ranking);
    Select_Fitness = cumsum(weights);
    choice=[]; 
    while length(choice) < N
        Random_Fitness=rand(1,n);
        selected = ranking(Select_Fitness <= Random_Fitness);
        choice=union(choice,selected,'stable')';
    end
end

function new_nest=get_Cuckoos(nest,lb,ub,t,maxiter)
    pa=t/2/maxiter; 
    n=size(nest,1);
    % probability  Eq. (12)
    K=rand(size(nest)) > pa; % Discovered or not -- a status vector
    stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
    new_nest=nest+stepsize.*K;  % Eq. (10)
    % Apply simple bounds/limits
    new_nest=simplebounds(new_nest,lb,ub);
end

function Sol=simplebounds(Sol,lb,ub)
for i=1:size(Sol,1)% return back populations
   Vub=Sol(i,:)>ub;
   Vlb=Sol(i,:)<lb;
   Sol(i,:)=(Sol(i,:).*(~(Vub+Vlb)))+(lb+rand(size(lb)).*(ub-lb)).*(Vub+Vlb);
end  
end