function Data_o = SSO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
X = Data_i.X;                        
Ffun=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
befit=zeros(1,Data_i.maxIter);
spidn=Data_i.pop;        itern=Data_i.maxIter;      f=Data_i.fobj;  lb=Data_i.lb;
ub=Data_i.ub;       dims=Data_i.dim;
rand('state',0');
fpl = 0.65;     % Lower Female Percent
fpu = 0.9;      % Upper Female Percent
fp = fpl+(fpu-fpl)*rand;	% Aleatory Percent
fn = round(spidn*fp);   % Number of females
mn = spidn-fn;
pm = exp(-(0.1:(3-0.1)/(itern-1):3));

fsp = zeros(fn,dims);   % Initlize females
msp = zeros(mn,dims);   % Initlize males

fefit = zeros(fn,1);    % Initlize fitness females
mafit = zeros(mn,1);    % Initlize fitness males

spwei = zeros(spidn,1); % Initlize weigth spiders
fewei = zeros(fn,1);    % Initlize weigth spiders
mawei = zeros(mn,1);    % Initlize weigth spiders

for i=1:fn
    fsp(i,:)=X(i,:);
end

for i=fn+1:Data_i.pop
    msp(i-fn,:)=X(i,:);
end

for i=1:fn
    fefit(i)=f(fsp(i,:));
end

for i=1:mn
    mafit(i)=f(msp(i,:));
end

spfit = [fefit' mafit']';    % Mix Females and Males
bfitw = min(spfit);          % best fitness
wfit = max(spfit);           % worst fitness
for i=1:spidn
    spwei(i) = 0.001+((spfit(i)-wfit)/(bfitw-wfit));
end
fewei = spwei(1:fn);      % Separate the female mass
mawei = spwei(fn+1:spidn);

[~,Ibe] = max(spwei);
% Check if female or male
if Ibe > fn
    % Is Male
    spbest=msp(Ibe-fn,:);   % Asign best position to spbest
    bfit = mafit(Ibe-fn);      % Get best fitness for memory
else
    % Is Female
    spbest=fsp(Ibe,:);      % Asign best position to spbest
    bfit = fefit(Ibe);      % Get best fitness for memory
end

for i=1:itern
    [fsp] = FeMove(spidn,fn,fsp,msp,spbest,Ibe,spwei,dims,lb,ub,pm(i));
    [msp] = MaMove(fn,mn,fsp,msp,fewei,mawei,dims,lb,ub,pm(i));

    for j=1:fn                      
        fefit(j)=f(fsp(j,:));
    end

    for j=1:mn
        mafit(j)=f(msp(j,:));
    end

    spfit = [fefit' mafit']';    % Mix Females and Males
    bfitw = min(spfit);          % best fitness
    wfit = max(spfit);           % worst fitness
    % Obtain weight for every spider
    for j=1:spidn
        spwei(j) = 0.001+((spfit(j)-wfit)/(bfitw-wfit));
    end
    
    fewei = spwei(1:fn);      % Separate the female mass
    mawei = spwei(fn+1:spidn);% Separate the male mass
    [ofspr] = Mating(fewei,mawei,fsp,msp,dims);
    if isempty(ofspr)
    else
        [fsp,msp,fefit,mafit] = Survive(fsp,msp,ofspr,fefit,mafit,spfit,f,fn,dims);
        % ***** Recalculate the weigth or sort ***********
        spfit = [fefit' mafit']';    % Mix Females and Males
        bfitw = min(spfit);          % best fitness
        wfit = max(spfit);           % worst fitness
        % Obtain weight for every spider
        for j=1:spidn
            spwei(j) = 0.001+((spfit(j)-wfit)/(bfitw-wfit));
        end
        fewei = spwei(1:fn);      % Separate the female mass
        mawei = spwei(fn+1:spidn);% Separate the male mass
    end

    % Check if best position belongs to male or female
    [~,Ibe2] = max(spwei);
    if Ibe2 > fn
        % Is Male
        spbest2=msp(Ibe2-fn,:);      % Asign best position to spbest
        bfit2 = mafit(Ibe2-fn);      % Get best fitness for memory
    else
        % Is Female
        spbest2 = fsp(Ibe2,:);  % Asign best position to spbest
        bfit2 = fefit(Ibe2);    % Get best fitness for memory
    end
    %% Global Memory
    if bfit<=bfit2
        bfit = bfit;
        spbest = spbest;      % Asign best position to spbest
        befit(i) = bfit;
    else
        bfit = bfit2;
        spbest = spbest2;      % Asign best position to spbest
        befit(i) = bfit;
    end
end
IterCurve=befit;
if min(befit)<Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)
    [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),~]=min(befit);
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=i;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [fsp] = FeMove(spidn,fn,fsp,msp,spbest,Ibe,spmass,d,lb,ub,pm)
%FEMOVE Summary of this function goes here
%   Detailed explanation goes here
    % Preliminaries
    dt1=zeros(1,fn);
    dt2=zeros(1,spidn-fn);
    % Scale for distance
    scale=(-lb(1)+ub(1));
    % Start looking for any stronger vibration
    for i=1:fn          % Move the females
        for j=1:fn      % Check all the spiders
            if spmass(j)>spmass(i)  % If theres someone more atrractive
                % Calculate the distance
%                 dt1(j)=sqrt(sum((fsp(i,:)-fsp(j,:)).^2));
                dt1(j)=norm(fsp(i,:)-fsp(j,:));
            else
                dt1(j)=0;   % Make sure the value is zero
            end
        end
        for j=1:spidn-fn
            if spmass(fn+j)>spmass(i)
                % Calculate the distance
                dt2(j)=norm(fsp(i,:)-msp(j,:));
                %sqrt(sum((fsp(i,:)-msp(j,:)).^2));
            else
                dt2(j)=0;   % Make sure the value is zero
            end
        end
        % Scaled Distance for the system
        dt=[dt1 dt2]./scale;%%%%&%%
        % Choose the shortest distance
        [~,Ind,val] = find(dt);    % Choose where the distance is non zero
        [~,Imin] = min(val);       % Get the shortest distance
        Ish=Ind(Imin);              % Index of the shortest distance
        %% Check if male or female
        if Ish > fn
        % Is Male
            spaux=msp(Ish-fn,:);   % Assign the shortest distance to spaux
        else
            % Is Female
            spaux=fsp(Ish,:);   % Assign the shortest distance to spaux
        end
        % Calculate the Vibrations
        if isempty(val)     % Check if is the same element
           Vibs=0;          % Vib for the shortest
           spaux=zeros(1,d);
        else
            Vibs=2*(spmass(Ish)*exp(-(rand*dt(Ish).^2)));   % Vib for the shortest
        end
        %% Check if male or female
        if Ibe > fn
            % Is Male
            dt2=norm(fsp(i,:)-msp(Ibe-fn,:));
        else
            % Is Female
            % Vibration of the best
            dt2=norm(fsp(i,:)-fsp(Ibe,:));
        end
        dtb=dt2./scale;
        Vibb=2*(spmass(Ibe)*exp(-(rand*dtb.^2)));
        if rand>=pm
            % Do an atracction
            betha = rand(1,d);
            gamma = rand(1,d);
            tmpf = 2*pm.*(rand(1,d)-0.5);
            fsp(i,:)=fsp(i,:)+(Vibs*(spaux-fsp(i,:)).*betha)+(Vibb*(spbest-fsp(i,:)).*gamma)+tmpf;
        else
            % Do a repulsion
            betha = rand(1,d);
            gamma = rand(1,d);
            tmpf = 2*pm.*(rand(1,d)-0.5);
            fsp(i,:)=fsp(i,:)-(Vibs*(spaux-fsp(i,:)).*betha)-(Vibb*(spbest-fsp(i,:)).*gamma)+tmpf;
        end
    end
    % Check limits
    for i=1:d
        for j=1:fn
            if fsp(j,i)<lb(i), fsp(j,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
            if fsp(j,i)==lb(i), fsp(j,i)=lb(i); end

            if fsp(j,i)>ub(i), fsp(j,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
            if fsp(j,i)==ub(i), fsp(j,i)=ub(i); end
        end
    end
end

function [msp] = MaMove(fn,mn,fsp,msp,femass,mamass,d,lb,ub,pm)
%MAMOVE Summary of this function goes here
%   Detailed explanation goes here
    %Preliminaries
    dt=zeros(1,mn);
    % Scale for distance
    scale=(-lb(1)+ub(1));%/2;
    [Indb,~] = find(mamass>=median(mamass));  % Male spiders above median
    for i=1:mn
        if ismember(i,Indb)     % Spider above the median
            % Start looking for a female with stronger vibration
            for j=1:fn
                if femass(j)>mamass(i)
                    % Calculate the distance
                    dt(j)=norm(msp(i,:)-fsp(j,:));
                else
                    dt(j)=0;
                end
            end
            % Choose the shortest distance
            [~,Ind,val] = find(dt);   % Choose where the distance in non zero
            [~,Imin] = min(val);      % Get the shortest distance
            Ish = Ind(Imin);
            % Update moves
            if isempty(val)
                Vib=0;
                spaux=zeros(1,d);
            else
                dt=dt./scale;
                Vib = 2*femass(Ish)*exp(-(rand*dt(Ish).^2));
                spaux=fsp(Ish,:);
            end
            delta = 2*rand(1,d)-.5;
            tmpf = 2*pm.*(rand(1,d)-0.5);
            msp(i,:) = msp(i,:)+Vib*(spaux-msp(i,:)).*delta+tmpf;
        else % de aqui para abajo falta
            %% Spider below median, go to weigthed mean
            % Generate the weighted mean
            spdpos = [fsp' msp']';
            spdwei = [femass' mamass']';
            weigth = repmat(spdwei,1,d);
            dim = find(size(spdpos)~=1,1);
            wmean = sum(weigth.*spdpos,dim)./sum(weigth,dim);
            %% Move
            delta = 2*rand(1,d)-.5;
            tmpf = 2*pm.*(rand(1,d)-0.5);
            msp(i,:) = msp(i,:)+(wmean-msp(i,:)).*delta+tmpf;
        end
    end
    % Check limits
    for i=1:d
        for j=1:mn
            if msp(j,i)<lb(i), msp(j,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
            if msp(j,i)==lb(i), msp(j,i)=lb(i); end

            if msp(j,i)>ub(i), msp(j,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
            if msp(j,i)==ub(i), msp(j,i)=ub(i); end
        end
    end
end

function [ofsp] = Mating(femass,mamass,fsp,msp,dims)
%CROSSOVER Summary of this function goes here
%   Detailed explanation goes here
    % Generate the offsprings
    ofsp = [];
    cont = 1;
    % Check whether a spider is good or not (above median)
    [Indf,~] = find(femass);                % Female spiders
    [Indm,~] = find(mamass>median(mamass)); % Male spiders above median
    fespid=fsp(Indf,:);             % Female spiders
    maspid=msp(Indm,:);             % Only the Male spiders above median
    sp2mate = [];
    %% Calculate the radio
	rad = zeros(1,dims);
    spid = [fsp' msp']';
    for i=1:dims
        rad(i) = max(spid(:,i))-min(spid(:,i));
    end
    r=(sum(rad)/2)/(dims);
    %% Start looking if there's a good female near
    [sz, ~] = size(Indf);
    dist=zeros(1,sz);
    for i=1:size(Indm)
        iaux = 1;       % Aux to form the elements to mate
        for j=1:size(Indf)
            dist(j)=norm(msp(Indm(i),:)-fsp(Indf(j),:));
        end
        for k=1:size(Indf)
            if dist(k)<r
                mate(iaux,:) = fsp(Indf(k),:);
                mass(iaux) = femass(Indf(k));
                iaux = iaux+1;
            % Form the matrix with elements to mate
            sp2mate = [msp(Indm(i),:)' mate']';
            masmate = [mamass(Indm(i)) mass];                
            end
        end
        % Realizo el mate
        if isempty(sp2mate)
            % do nothing
        else
            [num2,n] = size(sp2mate);
            for k=1:num2
            for j=1:n
                accumulation = cumsum(masmate);
                p = rand() * accumulation(end);
                chosen_index = -1;
                for index = 1 : length(accumulation)
                    if (accumulation(index) > p)
                        chosen_index = index;
                        break;
                    end
                end
                choice = chosen_index;
                % Forma the new element
                ofsp(k,j)=sp2mate(choice,j);
            end
            end
            cont = cont+1;
        end
    end
end

function [fsp,msp,fefit,mafit] = Survive(fsp,msp,ofspr,fefit,mafit,spfit,fun,fn,dims)
%SURVIVE Summary of this function goes here
%   Detailed explanation goes here
    [n1, ~] = size(ofspr);
    %Evalute the offspring
    for j=1:n1
        offit(j)=fun(ofspr(j,:));
    end
    for i=1:n1
        %Calculate the worst spider
        [w1, w2]=max(spfit);
        %If the offspring is better than the worst spider
        if offit(i)<w1
            %Check if is male or female
            if w2>fn
                %Male
                msp(w2-fn,:)=ofspr(i,:);
                mafit(w2-fn)=offit(i);
            else
                %Female
                fsp(w2,:)=ofspr(i,:);
                fefit(w2)=offit(i);
            end
            spfit(w2)=offit(i);
        end  
    end
end