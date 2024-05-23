function Data_o = CSAE(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
chameleonPositions = Data_i.X;                        
fit=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
fitness=fit;
[fmin0,index]=min(fit);
chameleonBestPosition = chameleonPositions; % Best position initialization
gPosition = chameleonPositions(index,:); % initial global position
fobj=Data_i.fobj;   iteMax=Data_i.maxIter;   searchAgents=Data_i.pop;  dim=Data_i.dim; ub=Data_i.ub;  lb=Data_i.lb;
v=0.1*chameleonBestPosition;% initial velocity
v0=0.0*v;
rho=1.0;
p1=2.0;  
p2=2.0;  
c1=2.0; 
c2=1.80;  
gamma=2.0; 
alpha = 4.0;  
beta=3.0; 
for t=1:Data_i.maxIter
    a = 2590*(1-exp(-log(t))); 
    omega=(1-(t/iteMax))^(rho*sqrt(t/iteMax)) ; 
    p1 = 2* exp(-2*(t/iteMax)^2);  % 
    p2 = 2/(1+exp((-t+iteMax/2)/100)) ;       
    mu= gamma*exp(-(alpha*t/iteMax)^beta);
    ch=ceil(searchAgents*rand(1,searchAgents));
    for i=1:searchAgents  
        if rand>=0.1
            chameleonPositions(i,:)= chameleonPositions(i,:)+ p1*(chameleonBestPosition(ch(i),:)-chameleonPositions(i,:))*rand()+... 
            + p2*(gPosition -chameleonPositions(i,:))*rand();
        else 
            for j=1:dim
                chameleonPositions(i,j)= gPosition(j)+mu*((ub(j)-lb(j))*rand+lb(j))*sign(rand-0.50) ;
            end 
        end   
    end
    for i=1:searchAgents   
        v(i,:)= omega*v(i,:)+ p1*(chameleonBestPosition(i,:)-chameleonPositions(i,:))*rand +.... 
               + p2*(gPosition-chameleonPositions(i,:))*rand;        
         chameleonPositions(i,:)=chameleonPositions(i,:)+(v(i,:).^2 - v0(i,:).^2)/(2*a);
    end
    v0=v;
    for i=1:searchAgents
         if chameleonPositions(i,:)<lb
            chameleonPositions(i,:)=lb;
         elseif chameleonPositions(i,:)>ub
                chameleonPositions(i,:)=ub;
         end
    end
    for i=1:searchAgents
        chameleonPositions(i,:)=BoundaryCheck(chameleonPositions(i,:),Data_i.ub,Data_i.lb,Data_i.dim);  %%%%%*2
        fit(i,1)=fobj(chameleonPositions(i,:)) ;
        if fit(i)<fitness(i)
            chameleonBestPosition(i,:) = chameleonPositions(i,:); % Update the best positions  
            fitness(i)=fit(i); % Update the fitness
        end
    end
    [fmin,index]=min(fitness);
    if fmin < fmin0
        gPosition = chameleonBestPosition(index,:); % Update the global best positions
        fmin0 = fmin;
        if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fmin
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fmin;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
        end
    end
    IterCurve(t)=fmin0;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function answer = get_orthonormal(m,n)
    if ( (nargin==2) && (m>n) && (isnumeric(m)*isnumeric(n)) )
    elseif ( nargin==1 && isnumeric(m) && length(m)==1 )
        n=m;  
    else
        error('Incorrect Inputs. Please read help text in m-file.')
    end
    count=0;
    while (count==0)
        [P,D] = eig(B) ;
        if ((P'*P - eye(m))>eps) 
            % error, vectors not orthonormal, repeat the random matrix draw again
            count=0;
        else
            % we want the first n of these orthonormal columns
            answer=P(:,1:n) ;
            count=1;
        end
    end
end


function [chameleonPositions]=rotation(chameleonPosition, searchAgents, dim)
for i=1:searchAgents      
      if (dim>2) 
          xmax=1;xmin=-1;
          th=round(xmin+rand(1,1)*(xmax-xmin));
          vec=get_orthonormal(dim,2);
          vecA=vec(:,1);
          vecB=vec(:,2);
          theta=(th*rand()*180)*(pi/180) ;
          Rot = RotMatrix(theta,vecA, vecB) ;
         if (theta~=0)
            V=[chameleonPosition(i,:) ]; 
            V_centre=mean(V,1); %Centre, of line
            Vc=V-ones(size(V,1),1)*V_centre; %Centering coordinates
            Vrc=[Rot*Vc']'; %Rotating centred coordinates
            Vr=Vrc+ones(size(V,1),1)*V_centre; %Shifting back to original location
            chameleonPosition(i,:)=((Vr)/1); 
         end
     else
          xmax=1;xmin=-1;
          th=round(xmin+rand(1,1)*(xmax-xmin));
          theta=th*rand()*180*(pi/180);
          Rot = RotMatrix(theta);
           if (theta~=0)
            V=[chameleonPosition(i,:) ];  
            V_centre=mean(V,1); %Centre, of line
            Vc=V-ones(size(V,1),1)*V_centre; %Centering coordinates
            Vrc=[Rot*Vc']'; %Rotating centred coordinates
            Vr=Vrc+ones(size(V,1),1)*V_centre; %Shifting back to original location
            chameleonPosition(i,:)=((Vr)/1);
           end
      end
end  
    chameleonPositions=chameleonPosition;
end

function R = RotMatrix(alpha, u, v)
    if numel(alpha) ~= 1
   error('JSimon:RotMatrrix:BadInput1', ...
      'Angle of rotation must be a scalar.');
    s = sin(alpha);
    c = cos(alpha);
    switch nargin
   case 1
      % 2D rotation matrix:
      R = [c, -s;  s, c];
   case 2
      if numel(u) ~= 3
         error('JSimon:RotMatrrix:BadAxis2D', ...
            '3D: Rotation axis must have 3 elements.');
      end
      u = u(:);
      u = u ./ sqrt(u.' * u);
      x  = u(1);
      y  = u(2);
      z  = u(3);
      mc = 1 - c;
      R  = [c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s; ...
            x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
            x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];       
   case 3
      n = numel(u);
      if n ~= numel(v)
         error('JSimon:RotMatrrix:BadAxes3D', ...
            'ND: Axes to define plane of rotation must have the same size.');
      end
      u = u(:);
      u = u ./ sqrt(u.' * u);
      v = v(:);
      v = v - (u.' * v) * u;
      v = v ./ sqrt(v.' * v);
      R = eye(n) + ...
         (v * u.' - u * v.') * s + ...
         (u * u.' + v * v.') * (c - 1);
   otherwise
      error('JSimon:RotMatrrix:BadNInput', ...
            '1 to 3 inputs required.');
    end

    end
end


