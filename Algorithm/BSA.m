function Data_o = BSA(Data_i,Data_o)
t_c=clock;                                                              
M = Data_i.maxIter;       
pop = Data_i.pop;         
dim = Data_i.dim;               
x = Data_i.X;                   
fit = Data_i.F_value;     
lb= Data_i.lb;   
ub= Data_i.ub;    
FQ = 10;   
c1 = 1.5;
c2 = 1.5;
a1 = 1;
a2 = 1;

pFit = fit; % The individual's best fitness value
pX = x;     % The individual's best position corresponding to the pFit
CNVG=inf*zeros(1,Data_i.maxIter);
[ fMin, bestIndex ] = min( fit );  % fMin denotes the global optimum
% bestX denotes the position corresponding to fMin
bestX = x( bestIndex, : );   
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fMin;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=bestIndex;

 for iteration = 1 : M
     
    prob = rand( pop, 1 ) .* 0.2 + 0.8;%The probability of foraging for food
    
    if( mod( iteration, FQ ) ~= 0 )         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Birds forage for food or keep vigilance
        sumPfit = sum( pFit );
        meanP = mean( pX );
        for i = 1 : pop
            if rand < prob(i)
                x( i, : ) = x( i, : ) + c1 * rand.*(bestX - x( i, : ))+ ...
                    c2 * rand.*( pX(i,:) - x( i, : ) );
            else
                person = randiTabu( 1, pop, i, 1 );
                
                x( i, : ) = x( i, : ) + rand.*(meanP - x( i, : )) * a1 * ...
                    exp( -pFit(i)/( sumPfit + realmin) * pop ) + a2 * ...
                    ( rand*2 - 1) .* ( pX(person,:) - x( i, : ) ) * exp( ...
                    -(pFit(person) - pFit(i))/(abs( pFit(person)-pFit(i) )...
                    + realmin) * pFit(person)/(sumPfit + realmin) * pop ); 
            end
            
            x( i, : ) = Bounds( x( i, : ), lb, ub );  
            fit( i ) = Data_i.fobj( x( i, : ) );
        end
        
    
    else
        FL = rand( pop, 1 ) .* 0.4 + 0.5;    
        
       
        [~, minIndex ] = min( pFit );
        [~, maxIndex ] = max( pFit );
        choose = 0;
        if ( minIndex < 0.5*pop && maxIndex < 0.5*pop )
            choose = 1;
        end
        if ( minIndex > 0.5*pop && maxIndex < 0.5*pop )
            choose = 2;
        end
        if ( minIndex < 0.5*pop && maxIndex > 0.5*pop )
            choose = 3;
        end
        if ( minIndex > 0.5*pop && maxIndex > 0.5*pop )
            choose = 4;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if choose < 3
            for i = (pop/2+1) : pop
                x( i, : ) = x( i, : ) * ( 1 + randn );
                x( i, : ) = Bounds( x( i, : ), lb, ub );
                fit( i ) = Data_i.fobj( x( i, : ) );
            end
            if choose == 1 
                x( minIndex,: ) = x( minIndex,: ) * ( 1 + randn );
                x( minIndex, : ) = Bounds( x( minIndex, : ), lb, ub );
                fit( minIndex ) = Data_i.fobj( x( minIndex, : ) );
            end
            for i = 1 : 0.5*pop
                if choose == 2 || minIndex ~= i
                    person = randi( [(0.5*pop+1), pop ], 1 );
                    x( i, : ) = x( i, : ) + (pX(person, :) - x( i, : )) * FL( i );
                    x( i, : ) = Bounds( x( i, : ), lb, ub );
                    fit( i ) = Data_i.fobj( x( i, : ) );
                end
            end
        else
            for i = 1 : 0.5*pop
                x( i, : ) = x( i, : ) * ( 1 + randn );
                x( i, : ) = Bounds( x( i, : ), lb, ub );
                fit( i ) = Data_i.fobj( x( i, : ) );
            end
            if choose == 4 
                x( minIndex,: ) = x( minIndex,: ) * ( 1 + randn );
                x( minIndex, : ) = Bounds( x( minIndex, : ), lb, ub );
                fit( minIndex ) = Data_i.fobj( x( minIndex, : ) );
            end
            for i = (0.5*pop+1) : pop
                if choose == 3 || minIndex ~= i
                    person = randi( [1, 0.5*pop], 1 );
                    x( i, : ) = x( i, : ) + (pX(person, :) - x( i, : )) * FL( i );
                    x( i, : ) = Bounds( x( i, : ), lb, ub );
                    fit( i ) = Data_i.fobj( x( i, : ) );
                end
            end   
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the individual's best fitness vlaue and the global best one
   
    for i = 1 : pop 
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )
            fMin = pFit( i );
            bestX = pX( i, : );
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    CNVG(iteration) = fMin;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fMin;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iteration;                    
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=CNVG;
end
% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  
  % Apply the upper bound vector 
  J = temp > Ub;
  temp(J) = Ub(J);
  % Update this new move 
  s = temp;
end
%--------------------------------------------------------------------------
% This function generate "dim" values, all of which are different from
%  the value of "tabu"
function value = randiTabu( min, max, tabu, dim )
    value = ones( dim, 1 ) .* max .* 2;
    num = 1;
    while ( num <= dim )
        temp = randi( [min, max], 1, 1 );
        if( length( find( value ~= temp ) ) == dim && temp ~= tabu )
            value( num ) = temp;
            num = num + 1;
        end
    end
end