function Data_o = CSO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
x = Data_i.X;                        
fit=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
FitFunc = Data_i.fobj;      M = Data_i.maxIter;         pop = Data_i.pop;       dim = Data_i.dim;  
lb=Data_i.lb;       ub=Data_i.ub;
G = 10;       
rPercent = 0.15;        hPercent = 0.7;         mPercent = 0.5;
pFit = fit; % The individual's best fitness value
pX = x;
rNum = round( pop * rPercent );    % The population size of roosters
hNum = round( pop * hPercent );    % The population size of hens
cNum = pop - rNum - hNum;          % The population size of chicks
mNum = round( hNum * mPercent );   % The population size of mother hens

[ fMin, bestIndex ] = min( fit );
bestX = x( bestIndex, : ); 
for t = 1 : M
    FL = rand( pop, 1 ) .* 0.4 + 0.5;  
    if mod( t, G ) == 1 || t == 1
        [ ans, sortIndex ] = sort( pFit );
        motherLib = randperm( hNum, mNum ) + rNum;
        mate = randpermF( rNum, hNum );
        mother = motherLib( randi( mNum, cNum, 1 ) );
    end   
    for i = 1 : rNum
        anotherRooster = randiTabu( 1, rNum, i, 1 );  
        if( pFit( sortIndex( i ) ) <= pFit( sortIndex( anotherRooster ) ) )
            tempSigma = 1;
        else
            tempSigma = exp( ( pFit( sortIndex( anotherRooster ) ) - ...
                pFit( sortIndex( i ) ) ) / ( abs( pFit( sortIndex(i) ) )...
                + realmin ) );
        end
        x( sortIndex( i ), : ) = pX( sortIndex( i ), : ) .* ( 1 + ...
            tempSigma .* randn( 1, dim ) );
        x( sortIndex( i ), : ) = BoundaryCheck( x( sortIndex( i ), : ), ub, lb,dim );
        fit( sortIndex( i ) ) = FitFunc( x( sortIndex( i ), : ) );
    end    

    for i = ( rNum + 1 ) : ( rNum + hNum )
        other = randiTabu( 1,  i,  mate( i - rNum ), 1 );
        c1 = exp( ( pFit( sortIndex( i ) ) - pFit( sortIndex( mate( i - ...
            rNum ) ) ) )/ ( abs( pFit( sortIndex(i) ) ) + realmin ) );
            
        c2 = exp( ( -pFit( sortIndex( i ) ) + pFit( sortIndex( other ) )));

        x( sortIndex( i ), : ) = pX( sortIndex( i ), : ) + ( pX(...
            sortIndex( mate( i - rNum ) ), : )- pX( sortIndex( i ), : ) )...
             .* c1 .* rand( 1, dim ) + ( pX( sortIndex( other ), : ) - ...
             pX( sortIndex( i ), : ) ) .* c2 .* rand( 1, dim ); 
        x( sortIndex( i ), : ) = BoundaryCheck( x( sortIndex( i ), : ), ub, lb,dim );
        fit( sortIndex( i ) ) = FitFunc( x( sortIndex( i ), : ) );
    end

    for i = ( rNum + hNum + 1 ) : pop    % Update the cNum chicks' values.
        x( sortIndex( i ), : ) = pX( sortIndex( i ), : ) + ( pX( ...
            sortIndex( mother( i - rNum - hNum ) ), : ) - ...
            pX( sortIndex( i ), : ) ) .* FL( i );
        x( sortIndex( i ), : ) = BoundaryCheck( x( sortIndex( i ), : ), ub, lb,dim );
        fit( sortIndex( i ) ) = FitFunc( x( sortIndex( i ), : ) );
    end

    for i = 1 : pop 
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )
            fMin = pFit( i );
            bestX = pX( i, : );
            if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>fMin
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=fMin;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
    end
    IterCurve(t)=fMin;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

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
%--------------------------------------------------------------------------
function result = randpermF( range, dim )
% The original function "randperm" in Matlab is only confined to the
% situation that dimension is no bigger than dim. This function is 
% applied to solve that situation.

temp = randperm( range, range );
temp2 = randi( range, dim, 1 );
index = randperm( dim, ( dim - range ) );
result = [ temp, temp2( index )' ];
end