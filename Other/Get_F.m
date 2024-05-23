function [LB,UB,Dim,F_obj] = Get_F(F)
F_obj=str2func(F);
switch F
    case 'F1'
        LB=-100;
        UB=100;
        Dim =30;             
    case 'F2'
        LB=-10;
        UB=10;
        Dim = 30;
    case 'F3'
        LB=-100;
        UB=100;
        Dim = 30;
    case 'F4'
        LB=-100;
        UB=100;
        Dim = 30;
    case 'F5'
        LB=-30;
        UB=30;
        Dim = 30;
    case 'F6'
        LB=-100;
        UB=100;
        Dim = 30;
    case 'F7'
        LB=-1.28;
        UB=1.28;
        Dim = 30;
    case 'F8'
       LB=-500;
        UB=500;
        Dim = 30;
    case 'F9'
        LB=-5.12;
        UB=5.12;
        Dim = 30;
    case 'F10'
        LB=-32;
        UB=32;
        Dim = 30;
    case 'F11'
        LB=-600;
        UB=600;
        Dim = 30;
    case 'F12'
        LB=-50;
        UB=50;
        Dim = 30;
    case 'F13'
        LB=-50;
        UB=50;
        Dim = 30;
    case 'F14'
        LB=-65;
        UB=65;
        Dim=2;
    case 'F15'
        LB=-5;
        UB=5;
        Dim=4;
    case 'F16'
        LB=-5;
        UB=5;
        Dim=2;
    case 'F17'
        LB=[-5,0];
        UB=[10,15];
        Dim=2;
    case 'F18'
        LB=-2;
        UB=2;
        Dim=2;
    case 'F19'
        LB=1;
        UB=3;
        Dim=3;
    case 'F20'
        LB=0;
        UB=1;
        Dim=6;     
    case 'F21'
        LB=0;
        UB=10;
        Dim=4;    
    case 'F22'
        LB=0;
        UB=10;
        Dim=4;    
    case 'F23'
        LB=0;
        UB=10;
        Dim=4;
    case 'F24'   
        LB=-6.4;
        UB=6.35;
        Dim=6;
    case 'F25'  
        LB=[0,0,0,-4,-4-(1/4)*floor(1/3),-4-(1/4)*floor(2/3),-4-(1/4)*floor(3/3),...
            -4-(1/4)*floor(4/3),-4-(1/4)*floor(5/3),-4-(1/4)*floor(6/3),-4-(1/4)*floor(7/3)...
            ,-4-(1/4)*floor(8/3),-4-(1/4)*floor(9/3),-4-(1/4)*floor(10/3),-4-(1/4)*floor(11/3)...
            ,-4-(1/4)*floor(12/3),-4-(1/4)*floor(13/3),-4-(1/4)*floor(14/3),-4-(1/4)*floor(15/3)...
            ,-4-(1/4)*floor(16/3),-4-(1/4)*floor(17/3),-4-(1/4)*floor(18/3),-4-(1/4)*floor(19/3)...
            ,-4-(1/4)*floor(20/3),-4-(1/4)*floor(21/3),-4-(1/4)*floor(22/3),-4-(1/4)*floor(23/3)...
            ,-4-(1/4)*floor(24/3),-4-(1/4)*floor(25/3),-4-(1/4)*floor(26/3)];
        UB=[4,4,pi,4,4+(1/4)*floor(1/3),4+(1/4)*floor(2/3),4+(1/4)*floor(3/3),4+(1/4)*floor(4/3)...
            ,4+(1/4)*floor(5/3),4+(1/4)*floor(6/3),4+(1/4)*floor(7/3),4+(1/4)*floor(8/3)...
            ,4+(1/4)*floor(9/3),4+(1/4)*floor(10/3),4+(1/4)*floor(11/3),4+(1/4)*floor(12/3)...
            ,4+(1/4)*floor(13/3),4+(1/4)*floor(14/3),4+(1/4)*floor(15/3),4+(1/4)*floor(16/3)...
            ,4+(1/4)*floor(17/3),4+(1/4)*floor(18/3),4+(1/4)*floor(19/3),4+(1/4)*floor(20/3)...
            ,4+(1/4)*floor(21/3),4+(1/4)*floor(22/3),4+(1/4)*floor(23/3),4+(1/4)*floor(24/3)...
            ,4+(1/4)*floor(25/3),4+(1/4)*floor(26/3)];
        Dim=30;
    case 'F26'  
        LB=0.6;   
        UB=0.9;
        Dim=1;
    case 'F27'  
        LB=0;
        UB=5;
        Dim=1;
    case 'F28' 
        LB=[0,0,0,-4,-4-(1/4)*floor(1/3),-4-(1/4)*floor(2/3),-4-(1/4)*floor(3/3),...
            -4-(1/4)*floor(4/3),-4-(1/4)*floor(5/3),-4-(1/4)*floor(6/3),-4-(1/4)*floor(7/3)...
            ,-4-(1/4)*floor(8/3),-4-(1/4)*floor(9/3),-4-(1/4)*floor(10/3),-4-(1/4)*floor(11/3)...
            ,-4-(1/4)*floor(12/3),-4-(1/4)*floor(13/3),-4-(1/4)*floor(14/3),-4-(1/4)*floor(15/3)...
            ,-4-(1/4)*floor(16/3),-4-(1/4)*floor(17/3),-4-(1/4)*floor(18/3),-4-(1/4)*floor(19/3)...
            ,-4-(1/4)*floor(20/3),-4-(1/4)*floor(21/3),-4-(1/4)*floor(22/3),-4-(1/4)*floor(23/3)...
            ,-4-(1/4)*floor(24/3),-4-(1/4)*floor(25/3),-4-(1/4)*floor(26/3)];
        UB=[4,4,pi,4,4+(1/4)*floor(1/3),4+(1/4)*floor(2/3),4+(1/4)*floor(3/3),4+(1/4)*floor(4/3)...
            ,4+(1/4)*floor(5/3),4+(1/4)*floor(6/3),4+(1/4)*floor(7/3),4+(1/4)*floor(8/3)...
            ,4+(1/4)*floor(9/3),4+(1/4)*floor(10/3),4+(1/4)*floor(11/3),4+(1/4)*floor(12/3)...
            ,4+(1/4)*floor(13/3),4+(1/4)*floor(14/3),4+(1/4)*floor(15/3),4+(1/4)*floor(16/3)...
            ,4+(1/4)*floor(17/3),4+(1/4)*floor(18/3),4+(1/4)*floor(19/3),4+(1/4)*floor(20/3)...
            ,4+(1/4)*floor(21/3),4+(1/4)*floor(22/3),4+(1/4)*floor(23/3),4+(1/4)*floor(24/3)...
            ,4+(1/4)*floor(25/3),4+(1/4)*floor(26/3)];
        Dim=30;
    case 'F29'  
        LB=[0,0,0,-4,-4-(1/4)*floor(1/3),-4-(1/4)*floor(2/3),-4-(1/4)*floor(3/3),...
            -4-(1/4)*floor(4/3),-4-(1/4)*floor(5/3),-4-(1/4)*floor(6/3),-4-(1/4)*floor(7/3)...
            ,-4-(1/4)*floor(8/3),-4-(1/4)*floor(9/3),-4-(1/4)*floor(10/3),-4-(1/4)*floor(11/3)...
            ,-4-(1/4)*floor(12/3),-4-(1/4)*floor(13/3),-4-(1/4)*floor(14/3),-4-(1/4)*floor(15/3)...
            ,-4-(1/4)*floor(16/3),-4-(1/4)*floor(17/3),-4-(1/4)*floor(18/3),-4-(1/4)*floor(19/3)...
            ,-4-(1/4)*floor(20/3),-4-(1/4)*floor(21/3),-4-(1/4)*floor(22/3),-4-(1/4)*floor(23/3)...
            ,-4-(1/4)*floor(24/3),-4-(1/4)*floor(25/3),-4-(1/4)*floor(26/3)];
        UB=[4,4,pi,4,4+(1/4)*floor(1/3),4+(1/4)*floor(2/3),4+(1/4)*floor(3/3),4+(1/4)*floor(4/3)...
            ,4+(1/4)*floor(5/3),4+(1/4)*floor(6/3),4+(1/4)*floor(7/3),4+(1/4)*floor(8/3)...
            ,4+(1/4)*floor(9/3),4+(1/4)*floor(10/3),4+(1/4)*floor(11/3),4+(1/4)*floor(12/3)...
            ,4+(1/4)*floor(13/3),4+(1/4)*floor(14/3),4+(1/4)*floor(15/3),4+(1/4)*floor(16/3)...
            ,4+(1/4)*floor(17/3),4+(1/4)*floor(18/3),4+(1/4)*floor(19/3),4+(1/4)*floor(20/3)...
            ,4+(1/4)*floor(21/3),4+(1/4)*floor(22/3),4+(1/4)*floor(23/3),4+(1/4)*floor(24/3)...
            ,4+(1/4)*floor(25/3),4+(1/4)*floor(26/3)];
        Dim=30;  
    case 'F30'  
        LB=0;
        UB=2*pi;
        Dim=20;
    case 'F31'  
        LB=1;   
        UB=15;
        Dim=7;
    case 'F32'  
        LB=zeros(1,126);
        UB=GD_max();
        Dim=126;
    case 'F33'  
        LB=[0.2,0.2,0.2,0.2,0.2,0.2,-180,-180,-180,-180,-180,-180];
        UB=[1,1,1,1,1,1,180,180,180,180,180,180];
        Dim=12;
    case 'F34' 
        min=[10,20,30,40,50];           max=[75,125,175,250,300];
        LB=repmat(min,1,24);
        UB=repmat(max,1,24);
        Dim=120;
    case 'F35'  
        min=[150,135,73,60,73,57,20,47,20];     max=[470,460,340,300,243,160,130,120,80];
        LB=repmat(min,1,24);
        UB=repmat(max,1,24);
        Dim=216;
        F_obj=str2func('F34');          
    case 'F36'  
        LB=[100,50,80,50,50,50];
        UB=[500,200,300,150,200,120];
        Dim=6;
    case 'F37'  
        LB=[0,0,0,60,60,60,60,60,60,40,40,55,55];
        UB=[680,360,360,180,180,180,180,180,180,120,120,120,120];
        Dim=13;
        F_obj=str2func("F36");
    case 'F38'  
        LB=[150,150,20,20,150,135,135,60,25,25,20,20,25,15,15];
        UB=[455,455,130,130,470,460,465,300,162,160,80,80,85,55,55];
        Dim=15;
        F_obj=str2func("F36");
    case 'F39'  
        LB=[36,36,60,80,47,68,110,135,135,130,94,94,125,125,125,125....
            ,220,220,242,242,254,254,254,254,254,254,10,10,10,47,60....
            ,60,60,90,90,90,25,25,25,242];
        UB=[114,114,120,190,97,140,300,300,300,300,375,375,500,500....
            ,500,500,500,500,550,550,550,550,550,550,550,550,150,....
            150,150,97,190,190,190,200,200,200,110,110,110,550];
        Dim=40;
        F_obj=str2func("F36");
    case 'F40'  
        LB=[71,120,125,125,90,90,280,280,260,260,260,260,260,260,260....
            ,260,260,260,260,260,260,60,260,260,280,280,280,280,260,260....
            ,260,260,260,260,260,260,120,120,423,423,3,3,160,160,160,160....
            ,160,160,160,160,165,165,165,165,180,180,103,198,100,153,163....
            ,95,160,160,196,196,196,196,130,130,137,137,195,175,175,175,175....
            ,330,160,160,200,56,115,115,115,207,207,175,175,175,175,360,415....
            ,795,795,578,615,612,612,758,755,750,750,713,718,791,786,795,795....
            ,795,795,94,94,94,244,244,244,95,95,116,175,2,4,15,9,12,10,112....
            ,4,5,5,50,5,42,42,41,17,7,7,26];
        UB=[119,189,190,190,190,190,490,490,496,496,496,496,506,509,506,505....
            ,506,506,505,505,505,505,505,505,537,537,549,549,501,501,506,506....
            ,506,506,500,500,241,241,774,769,19,28,250,250,250,250,250,250....
            ,250,250,504,504,504,504,471,561,341,617,312,471,500,302,511,511....
            ,490,490,490,490,432,432,455,455,541,536,540,538,540,574,531,531....
            ,542,132,245,245,245,307,307,345,345,345,345,580,645,984,978,682,720....
            ,718,720,964,958,1007,1006,1013,1020,954,952,1006,1013,1021,1015....
            ,203,203,203,379,379,379,190,189,194,321,19,59,83,53,37,34,373,20....
            ,38,19,98,10,74,74,105,51,19,19,40];
        Dim=140;
        F_obj=str2func("F36");
    case 'F41'  
        min=[5,6,10,13];     max=[15,15,30,25];
        LB=repmat(min,1,24);
        UB=repmat(max,1,24);
        Dim=96;
    case 'F42' 
        min=[5,6,10,13];     max=[15,15,30,25];
        LB=repmat(min,1,24);
        UB=repmat(max,1,24);
        Dim=96;
    case 'F43' 
        min=[5,6,10,13];     max=[15,15,30,25];
        LB=repmat(min,1,24);
        UB=repmat(max,1,24);
        Dim=96;
    case 'F44' 
        LB=[1900,2.5,0,0,100,100,100,100,100,100,0.01,0.01,0.01,0.01,0.01,0.01...
            ,1.1,1.1,1.05,1.05,1.05,-pi,-pi,-pi,-pi,-pi];
        UB=[2300,4.05,1,1,500,500,500,500,500,600,0.99,0.99,0.99,0.99,0.99,0.99...
            ,6,6,6,6,6,pi,pi,pi,pi,pi];
        Dim=26;
    case 'F45' 
        LB=[-1000,3,0,0,100,100,30,400,800,0.01,0.01,0.01,0.01,0.01,1.05,1.05...
            ,1.15,1.7,-pi,-pi,-pi,-pi];
        UB=[0,5,1,1,400,500,300,1600,2200,0.9,0.9,0.9,0.9,0.9,6,6,6.5,291,pi...
            ,pi,pi,pi];
        Dim=22;
    %----------------CEC2019----------------
    case 'F46'                  
        LB=-8192;
        UB=8192;
        Dim=9;
    case 'F47'                  
        LB=-16384;
        UB=16384;
        Dim=16;
    case 'F48'                  
        LB=-4;
        UB=4;
        Dim=18;
    case 'F49'                  
        LB=-100;
        UB=100;
        Dim=10;
    case 'F50'                  
        LB=-100;
        UB=100;
        Dim=10;
    case 'F51'                  
        LB=-100;
        UB=100;
        Dim=10;
    case 'F52'                  
        LB=-100;
        UB=100;
        Dim=10;
    case 'F53'                  
        LB=-100;
        UB=100;
        Dim=10;
    case 'F54'                  
        LB=-100;
        UB=100;
        Dim=10;
    case 'F55'                  
        LB=-100;
        UB=100;
        Dim=10;
    %----------------CEC2022----------------
    case 'F56'               
        LB=-100;
        UB=100;
        Dim=20;
    case 'F57'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F58'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F59'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F60'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F61'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F62'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F63'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F64'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F65'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F66'           
        LB=-100;
        UB=100;
        Dim=20;
    case 'F67'           
        LB=-100;
        UB=100;
        Dim=20;
    %----------------application----------------
    case 'F68'
        LB=[0.0625,0.0625,10,10];              
        UB=[0.0625*99,0.0625*99,200,200];
        Dim=4;
    case 'F69'                                  
        LB = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5];
        UB = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5];
        Dim = 7;
    case 'F70'                                  
        LB = [0.1, 0.1, 0.1, 0.1];
        UB = [10, 10, 10, 2];
        Dim = 4;
    case 'F71'                                 
        D=160;
        d=90;
        LB=[0.5*(D+d), 0.15*(D-d), 4, 0.515, 0.515, 0.4, 0.6, 0.3, 0.02, 0.6];
        UB=[0.6*(D+d), 0.45*(D-d), 50, 0.6, 0.6, 0.5, 0.7, 0.4, 0.1, 0.85];
        Dim=10;
    case 'F72'                                  
        LB=[0.05, 0.25, 2];   
        UB=[2, 1.3, 15];
        Dim=3;
    case 'F73'                                  
        LB = [12, 12, 12, 12];
        UB= [60, 60, 60, 60];
        Dim=4;
    case 'F74'                                 
        LB=[0,0,0,0,0];
        UB=[60,60,90,90,90];
        Dim=5;
    case 'F75'                                  
        LB = [0, 0];
        UB = [1, 1];
        Dim=2;
    case 'F76'                                  
        LB = [0.01,0.01,0.01,0.01,0.01];
        UB= [100,100,100,100,100];
        Dim=5;
    case 'F77'                                  
        LB=[60,90,1,600,2];
        UB=[80,110,3,1000,9];
        Dim=5;
    otherwise
        disp('数据越界！');
    end
end

% F1

function o = F1(x)
    o=sum(x.^2);
end

% F2

function o = F2(x)
    o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
    dim=size(x,2);
    o=0;
    for i=1:dim
        o=o+sum(x(1:i))^2;
    end
end

% F4

function o = F4(x)
    o=max(abs(x));
end

% F5

function o = F5(x)
    dim=size(x,2);
    o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
    o=sum(abs((x+.5)).^2);
end

% F7

function o = F7(x)
    dim=size(x,2);
    o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
    o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
    dim=size(x,2);
    o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
    dim=size(x,2);
    o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
    dim=size(x,2);
    o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
    dim=size(x,2);
    o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
    (1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end
function o=Ufun(x,a,k,m)
    o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end
% F13

function o = F13(x)
    dim=size(x,2);
    o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
    ((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function o = F14(x)
    aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
    -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
    
    for j=1:25
        bS(j)=sum((x'-aS(:,j)).^6);
    end
    o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
    aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
    bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
    o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16

function o = F16(x)
    o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17

function o = F17(x)
    o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18

function o = F18(x)
    o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
        (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19

function o = F19(x)
    aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
    pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
    o=0;
    for i=1:4
        o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
    end
end

% F20

function o = F20(x)
    aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
    cH=[1 1.2 3 3.2];
    pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
    .2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
    o=0;
    for i=1:4
        o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
    end
end

% F21

function o = F21(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    o=0;
    for i=1:5
        o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

% F22

function o = F22(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    o=0;
    for i=1:7
        o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

% F23

function o = F23(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    o=0;
    for i=1:10
        o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end


function o=F24(x)
    d=length(x);
    if d<6
        disp('dimension-size should be six.')
    else
        if d>6
            disp('dimension-size is more than 6.')
            disp('function has been evaluated on first six dimensions.')
        end
    end
    theta=2*pi/100;
    o=0;
    for t=0:100
        y_t=x(1)*sin(x(2)*t*theta+x(3)*sin(x(4)*t*theta+x(5)*sin(x(6)*t*theta)));
        y_0_t=1*sin(5*t*theta-1.5*sin(4.8*t*theta+2*sin(4.9*t*theta)));
        o=o+(y_t-y_0_t)^2;
    end
end


function o=F25(x)
    r=[];
    p=size(x);
    if rem(p(2),3)~=0
        disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
    end
    n=p(2)/3;
    x=reshape(x,3,n)';
    v=0;
    a=ones(n,n);
    b=repmat(2,n,n);
    for i=1:(n-1)
        for j=(i+1):n
            r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
            v=v+(a(i,j)/r(i,j)^12-b(i,j)/r(i,j)^6);
        end
    end
    o=v;
end

function o = F26(x)
   
    tol=1.0e-01; 
    tspan=[0 0.78];
    yo =[1 0 0 0 0 0 0];
    u=x;          
    options = odeset('RelTol',tol);
    [T,Y] = ode45(@(t,y) diffsolv(t,y,u),tspan,yo,options);
    w=size(Y);
    o=Y(w(1),w(2))*1e+003;
end

function o = F27(x)
        
    tol=1.0e-01;
    tspan=[0 0.78];
    yo =[ 0.09 0.09]';
    u=x;%u should be passed here.
    options = odeset('RelTol',tol);
    [T,Y] = ode45(@(t,y) intgrl(t,y,u),tspan,yo,options);
    o=sum(sum(Y.^2,2)+(0.1)*(u).*(u));
end


function o = F28(x)
    p=size(x);
    if rem(p(2),3)~=0
        disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
    end
    NP=p(2)/3;
    x=reshape(x,3,NP)';
    R1=3.0;
    R2=0.2;
    A=3.2647e+3;
    B=9.5373e+1;
    lemda1=3.2394;
    lemda2=1.3258;
    lemda3=1.3258;
    c=4.8381;
    d=2.0417;
    n1=22.956;
    gama=0.33675;
    h=0;
    E=zeros(1,NP);
    r=zeros(NP);
    for i=1:NP
        for j=1:NP
            r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
            if r(i,j)<(R1-R2)
                fcr(i,j)=1;
            elseif  r(i,j)>(R1+R2)
                fcr(i,j)=0;
            else
                fcr(i,j)=0.5-0.5*sin(pi/2*(r(i,j)-R1)/R2);
            end
            VRr(i,j)=A*exp(-lemda1*r(i,j));
            VAr(i,j)=B*exp(-lemda2*r(i,j));
        end
    end
    for i=1:NP
        for j=1:NP
            if i==j
                continue
            end
            jeta=zeros(NP,NP);
            for k=1:NP
                if i==k || j==k        
                    continue;  
                end
                rd1=sqrt(sum((x(i,:)-x(k,:)).^2));
                rd3=sqrt(sum((x(k,:)-x(j,:)).^2));
                rd2=sqrt(sum((x(i,:)-x(j,:)).^2));
                ctheta_ijk=(rd1^2+rd2^2-rd3^3)/(2*rd1*rd2);
                G_th_ijk =1+(c^2)/(d^2)-(c^2)/(d^2+(h-ctheta_ijk)^2);
                jeta(i,j)=jeta(i,j)+fcr(i,k)*G_th_ijk*exp(lemda3^3*(r(i,j)-r(i,k))^3);
            end
            Bij=(1+(gama*jeta(i,j))^n1)^(-0.5/n1);
            E(i)=E(i)+fcr(i,j)*(VRr(i,j)-Bij*VAr(i,j))/2;
        end    
    end
    o=sum(E);
end


function o = F29(x)
    
    p=size(x);
    if rem(p(2),3)~=0
        disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
    end
    NP=p(2)/3;
    x=reshape(x,3,NP)';
    R1=2.85;
    R2=0.15;
    A=1.8308e+3;
    B=4.7118e+2;
    lemda1=2.4799;
    lemda2=1.7322;
    lemda3=1.7322;
    c=1.0039e+05;
    d=1.6218e+01;
    n1=7.8734e-01;
    gama=1.0999e-06;
    h=-5.9826e-01;
    E=zeros(1,NP);
    r=zeros(NP);
    for i=1:NP
        for j=1:NP
            r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
            if r(i,j)<(R1-R2)
                fcr(i,j)=1;
            elseif  r(i,j)>(R1+R2)
                fcr(i,j)=0;
            else
                fcr(i,j)=0.5-0.5*sin(pi/2*(r(i,j)-R1)/R2);
            end
            VRr(i,j)=A*exp(-lemda1*r(i,j));
            VAr(i,j)=B*exp(-lemda2*r(i,j));
        end
    end
    for i=1:NP
        for j=1:NP
            if i==j
                continue
            end
            jeta=zeros(NP,NP);
            for k=1:NP
                if i==k || j==k        
                    continue  
                end
                rd1=sqrt(sum((x(i,:)-x(k,:)).^2));
                rd3=sqrt(sum((x(k,:)-x(j,:)).^2));
                rd2=sqrt(sum((x(i,:)-x(j,:)).^2));
                ctheta_ijk=(rd1^2+rd2^2-rd3^3)/(2*rd1*rd2);
                G_th_ijk =1+(c^2)/(d^2)-(c^2)/(d^2+(h-ctheta_ijk)^2);
                jeta(i,j)=jeta(i,j)+fcr(i,k)*G_th_ijk*exp(lemda3^3*(r(i,j)-r(i,k))^3);
            end
            Bij=(1+(gama*jeta(i,j))^n1)^(-0.5/n1);
            E(i)=E(i)+fcr(i,j)*(VRr(i,j)-Bij*VAr(i,j))/2;
        end    
    end
    o=sum(E);
end

function o = F30(x)
    
    hsum=[];
    d=length(x);
    var=2*d-1;
    for kk=1:2*var
        if rem(kk,2)~=0
            i=(kk+1)/2;
            hsum(kk)=0;
            for j=i:d    
                summ=0;
                for i1=(abs(2*i-j-1)+1):j
                    summ=x(i1)+summ;
                end
                hsum(kk)= cos(summ)+ hsum(kk);
            end
        else 
            i=kk/2;
            hsum(kk)=0;
            for j=(i+1):d    
                summ=0;
                for i1=(abs(2*i-j)+1):j
                    summ=x(i1)+summ;
                end
                hsum(kk)= cos(summ)+ hsum(kk);
            end
            hsum(kk)=hsum(kk)+0.5;
        end
    end
    o=max(hsum); 
end

function o = F31(x)

    sw=ceil(x);
    data6Bus;
    n1=length(Linedata(:,1));
    sw1=sw;
    for k=1:length(sw1)
        if isnan(sw1(k))
            disp('error!');
        end    
        Linedata(n1+k,:)=Candidate(sw1(k),:);
    end
    n_orginalLine=n1;
    n=length(Pgen);
    B=zeros(n,n);
    Nline=length(Linedata(:,1));

    Xline=Linedata(:,4);
    pijmax=Linedata(:,6);
    Tap=ones(n);
    for C=1:Nline
        bline(C)=1/Xline(C);
        k=Linedata(C,2);
        m=Linedata(C,3);
        B(k,m)=B(k,m)-(bline(C));
        B(m,k)=B(k,m);
        B(k,k)=B(k,k)+(bline(C));
        B(m,m)=B(m,m)+(bline(C));
    end    
    B(1,1)=10000000;
    X=inv(B);
    delP= Pgen-Pload;
    delP=(delP');

    delta=X*(delP);
    pij=zeros(Nline,1);
    for k=1:Nline
        i=Linedata(k,2);
        j=Linedata(k,3);
        pij(k)=(delta(i)-delta(j))/Xline(k);
    end
    PIPbase=0.0;
    o=sum(Linedata(n_orginalLine+1:end,7))+30;
    pen=0;
    for i=1:length(Linedata(:,1))
        pen=pen+5000*max((abs(pij(i))-Linedata(i,6)),0);
    end
    
    for i=1:length(Candidate(:,1))
        [a ]=find(sw==i);
        if length(a)>3
            pen=pen+1000;
        end
    end
    o=o+pen;
end

function o = F32(x)
    Kp=100;
    EBEinputfile;
    n=length(bus_spec(:,1));
    Pg=(bus_spec(:,7))/100;
    Pd=(bus_spec(:,5))/100;
    na=linedata(:,1); nb=linedata(:,2);
    g=find(Pg>0);
    d=find(Pd>0); 
    Pg = Pg(g); 
    Pd = Pd(d);
    BT=zeros(length(g),length(d));

    BT(1,4)=5;BT(1,5)=10;BT(1,6)=5;
    BT(2,3)=5;
    BT(3,21)=2.5;
    BT(4,21)=2.5;BT(4,16)=15;
    BT(5,12)=2.5;BT(6,8)=2.5;
    
    BT=BT/100;

    Pg2 = sum(BT,2);  %  generations involved in BT 
    Pd2 = sum(BT,1);  %  loads involved in BT

    GD = zeros(length(g),length(d));
    for i=1:length(g)
        GD(i,:)=x(((i-1)*length(d)+1):(i*length(d)));
    end
    %% calculation of PTDF
    line_data=linedata;
    [YIbus] = EBEformybus(line_data,n);
    YIbus(1,:)=[];
    YIbus(:,1)=[];
    Xt=inv(YIbus);
    X=zeros(n,n);
    for i=0:(n-1)
        for j=0:(n-1)
        if (i~=0) && (j~=0)
            X(i+1,j+1)=Xt(i,j);
        else
            X(i+1,j+1)=0;
        end
        end
    end
    % PTDF 
    xij=imag(line_data(:,3));
    for i=1:length(g)
        for j=1:length(d)
            for k=1:length(na)
            PTDF(i,j,k)= [X(na(k),g(i))-X(nb(k),g(i))-X(na(k),d(j))+X(nb(k),d(j))]/xij(k);
            end
        end
    end
    %% calculation of objective fn.
    Rg = [32.7290   32.1122   30.3532   33.6474   64.1156   66.8729];
    Rd = [7.5449 10.7964 10.9944 11.0402 11.7990 15.3803 42.6800 41.4551 ...
          73.1939 57.0430 45.5920 43.6553 61.8002 59.6409 57.0279 ...
          51.0749 67.1070 60.6623 198.6744  178.9956  199.9483];
    FC=100*xij/sum(xij);
    flows=zeros(length(na));
    cost_line=zeros(length(na));
    for k=1:length(na)
        for i=1:length(g)
            for j=1:length(d)
                flows(k)=flows(k) + abs(PTDF(i,j,k)*GD(i,j)) + abs(PTDF(i,j,k)*BT(i,j));
            end
        end
        cost_line(k)=FC(k)/flows(k);
    end
    pg = Pg - Pg2;
    pd = Pd - Pd2';
    rate_ebe1=0;
    cost_l=zeros(length(g),length(d));
    cost_gen=zeros(length(g));
    for i=1:length(g)
        for j=1:length(d)
            for k=1:length(na)
                cost_l(i,j) = cost_l(i,j) + abs(cost_line(k)*PTDF(i,j,k));
            end
            cost_gen(i) = cost_gen(i) + GD(i,j)*cost_l(i,j);
        end
        rate_ebe1=rate_ebe1 + (cost_gen(i)/pg(i) - Rg(i))^2;
    end
    rate_ebe2=0;
    cost_load=zeros(length(d));
    for j=1:length(d)
        for i=1:length(g)   
            cost_load(j) = cost_load(j) + GD(i,j)*cost_l(i,j);
        end
        rate_ebe2=rate_ebe2 + (cost_load(j)/pd(j) - Rd(j))^2;
    end    
    rate_d = rate_ebe1 + rate_ebe2;
    %%    CONSTRAINT  VIOLATIONS  
    Pg_x = sum(GD,2)+ sum(BT,2);
    Pd_x = sum(GD,1)+ sum(BT,1);
    Gpen = 0;
    for i=1:length(g)
        Gpen = Gpen + 100*abs(Pg_x(i)-Pg(i));
    end
    LDpen = 0;
    for i=1:length(d)
       LDpen = LDpen + 100*abs(Pd_x(i)-Pd(i));
    end
    PENALTY  = Gpen + LDpen;
    o = rate_d + 50*Kp*PENALTY;
end

function o = F33(x)   
    null=[50,120];
    phi_desired=180;
    distance=0.5;
    o=antennafunccircular(x,null,phi_desired,distance);
end

function o = F34(x)
    x=round(x*10000)/10000; %% For fixing the 4 digit precision
    if length(x)==120
        o = fn_DED_5(x);        %%  fn_DED_5         F34
    else
        o = fn_DED_10(x);        %% fn_DED_10        F35
    end
end

function o = F36(x)
    x=round(x*10000)/10000; %% For fixing the 4 digit precision
    if length(x)==6
        o = fn_ELD_6(x);        % fn_ELD_6      F35
    elseif length(x)==13
        o = fn_ELD_13(x);       % fn_ELD_13     F37
    elseif length(x)==15
        o = fn_ELD_15(x);       % fn_ELD_15     F38
    elseif length(x)==40
        o = fn_ELD_40(x);       % fn_ELD_40     F39
    else
        o = fn_ELD_140(x);      % fn_ELD_140    F40
    end    
end

function o = F41(x)        
   o = fn_HT_ELD_Case_1(x); 
end

function o = F42(x)         
   o = fn_HT_ELD_Case_2(x); 
end

function o = F43(x)        
   o = fn_HT_ELD_Case_3(x); 
end

function o = F44(x)        
    load(".\CEC\CEC2011\12-13\messengerfull.mat");
    problem=MGADSMproblem;
    o=messengerfull(x,problem);
end

function o = F45(x)        
    load(".\CEC\CEC2011\12-13\cassini2.mat");
    problem=MGADSMproblem;
    o=cassini2(x,problem);
end

function o = F46(x)
    fun=1;
    o= ZHX2019(x', fun);
end

function o = F47(x)
    fun=2;    
    o= ZHX2019(x', fun);
end

function o = F48(x)
    fun=3;  
    o= ZHX2019(x', fun);
end

function o = F49(x)
    fun=4;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F50(x)
    fun=5;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F51(x)
    fun=6;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F52(x)
    fun=7;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F53(x)
    fun=8;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F54(x)
    fun=9;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F55(x)
    fun=10;
    if length(x)<10
        x=[x,(2*rand(1,10-length(x))-1).*100];
    end
    o= ZHX2019(x', fun);
end

function o = F56(x)
    fun=1;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F57(x)
    fun=2;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F58(x)
    o=0;
    d = length(x);
    temp1=0;
    temp2=0;
    for i=1:d-1
        temp1=sin(sqrt(x(i)^2+x(i+1)^2))^2;
        temp2=1+0.001*(x(i)^2+x(i+1)^2);
        o=o+0.5+(temp1-0.5)/temp2^2;
    end
    temp1=sin(sqrt(x(1)^2+x(d)^2))^2;
    temp2=1+0.001*(x(1)^2+x(d)^2);
    o=o+0.5+(temp1-0.5)/(temp2^2);
end

function o = F59(x)
    o=0;
    d = length(x);
    for i=1:d
       o=o+x(i)^2-10*cos(2*pi*x(i))+10;
    end
end

function o = F60(x)
    fun=5;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F61(x)
    fun=6;
    if length(x)<10
        
        x=[x,(2*rand(1,20-length(x))-1).*100];
    end    
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F62(x)
    fun=7;
    if length(x)<10
        x=[x,(2*rand(1,20-length(x))-1).*100];
    end
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F63(x)
    fun=8;
    if length(x)<10
        x=[x,(2*rand(1,20-length(x))-1).*100];
    end
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F64(x)
    fun=9;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F65(x)
    fun=10;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end
  
function o = F66(x)
    fun=11;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o = F67(x)
    fun=12;
    fhd=@ZHX2022;
    o= feval(fhd, x', fun);
end

function o=F68(x)                                
x1= x(1); x2 = x(2); x3 = x(3); x4 = x(4);
f = 0.6224*x1*x3*x4 + 1.7781*x2*x3^2+3.1661*x1^2*x4+19.84*x1^2*x3;
panaty_factor = 10e100; 

g1 = -x1+0.0193*x3;
panalty_1 = panaty_factor*(max(0,g1))^2;
g2 = -x2+0.00954*x3;
panalty_2 = panaty_factor*(max(0,g2))^2;
g3 = -pi*x3^2*x4 - (4/3)*pi*x3^3 + 1296000;
panalty_3 = panaty_factor*(max(0,g3))^2;
g4 = x4 - 240;
panalty_4 = panaty_factor*(max(0,g4))^2;
o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4;
end

function o=F69(x)                              
    x1= x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    x6 = x(6);
    x7 = x(7);
    
    f = 0.7854*x1*x2^2*(3.3333*x3^2 + 14.9334*x3 - 43.0934) -...
        1.508*x1*(x6^2 + x7^2) + 7.4777*(x6^3+x7^3) + 0.7854*(x4*x6^2 + x5*x7^2);
    
    panaty_factor = 10e100; 
    
    g1 = (27/(x1*x2^2*x3))-1 ;
    penalty_g1 = panaty_factor*(max(0,g1))^2;
    %  g2
    g2 = (397.5/(x1*x2^2*x3^2)) - 1;
    penalty_g2 = panaty_factor*(max(0,g2))^2;
    %  g3
    g3 = (1.93*x4^3/(x2*x3*x6^4)) - 1;
    penalty_g3 = panaty_factor*(max(0,g3))^2;
    %  g4
    g4 = (1.93*x5^3/(x2*x3*x7^4)) - 1;
    penalty_g4 = panaty_factor*(max(0,g4))^2;
    %  g5
    g5 = (1/(110*x6^3))*(((745*x4/(x2*x3))^2 + 16.9*1e6)^0.5)-1;
    penalty_g5 = panaty_factor*(max(0,g5))^2;
    %  g6
    g6 = (1/(85*x7^3))*(((745*x5/(x2*x3))^2 + 157.5*1e6)^0.5)-1;
    penalty_g6 = panaty_factor*(max(0,g6))^2;
    %  g7
    g7 = (x2*x3/40) - 1;
    penalty_g7 = panaty_factor*(max(0,g7))^2;
    %  g8
    g8 = (5*x2/x1) - 1;
    penalty_g8 = panaty_factor*(max(0,g8))^2;
    %  g9
    g9 = (x1/(12*x2)) - 1;
    penalty_g9 = panaty_factor*(max(0,g9))^2;
    % g10
    g10 = ((1.5*x6+1.9)/(x4)) - 1;
    penalty_g10 = panaty_factor*(max(0,g10))^2;
    % g11
    g11 = ((1.1*x7+1.9)/(x5)) - 1;
    penalty_g11 = panaty_factor*(max(0,g11))^2;
    %
     o = f+penalty_g1 + penalty_g2 + penalty_g3 + penalty_g4 + penalty_g5 + penalty_g6 + penalty_g7 +...
        penalty_g8 + penalty_g9 + penalty_g10 + penalty_g11;
end

function o = F70(x)            
     P = 6000;
    L=14;
    E=30e6;
    G = 12e6;
    tmax=13600;
    sigma_max=30000;
    deta_max=0.25;
    x1= x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    % M
    M=P*(L+0.5*x2);
    % R
    r1 = (x2^2)/4;
    r2 = ((x1+x3)/2)^2;
    R = (r1+r2)^0.5;
    % J
    j1=sqrt(2)*x1*x2;
    j2 = (x2^2)/12;
    j3 = ((x1+x3)/2)^2;
    J=2*(j1*(j2+j3));
    %
    sigma_x=6*P*L/(x4*x3^2);
    deta_x = 4*P*L^3/((E*x4)*(x3^3));
    % Pc
    p1 = (4.013*E*((x3^2)*(x4^6)/36)^0.5)/(L^2);
    p2 = (x3/(2*L))*(E/(4*G))^0.5;
    Pc = p1*(1-p2);
    %
    t_1 =  P/(sqrt(2*x1*x2));
    t_2 = M*R/J;
    t=((t_1)^2 + 2*t_1*t_2*(x2/(2*R))+(t_2)^2)^0.5;
    
    f = 1.10471*(x1^2)*x2+0.04811*x3*x4*(14+x2);
    panaty_factor = 10e100; 
    % g1
    g1 = t-tmax;
    penalty_g1 = panaty_factor*(max(0,g1))^2;
    % g2
    g2 = sigma_x-sigma_max;
    penalty_g2 = panaty_factor*(max(0,g2))^2;
    % g3
    g3 = deta_x - deta_max;
    penalty_g3 = panaty_factor*(max(0,g3))^2;
    % g4
    g4 = x1-x4;
    penalty_g4 = panaty_factor*(max(0,g4))^2;
    % g5
    g5 = P - Pc;
    penalty_g5 = panaty_factor*(max(0,g5))^2;
    % g6
    g6 = 0.125-x1;
    penalty_g6 = panaty_factor*(max(0,g6))^2;
    % g7
    g7 =1.10471*(x1^2)*x2+0.04811*x3*x4*(14+x2)-5;
    penalty_g7 = panaty_factor*(max(0,g7))^2;
    %
    o = f+penalty_g1+penalty_g2+penalty_g3+penalty_g4+penalty_g5 +penalty_g6+penalty_g7;
end

function o = F71(x)              
    panaty_factor = 10e100;
    D=160;d=90;Bw=30;
    Dm=x(1); Db=x(2); Z=x(3); fi=x(4); f0=x(5); KDmin=x(6); KDmax=x(7);
    ep = x(8); ee=x(9); xi=x(10);
    Z=round(Z);
    % ri=11.033; r0=11.033; fi=ri/Db; f0=r0/Db;
    T=D-d-2*Db;
    phio=2*pi-acos(((((D-d)/2)-3*(T/4))^2+(D/2-T/4-Db)^2-(d/2+T/4)^2)...
        /(2*((D-d)/2-3*(T/4))*(D/2-T/4-Db)));
    %
    g(1)=1+phio/(2*asin(Db/Dm))-Z;
    g(2)=-2*Db+KDmin*(D-d);
    g(3)= -KDmax*(D-d)+2*Db;
    g(4)=xi*Bw-Db;
    g(5)=-Dm+0.5*(D+d);
    g(6)=-(0.5+ee)*(D+d)+Dm;
    g(7)=-0.5*(D-Dm-Db)+ep*Db;
    g(8)=0.515-fi;
    g(9)=0.515-f0;
    
    penalty=panaty_factor*sum(g(g>0).^2);
    
    gama=Db/Dm;
    fc=37.91*((1+(1.04*((1-gama/1+gama)^1.72)*((fi*(2*f0-1)/f0*...
        (2*fi-1))^0.41))^(10/3))^-0.3)*((gama^0.3*(1-gama)^1.39)/...
        (1+gama)^(1/3))*(2*fi/(2*fi-1))^0.41;
    if Db<=25.4
        f=-fc*Z^(2/3)*Db^1.8;
    else
        f=-3.647*fc*Z^(2/3)*Db^1.4;
    end
    o = f+penalty;
end

function o=F72(x)     
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    f = (x3+2)*x2*(x1^2);
    % 
    g1 = 1-((x2^3)*x3)/(71785*(x1^4));
    g2 = (4*(x2^2)-x1*x2)/(12566*(x2*(x1^3)-(x1^4))) + 1/(5108*(x1^2))-1;
    g3 = 1-(140.45*x1)/((x2^2)*x3);
    g4 = ((x1+x2)/1.5)-1;
    panaty_1 = 10e100*(max(0,g1))^2; 
    panaty_2 = 10e100*(max(0,g2))^2; 
    panaty_3 = 10e100*(max(0,g3))^2; 
    panaty_4 = 10e100*(max(0,g4))^2; 
    o = f + panaty_1+panaty_2+panaty_3+panaty_4;
end
 
function o=F73(x)                            
    x = round(x);
    o = (1/6.931-(x(2)*x(3)/(x(1)*x(4))))^2;
end

function o=F74(x)
    panalty_factor = 10e100; 
    d1 = x(1)*1e-3; d2 = x(2)*1e-3; d3 = x(3)*1e-3; d4 = x(4)*1e-3; w = x(5)*1e-3;
    N = 350; N1 = 750; N2 = 450; N3 = 250; N4 = 150;
    rho = 7200; a = 3; mu = 0.35; s = 1.75*1e6; t = 8*1e-3;
   
    C1 = pi*d1/2*(1+N1/N)+(N1/N-1)^2*d1^2/(4*a)+2*a;
    C2 = pi*d2/2*(1+N2/N)+(N2/N-1)^2*d2^2/(4*a)+2*a;
    C3 = pi*d3/2*(1+N3/N)+(N3/N-1)^2*d3^2/(4*a)+2*a;
    C4 = pi*d4/2*(1+N4/N)+(N4/N-1)^2*d4^2/(4*a)+2*a;
    R1 = exp(mu*(pi-2*asin((N1/N-1)*d1/(2*a))));
    R2 = exp(mu*(pi-2*asin((N2/N-1)*d2/(2*a))));
    R3 = exp(mu*(pi-2*asin((N3/N-1)*d3/(2*a))));
    R4 = exp(mu*(pi-2*asin((N4/N-1)*d4/(2*a))));
    P1 = s*t*w*(1-exp(-mu*(pi-2*asin((N1/N-1)*d1/(2*a)))))*pi*d1*N1/60;
    P2 = s*t*w*(1-exp(-mu*(pi-2*asin((N2/N-1)*d2/(2*a)))))*pi*d2*N2/60;
    P3 = s*t*w*(1-exp(-mu*(pi-2*asin((N3/N-1)*d3/(2*a)))))*pi*d3*N3/60;
    P4 = s*t*w*(1-exp(-mu*(pi-2*asin((N4/N-1)*d4/(2*a)))))*pi*d4*N4/60;
    
    g(1) = -R1+2;
    g(2) = -R2+2;
    g(3) = -R3+2;
    g(4) = -R4+2;
    g(5) = -P1+(0.75*745.6998);
    g(6) = -P2+(0.75*745.6998);
    g(7) = -P3+(0.75*745.6998);
    g(8) = -P4+(0.75*745.6998);
    h(1) = C1-C2;
    h(2) = C1-C3;
    h(3) = C1-C4;
    
    panalty=panalty_factor*(sum(g(g>0).^2) + sum(h(h~=0).^2));
    
    f = rho*w*pi/4*(d1^2*(1+(N1/N)^2)+d2^2*(1+(N2/N)^2)+d3^2*(1+(N3/N)^2)+d4^2*(1+(N4/N)^2));
    o=f+panalty;
end

function o=F75(x)                                     
l = 100; P = 2; q = 2;
x1= x(1);
x2 = x(2);
f = l*(2*sqrt(2)*x1+x2);
panaty_factor = 10e100; 
%
g1 = P*(sqrt(2)*x1+x2)/(sqrt(2)*x1^2+2*x1*x2)-q;
penalty_g1 = panaty_factor*(max(0,g1))^2;
g2 = P*(x2)/(sqrt(2)*x1^2+2*x1*x2)-q;
penalty_g2 = panaty_factor*(max(0,g2))^2;
g3 = P/(sqrt(2)*x2+x1)-q;
penalty_g3 = panaty_factor*(max(0,g3))^2;
o = f+penalty_g1+penalty_g2+penalty_g3;
end

function o=F76(x)               
    panalty_factor = 10e100; 
    g(1)=61/x(1)^3+37/x(2)^3+19/x(3)^3+7/x(4)^3+1/x(5)^3-1;
    penalty=panalty_factor*sum(g(g>0).^2);
    o=0.0624*sum(x)+penalty;
end

function o=F77(x)                         
    panalty_factor = 10e100; 
    x = round(x); 
    
    Mf = 3; Ms = 40; Iz = 55; n = 250; Tmax = 15; s = 1.5; delta = 0.5;
    Vsrmax = 10; rho = 0.0000078; pmax = 1; mu = 0.5; Lmax = 30; delR = 20;
    Mh = 2/3*mu*x(4)*x(5)*(x(2)^3-x(1)^3)/(x(2)^2-x(1)^2);
    Prz = x(4)/(pi*(x(2)^2-x(1)^2));
    Vsr = (2*pi*n/90)*(x(2)^3-x(1)^3)/(x(2)^2-x(1)^2);
    T   = (Iz*pi*n/30)/(Mh+Mf);
    
    g(1) = -x(2)+x(1)+delR;
    g(2) = (x(5)+1)*(x(3)+delta)-Lmax;
    g(3) = Prz-pmax;
    g(4) = Prz*Vsr-pmax*Vsrmax;
    g(5) = Vsr-Vsrmax;
    g(6) = T-Tmax;
    g(7) = s*Ms-Mh;
    g(8) = -T;
    penalty=panalty_factor*sum(g(g>0).^2);
    
    f = pi*(x(2)^2-x(1)^2)*x(3)*(x(5)+1)*rho;
    o=f+penalty;
end