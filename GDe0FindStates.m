function Coefs = GDe0FindStates(n,dotsperN0)
%-------------
if n==0
    Coefs1=cat(2,transpose(linspace(0,pi,dotsperN0)),transpose(linspace(0,2*pi,dotsperN0)),zeros(dotsperN0,1));
    Coefs2=cat(2,zeros(dotsperN0,1),transpose(linspace(0,2*pi,dotsperN0)),transpose(linspace(0,1,dotsperN0)));
    Coefs=datasample(cat(1,Coefs1,Coefs2),dotsperN0);
    return;
elseif n==1
    Coefs=cat(2,pi*ones(dotsperN0,1),transpose(linspace(0,2*pi,dotsperN0)),ones(dotsperN0,1));
    return;
else
    birds=ceil(10*sqrt(dotsperN0));
    TestCoefs=GetStates(0,pi,birds,n);
    if size(TestCoefs,1)>1
        Coefs=datasample(GetStates(min(TestCoefs(:,1)),max(TestCoefs(:,1)),2*birds,n),dotsperN0);
    else
        Coefs=zeros(3,1);
    end
end
end
function n=GetN0(Theta,Phi,y)
   n=2.*(y.*sqrt(3+cos(Theta))).*(sin(Theta./2).^2)./(sqrt(2).*(1+y.^2+2.*y.*cos(Theta./2).^3.*cos(Phi)));%for Gamma
end
function y=GetY(Theta,Phi,n)
y1=(-2*n*cos(Theta/2)^3*cos(Phi)+sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2-sqrt(-4*n^2+(2*n*cos(Theta/2)^3*cos(Phi)-sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2)^2))/(2*n);
y2=(-2*n*cos(Theta/2)^3*cos(Phi)+sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2+sqrt(-4*n^2+(2*n*cos(Theta/2)^3*cos(Phi)-sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2)^2))/(2*n);
if isreal(y1)&&0<=y1&&y1<1&&abs(GetN0(Theta,Phi,y1)-n)<=10^(-16)
    y=y1;
    return;
elseif isreal(y2)&&0<=y2&&y2<1&&abs(GetN0(Theta,Phi,y2)-n)<=10^(-16)
    y=y2;
    return;
else
    y=-1;%fail to find y
end
end
function  Coefs=GetStates(minTheta,maxTheta,birds,n)
Theta=linspace(minTheta,maxTheta,birds);
Phi=linspace(0,2*pi-10^(-60),birds);
PossibleCoefs=zeros(birds*birds,3);
CoefsCount=0;
for i=1:birds
    for j=1:birds
        Possibley=GetY(Theta(i),Phi(j),n);
        if Possibley~=-1
            CoefsCount=CoefsCount+1;
            PossibleCoefs(CoefsCount,:)=[Theta(i),Phi(j),Possibley];
        end
    end
end
if CoefsCount>2
    Coefs=PossibleCoefs(1:CoefsCount,:);
else
    Coefs=zeros(1,3);
end
end