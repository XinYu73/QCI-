function [bestAvailability,optimised_parameters,BIRDS,Availability] = GDePSO(Bird_in_swarm, Number_of_quality_in_Bird, MinMaxRange,availability_type,max_iteration,n0,measurement,velocity_clamping_factor, cognitive_constant, social_constant, Min_Inertia_weight, Max_Inertia_weight)
% Bird_in_swarm=Number of particle=agents=candidate
% Number_of_quality_in_Bird=Number of Variable
% MinMaxRange: jx2 matrix; jth row contains minimum and maximum values of the jth variable
% say you have a variable N1 ,which can have maximum value M1 and minimum value m1,then your matrix will be [m1 M1] ,for more:[m1 M1; m2 M2; mj Mj]
% availability_type is string 'min' or 'max' to check depending upon need to minimize or maximize the Food_availability
% Inertia_weight=At the beginning of the search procedure, diversification is heavily weighted, while intensification is heavily weighted at the end of the search procedure.
%-----------Checking all parameteres are entered
format long
if nargin < 10
    velocity_clamping_factor=2;% velocity_clamping_factor (normally 2)
    cognitive_constant=2;
    social_constant=2;% cognitive_constant=c1=individual learning rate (normally 2),social_constant=c2=social parameter (normally 2),normally C1+C2>=4
    Min_Inertia_weight=0.4;% Min_Inertia_weight=min of inertia weight (normally 0.4)
    Max_Inertia_weight=0.9;% Max_Inertia_weight=max of inertia weight (normally 0.9)
end
disp(n0);
availability_type=lower(availability_type(1:3));%universalize availability type 
[row,col]=size(MinMaxRange);% Checking for proper boundary Values and entered Matrix
if row~=Number_of_quality_in_Bird || col~=2
    error('Not a proper MinMaxRange Matrix')
end
for i=1:Number_of_quality_in_Bird
    if MinMaxRange(i,1)>=MinMaxRange(i,2)
        error('Minimum value greater than Maximum value!!!')
    end
end
bird_min_range=MinMaxRange(:,1);%distinguishing min and max range
bird_max_range=MinMaxRange(:,2);
format long;
%------------------------initialize birds
dotsmultiple=2;
tryhits=5;
coefs=GDe0FindStates(n0,tryhits);
while length(coefs)<=Bird_in_swarm
        tryhits=dotsmultiple*tryhits;
        coefs=GDe0FindStates(n0,tryhits);
end
coefs=datasample(GDe0FindStates(n0,tryhits),Bird_in_swarm);
%----------------------delta<=10^{-16}->hard
%to find enough legal position so we us RS
%methods here.
if n0==1||n0==0
   switch availability_type
       case 'min'
            %-----------------------------down
            downp=GDe2P(coefs(1,1),coefs(1,2),coefs(1,3),measurement);%更新p
            downcoef=coefs(1,:);
            for j=2:length(coefs)
                if GDe1Entanglement(downp,coefs(j,1),coefs(j,2),coefs(j,3),measurement)<0%当前p太大
                    downp=GDe2P(coefs(j,1),coefs(j,2),coefs(j,3),measurement);%更新p
                    downcoef=coefs(j,:);
                end
            end
            bestAvailability=downp;
            optimised_parameters=downcoef;
            BIRDS=0;
            Availability=0;
            return;
       case 'max'
            upp=GDe2P(coefs(1,1),coefs(1,2),coefs(1,3),measurement);%更新p
            upcoef=coefs(1,:);
            %-----------------------------up
            for j=2:length(coefs)
                if GDe1Entanglement(upp,coefs(j,1),coefs(j,2),coefs(j,3),measurement)>0%当前p太小
                    upp=GDe2P(coefs(j,1),coefs(j,2),coefs(j,3),measurement);%更新p
                    upcoef=coefs(j,:);
                end
            end
            bestAvailability=upp;
            optimised_parameters=upcoef;
            BIRDS=0;
            Availability=0;
            return;
   end
end
%----------------
Bird_in_swarm=length(coefs);
bird=zeros(Bird_in_swarm,Number_of_quality_in_Bird,max_iteration);
bird(:,:,1)=coefs(:,1:2);
%-------------------initialize Velocity for each bird
Vmax=bird_max_range*velocity_clamping_factor;
Vmin=-Vmax;
Velocity=zeros(Bird_in_swarm,Number_of_quality_in_Bird,max_iteration);
for i=1:Number_of_quality_in_Bird
    Velocity(:,i,1)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(Bird_in_swarm,1);
end
%-----------------initialize,pbest,lbest,gbest 
pBest=zeros(Bird_in_swarm,Number_of_quality_in_Bird);%position
lBest=zeros(Bird_in_swarm,Number_of_quality_in_Bird);
gBest=zeros(max_iteration,Number_of_quality_in_Bird);
gBest_availability=zeros(max_iteration,1);
availability=zeros(Bird_in_swarm,max_iteration);%to store personal best
 for p=1:Bird_in_swarm
     availability(p,1)=Food_availability(bird(p,:,1),n0,availability_type,measurement);
     pBest(p,:)=bird(p,:,1);
 end
 switch availability_type
     case 'min'
                    [~,index]=min(availability(:,1));
                    gBest_availability(1)=availability(index,1);
                    gBest(1,:)=bird(index,:,1);
     case 'max'
                    [~,index]=max(availability(:,1));
                    gBest_availability(1)=availability(index,1);
                    gBest(1,:)=bird(index,:,1);
     otherwise
         error('availability_type mismatch')
 end
%------------------------r3PSO
for p=1:Bird_in_swarm
    if p==1
       r3positon=[pBest(1,:);pBest(2,:);pBest(Bird_in_swarm,:)];% 3x2
       r3availability=[availability(1,1),availability(2,1),availability(Bird_in_swarm,1)];
    elseif p==Bird_in_swarm
       r3positon=[pBest(Bird_in_swarm-1,:);pBest(Bird_in_swarm,:);pBest(1,:)];% 3x2
       r3availability=[availability(Bird_in_swarm-1,1),availability(Bird_in_swarm,1),availability(1,1)];
    else
       r3positon=[pBest(p-1,:);pBest(p,:);pBest(p+1,:)];% 3x2
       r3availability=[availability(p-1,1),availability(p,1),availability(p+1,1)];
    end
    switch availability_type
        case 'min'
            [~,lindex]=min(r3availability);
            lBest(p,:)=r3positon(lindex,:);
        case 'max'
            [~,lindex]=max(r3availability);
            lBest(p,:)=r3positon(lindex,:);
    end
end    
 %----------------------first flight
 w=zeros(max_iteration,1);
  for p=1:Bird_in_swarm
        w(1)=((max_iteration - 1)*(Max_Inertia_weight - Min_Inertia_weight))/(max_iteration-1) + Min_Inertia_weight;%Linearly decreasinginertia weight
        %-----------------------PSO
%       Velocity(p,:,(1+1))=w(1)*Velocity(p,:,1) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest(1,:)-bird(p,:,1)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:,1));
        %-----------------------r3PSO
        Velocity(p,:,2)=w(1)*Velocity(p,:,1) + social_constant*rand(1,Number_of_quality_in_Bird).*(lBest(p,:)-bird(p,:,1)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:,1)); 
        Velocity(p,:,2)=MinMaxCheck(Vmin, Vmax, Velocity(p,:,2));
        bird(p,:,2)= bird(p,:,1) + Velocity(p,:,2);
        bird(p,:,2)=MinMaxCheck(bird_min_range, bird_max_range, bird(p,:,2));
  end
%************************Start Searching
for itr=2:max_iteration
     fprintf('Completed  %d  %% ...\n', uint8(itr*100/max_iteration));
    %----------------calculate availability
    %and update personal best
    for p=1:Bird_in_swarm
        switch availability_type
            case 'min'
                availability(p,itr)=Food_availability(bird(p,:,itr),n0,availability_type,measurement,min(availability(p,1:(itr-1))));
                if availability(p,itr)==-1%if we get a outlaw position
                    dotsmultiple=2;
                    tryhits=5;
                    newPossiblePosition=GDe0FindStates(n0,tryhits);
                    while length(newPossiblePosition)<=10
                            tryhits=dotsmultiple*tryhits;
                            newPossiblePosition=GDe0FindStates(n0,tryhits);
                    end
                    newPossiblePosition=datasample(GDe0FindStates(n0,tryhits),1);
                    bird(p,:,itr)=newPossiblePosition(1:2);
                    availability(p,itr)=Food_availability(bird(p,:,itr),n0,availability_type,measurement,min(availability(p,1:(itr-1))));%calculate the new availiability
                end
                [~,index]=min(availability(p,1:itr));%min(availability(p,:));,1:itr is meaningful value.
                pBest(p,:)=bird(p,:,index);%personal best 
            case 'max'
                availability(p,itr)=Food_availability(bird(p,:,itr),n0,availability_type,measurement,max(availability(p,1:(itr-1))));
                if availability(p,itr)==-1%if we get a outlaw position
                    dotsmultiple=2;
                    tryhits=5;
                    newPossiblePosition=GDe0FindStates(n0,tryhits);
                    while length(newPossiblePosition)<=10
                            tryhits=dotsmultiple*tryhits;
                            newPossiblePosition=GDe0FindStates(n0,tryhits);
                    end
                    newPossiblePosition=datasample(GDe0FindStates(n0,tryhits),1);
                    bird(p,:,itr)=newPossiblePosition(1:2);
                    availability(p,itr)=Food_availability(bird(p,:,itr),n0,availability_type,measurement,max(availability(p,1:(itr-1))));%calculate the new availiability
                end
                [~,index]=max(availability(p,1:itr));
                pBest(p,:)=bird(p,:,index);
        end
    end
    %--------------------update global best
    switch availability_type
            case 'min'
                 if min(availability(:,itr))<gBest_availability(itr-1)
                    [~,gindex]=min(availability(:,itr));
                    gBest_availability(itr)=availability(gindex,itr);
                    gBest(itr,:)=bird(gindex,:,itr);
                 else
                     gBest_availability(itr)=gBest_availability(itr-1);
                     gBest(itr,:)=gBest(itr-1,:);
                 end
            case 'max'
                 if max(availability(:,itr))>gBest_availability(itr-1)
                    [~,gindex]=max(availability(:,itr));
                    gBest_availability(itr)=availability(gindex,itr);
                    gBest(itr,:)=bird(gindex,:,itr);
                 else
                     gBest_availability(itr)=gBest_availability(itr-1);
                     gBest(itr,:)=gBest(itr-1,:);
                 end             
    end
   %--------------------updata local best
   for p=1:Bird_in_swarm
        switch availability_type
            case 'min'
                if p==1
                    r3positon=[pBest(1,:);pBest(2,:);pBest(Bird_in_swarm,:)];% 3x2
                    r3availability=[min(availability(1,1:itr)),min(availability(2,1:itr)),min(availability(Bird_in_swarm,1:itr))];
                elseif p==Bird_in_swarm
                    r3positon=[pBest(Bird_in_swarm-1,:);pBest(Bird_in_swarm,:);pBest(1,:)];% 3x2
                    r3availability=[min(availability(Bird_in_swarm-1,1:itr)),min(availability(Bird_in_swarm,1:itr)),min(availability(1,1:itr))];
                else
                    r3positon=[pBest(p-1,:);pBest(p,:);pBest(p+1,:)];% 3x2
                    r3availability=[min(availability(p-1,1:itr)),min(availability(p,1:itr)),min(availability(p+1,1:itr))];
                end
                [~,lindex]=min(r3availability);
                lBest(p,:)=r3positon(lindex,:);
            case 'max'
                if p==1
                    r3positon=[pBest(1,:);pBest(2,:);pBest(Bird_in_swarm,:)];% 3x2
                    r3availability=[max(availability(1,1:itr)),max(availability(2,1:itr)),max(availability(Bird_in_swarm,1:itr))];
                elseif p==Bird_in_swarm
                    r3positon=[pBest(Bird_in_swarm-1,:);pBest(Bird_in_swarm,:);pBest(1,:)];% 3x2
                    r3availability=[max(availability(Bird_in_swarm-1,1:itr)),max(availability(Bird_in_swarm,1:itr)),max(availability(1,1:itr))];
                else
                    r3positon=[pBest(p-1,:);pBest(p,:);pBest(p+1,:)];% 3x2
                    r3availability=[max(availability(p-1,1:itr)),max(availability(p,1:itr)),max(availability(p+1,1:itr))];
                end
                [~,lindex]=max(r3availability);
                lBest(p,:)=r3positon(lindex,:);
        end
   end  
    %--------------------update position
    for p=1:Bird_in_swarm
        w(itr)=((max_iteration - itr)*(Max_Inertia_weight - Min_Inertia_weight))/(max_iteration-1) + Min_Inertia_weight;%Linearly decreasinginertia weight
        %-----------------------r3PSO and PSO
%         if itr<=floor(2*max_iteration/3)
%             Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(lBest(p,:)-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:,itr)); 
%         else
%             Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest(itr,:)-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:,itr));
%         end
        %-----------------------PSO
       Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest(itr,:)-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:,itr));
        %-----------------------r3PSO 
%        Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(lBest(p,:)-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest(p,:)-bird(p,:,itr)); 
        Velocity(p,:,(itr+1))=MinMaxCheck(Vmin, Vmax, Velocity(p,:,(itr+1)));
        bird(p,:,(itr+1))= bird(p,:,itr) + Velocity(p,:,(itr+1));
        bird(p,:,(itr+1))=MinMaxCheck(bird_min_range, bird_max_range, bird(p,:,(itr+1)));
    end
end
%----------------------pic
figure
    plot(transpose(availability),'c-.',"MarkerSize",0.5);
    hold on
    plot(1:max_iteration,gBest_availability,'b-.','LineWidth',2);
    text(floor(2*max_iteration/3),0.05*mean(availability(:,1)),"$\leftarrow r3PSO | gbPSO \rightarrow  $","Interpreter","latex","FontSize",10,'HorizontalAlignment','center')
    ylabel(strcat("$","P_","{",measurement,"}","$"),"Interpreter","latex","FontSize",15);
    title(strcat("$","P \sim Iter-" ,"N0:",num2str(n0)," ","-Birds: ",num2str(Bird_in_swarm),"$"),"Interpreter","latex","FontSize",15);
    xlabel("$Iteration$","Interpreter","latex","FontSize",15);

%----------------find coefs and negaP
[bestAvailability,optimised_parameters]=Food_availability(gBest(end,:),n0,availability_type,measurement);
BIRDS=bird;
Availability=availability;
end
function [ A2comp ] = MinMaxCheck( minimum, maximum, A2comp )
if nargin < 3
    error('Missing input parameter(s)!')
end

if (length(minimum(:))==length(maximum(:)) && length(maximum(:))==length(A2comp(:)))
    
    size=length(minimum);
    for l=1:size
        if maximum(l)<minimum(l)
            error('Maximum value must be greater than minimum value !!')
        end
    end
    
    for l=1:size
        if(maximum(l)<A2comp(l)||minimum(l)>A2comp(l))
            if(maximum(l)<A2comp(l))
                A2comp(l)=maximum(l);
            else
                A2comp(l)=minimum(l);
            end
        end
    end
else
    error('All arrays must be same in length...')
end
end
function [Pneeded ,state] = Food_availability(ThetaPhi,n,availability_type,measurement,personbestAvailability)
format long
if nargin<5
    switch availability_type
        case 'min'
            personbestAvailability=1;
        case 'max'
            personbestAvailability=0;
    end
end
y=GetY(ThetaPhi(1),ThetaPhi(2),n);
if y==-1% test if current point is a legal position.
    Pneeded=personbestAvailability;
    state=zeros(3,1);
    return;
else
    state=[ThetaPhi(1),ThetaPhi(2),y];
    switch availability_type
        case 'min'
            if GDe1Entanglement(personbestAvailability,state(1),state(2),state(3),measurement)<0
                Pneeded=GDe2P(state(1),state(2),state(3),measurement);
            else
                Pneeded=personbestAvailability;%pretend everything is the same,this result should not influence its velocity
            end
        case 'max'
            if GDe1Entanglement(personbestAvailability,state(1),state(2),state(3),measurement)>0
                Pneeded=GDe2P(state(1),state(2),state(3),measurement);
            else
                Pneeded=personbestAvailability;
            end
    end
end
%----
end
function y=GetY(Theta,Phi,n)
y1=(-2*n*cos(Theta/2)^3*cos(Phi)+sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2-sqrt(-4*n^2+(2*n*cos(Theta/2)^3*cos(Phi)-sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2)^2))/(2*n);
y2=(-2*n*cos(Theta/2)^3*cos(Phi)+sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2+sqrt(-4*n^2+(2*n*cos(Theta/2)^3*cos(Phi)-sqrt(2)*sqrt(3+cos(Theta))*sin(Theta/2)^2)^2))/(2*n);
if isreal(y1)&&0<=y1&&y1<=1&&abs(GetN0(Theta,Phi,y1)-n)<=10^(-16)
    y=y1;
    return;
elseif isreal(y2)&&0<=y2&&y2<=1&&abs(GetN0(Theta,Phi,y2)-n)<=10^(-16)
    y=y2;
    return;
else
    y=-1;%fail to find y
end
end
function n=GetN0(Theta,Phi,y)
   n=2.*(y.*sqrt(3+cos(Theta))).*(sin(Theta./2).^2)./(sqrt(2).*(1+y.^2+2.*y.*cos(Theta./2).^3.*cos(Phi)));%for Gamma
end