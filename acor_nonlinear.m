clc;
clear;
close all;

%% Problem Definition
model=CreateModel();
CostFunction=@(x) MyCost(x,model);

x=model.x;
t=model.t;
n=model.n;

LBepsilon=1;
UBepsilon=100;
LBC=10;
UBC=1000;
LBsigma=1;
UBsigma=50;

nVar=3;              % "epsilon" , "C" and "sigma" Variables

VarSize=[1 nVar];    % Variables Matrix Size

VarMinepsilon=LBepsilon;         % Decision Variables Lower Bound
VarMaxepsilon=UBepsilon;         % Decision Variables Upper Bound
VarMinC=LBC;                     % Decision Variables Lower Bound
VarMaxC=UBC;                     % Decision Variables Upper Bound
VarMinsigma=LBsigma;             % Decision Variables Lower Bound
VarMaxsigma=UBsigma;             % Decision Variables Upper Bound

%% ACOR Parameters
MaxIt=100;          % Maximum Number of Iterations
nPop=50;            % Population Size (Archive Size)
nSample=30;         % Sample Size
q=0.1;              % Intensification Factor (Selection Pressure)
zeta=0.8;             % Deviation-Distance Ratio

%% Initialization

% Create Empty Individual Structure
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.sol=[];

% Create Population Matrix
pop=repmat(empty_individual,nPop,1);

% Initialize Population Members
for i=1:nPop
    
    % Create Random Solution
    pop(i).Position(1,1)=unifrnd(VarMinepsilon,VarMaxepsilon,1);
    pop(i).Position(1,2)=unifrnd(VarMinC,VarMaxC,1);
    pop(i).Position(1,3)=unifrnd(VarMinsigma,VarMaxsigma,1);
    %pop(i).Position(nVar)=unifrnd(VarMinN,VarMaxN,1)
    
    % Evaluation
    [pop(i).Cost,pop(i).sol]=CostFunction(pop(i).Position);
    
end

% Sort Population
[~, SortOrder]=sort([pop.Cost]);
pop=pop(SortOrder);

% Update Best Solution Ever Found
BestSol=pop(1);
WorstSol=pop(end);

BestCost=zeros(MaxIt,1);       % Array to Hold Best Cost Values
MeanCost=zeros(MaxIt,1);       % Array to Hold Mean Cost Values
WorstCost=zeros(MaxIt,1);      % Array to Hold Worst Cost Values
Bestepsilon=zeros(MaxIt,1);
BestC=zeros(MaxIt,1);
Bestsigma=zeros(MaxIt,1);

% Solution Weights
w=1/(sqrt(2*pi)*q*nPop)*exp(-0.5*(((1:nPop)-1)/(q*nPop)).^2);

% Selection Probabilities
p=w/sum(w);

%% ACOR Main Loop

for it=1:MaxIt
 tic   
    % Means
    s=zeros(nPop,nVar);
    for l=1:nPop
        s(l,:)=pop(l).Position;
    end
    
    % Standard Deviations
    Sigma=zeros(nPop,nVar);
    for l=1:nPop
        D=0;
        for r=1:nPop
            D=D+abs(s(l,:)-s(r,:));
        end
        Sigma(l,:)=zeta*D/(nPop-1);
    end
    
    % Create New Population Array
    newpop=repmat(empty_individual,nSample,1);
    for T=1:nSample
        
        % Initialize Position Matrix
        newpop(T).Position=zeros(VarSize);
        
        %%% Solution Construction
        l=RouletteWheelSelection(p);
        
        % Generate Gaussian Random Variable
        newpop(T).Position(1,1)=s(l,1)+Sigma(l,1)*randn;
        
        if newpop(T).Position(1,1)<VarMinepsilon
            
            newpop(T).Position(1,1)=VarMinepsilon;
            
        elseif newpop(T).Position(1,1)>VarMaxepsilon
            
            newpop(T).Position(1,1)=VarMaxepsilon;
        end
        
        % Select Gaussian Kernel
        l=RouletteWheelSelection(p);
        
        % Generate Gaussian Random Variable
        newpop(T).Position(1,2)=s(l,2)+Sigma(l,2)*randn;
        
        if newpop(T).Position(1,2)<VarMinC
            
            newpop(T).Position(1,2)=VarMinC;
            
        elseif newpop(T).Position(1,2)>VarMaxC
            
            newpop(T).Position(1,2)=VarMaxC;
        end
        
        % Generate Gaussian Random Variable
        newpop(T).Position(1,3)=s(l,3)+Sigma(l,3)*randn;
        
        if newpop(T).Position(1,3)<VarMinsigma
            
            newpop(T).Position(1,3)=VarMinsigma;
            
        elseif newpop(T).Position(1,3)>VarMaxsigma
            
            newpop(T).Position(1,3)=VarMaxsigma;
        end
        
        % Evaluation
        [newpop(T).Cost,newpop(T).sol]=CostFunction(newpop(T).Position);
        
    end
    
    % Merge Main Population (Archive) and New Population (Samples)
    pop=[pop
         newpop]; %#ok
    
    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Delete Extra Members
    pop=pop(1:nPop);
    
    % Update Best Solution Ever Found
    BestSol=pop(1);
    WorstSol=pop(end);
    
    % Store Mean Cost
    MCPop=[pop.Cost];
    MeanCost(it)=mean(MCPop);
    
    % Store Best Cost
    BestCost(it)=BestSol.Cost;
    
    % Store Worst Cost
    WorstCost(it)=WorstSol.Cost;
    
    %BEST
    Bestepsilon(it)=BestSol.Position(1,1);
    BestC(it)=BestSol.Position(1,2);
    Bestsigma(it)=BestSol.Position(1,3);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    disp(['Best Value for epsilon = ' num2str(Bestepsilon(it)) ': Best Value for C = ' num2str(BestC(it)) ': Best Value for sigma = ' num2str(Bestsigma(it))]);
    toc
end
%% Results

y=BestSol.sol.y;
epsilon=BestSol.Position(1,1);

% mat=[x',y'];
% mat2=sortrows(mat,1);
mat_plus=[(1:n)',(y(1,:)+epsilon)'];
mat2_plus=sortrows(mat_plus,1);
mat_min=[(1:n)',(y(1,:)-epsilon)'];
mat2_min=sortrows(mat_min,1);
save('Rasht_T_max')
figure;
plot((1:n),t,'o');
hold on;
plot((1:n),y,'k','LineWidth',2);
plot(mat2_plus(:,1),mat2_plus(:,2),'r:','LineWidth',2);
plot(mat2_min(:,1),mat2_min(:,2),'r:','LineWidth',2);
grid on;