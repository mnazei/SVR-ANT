function [RMSE,sol]=MyCost(X,model)

n=model.n;
x=model.x;
t=model.t;

epsilon=X(1,1);
C=X(1,2);
sigma=X(1,3);

Kernel=@(xi,xj) exp(-1/(2*sigma^2)*norm(xi-xj)^2);

H=zeros(n,n);
for i=1:n
    for j=i:n
        H(i,j)=Kernel(x(:,i),x(:,j));
        H(j,i)=H(i,j);
    end
end

HH=[ H -H
    -H  H];

f=[-t'; t']+epsilon;

Aeq=[ones(1,n) -ones(1,n)];
beq=0;

lb=zeros(2*n,1);
ub=C*ones(2*n,1);

options=optimset('MaxIter',200,'Display','off');

alpha=quadprog(HH,f,[],[],Aeq,beq,lb,ub,[],options);

alpha=alpha';

AlmostZero=(abs(alpha)<max(abs(alpha))*1e-4);

alpha(AlmostZero)=0;

alpha_plus=alpha(1:n);
alpha_minus=alpha(n+1:end);

eta=alpha_plus-alpha_minus;

S=find(alpha_plus+alpha_minus>0 & alpha_plus+alpha_minus<C);

y=zeros(size(t));
for i=1:n
    y(i)=MySVRFunc(x(:,i),eta(S),x(:,S),Kernel);
end

b=mean(t(S)-y(S)-sign(eta(S))*epsilon);
y=y+b;
%% Results
% yy=zeros(size(x));
% for k=1:n
%     yy(k)=MySVRFunc(x(:,k),eta(S),x(:,S),Kernel)+b;
% end

SSE=sse(y-t);
MSE=mse(y-t);
RMSE=sqrt(mse(y-t));

sol.MSE=MSE;
sol.RMSE=RMSE;
sol.SSE=SSE;
sol.y=y;

end
