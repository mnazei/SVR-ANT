function model=CreateModel()

X=xlsread('Babolsar.xlsx','TOTAL','c2:d769');
x=X';

T=xlsread('Babolsar.xlsx','TOTAL','h2:h769');
t=T';
n=numel(t);

model.x=x;
model.t=t;
model.n=n;

end
