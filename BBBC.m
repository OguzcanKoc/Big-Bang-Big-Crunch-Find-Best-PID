function [bestparams, bestval]=BBBC(paramslb, paramsub, objfunc, popnum, ds, sc)

D=length(paramslb);
range0=paramsub-paramslb;
x0=zeros(D, popnum);
for i=1:popnum
    x0(:, i)=paramslb+(paramsub-paramslb)/2;
end
x=x0;
bestparams=x(:, 1);
fhd=str2func(objfunc);
fvals=feval(fhd, x, 0, 0);
% fvals=deneme_1_BBBC(x);
bestval=min(fvals);
range=range0;
normrange=norm(range);
stpcrtr=range*sc/100;
decscl=ds/100;

it=0;
while range>=stpcrtr
    range=range-range*decscl;
    it=it+1;
end
maxiter=it;

range=range0;
it=1;
while range>=stpcrtr
    str=sprintf('iter: %d/%d, bestval: %g, range norm: %g', it, maxiter, bestval, normrange);
    disp(str)
    
    for i=1:popnum
        x(:, i)=bestparams+2*(rand(D, 1)-0.5).*range;
        ind=find(x(:, i)<paramslb);
        x(ind, i)=paramslb(ind);
        ind=find(x(:, i)>paramsub);
        x(ind, i)=paramsub(ind);
    end
    fhd=str2func(objfunc);
    fvals=feval(fhd, x, it, maxiter);
    % fvals=deneme_1_BBBC(x);
    if min(fvals)<bestval
        bestval=min(fvals);
        ind=find(fvals==min(fvals));
        ind=min(ind);
        bestparams=x(:, ind);
    end
    
    range=range-range*decscl;
    normrange=norm(range);
    it=it+1;
end