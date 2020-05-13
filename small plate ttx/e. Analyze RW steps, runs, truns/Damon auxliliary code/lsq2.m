function [p,s]=lsq2(FUN,p0,x,y,e,varargin);


extravars{1}=FUN;
extravars{2}=x;
extravars{3}=y;
extravars{4}=e;

opts=optimset;
[p,chi2]=fminsearch(@sqerr,p0,opts,extravars);
s.chi2=chi2;
s.df=length(x)-length(p);
P_fit=chisqp(chi2,s.df);
P_fit=2*min(P_fit,1-P_fit);
s.p_fit=P_fit;

if length(varargin)>1
    NEWERROR=varargin{2};
else
    NEWERROR=0;
end

if NEWERROR
    for i=1:100
        [pstore(i,:),dum]=lsq2(FUN,p,x,y+e.*randn(size(e)),e,0,0);    % go through 100 times to find distribution of parameters
    end
    s.perr=pstore;
else
    warning off MATLAB:divideByZero
	perr=hessian2(@sqerr,p,p/100,extravars);
	perr=[sqrt(1./diag(perr))];
	s.perr=perr;
end


%% add in stuff here for the error in p -- hessian, etc.

if length(varargin)
    if varargin{1}>0
        if varargin{1}==2
            figure;
        end
        errorbar(x,y,e);
        hold on;
        yfit=feval(FUN,p,x);
        plot(x,yfit,'k-');
    end
end


function fiterr=sqerr(p,extra);

FUN=extra{1};
x=extra{2};
y=extra{3};
e=extra{4};

yfit=feval(FUN,p,x);
fiterr=sum((y-yfit).^2./e.^2);