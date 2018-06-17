function pop2tw(replica,mi,di)
% the code simulates the binned traveling-wave dynamics with two
% populations that are separating with less common change in fitness in t
% the code is adapted from the mean-field bin simulation mfdist.m
% the function takes in three variables
% replica labels different rounds of simulations with the same parameter
% mi gives which of the mutation rates is using
% di gives which of fitness changes is using
%replica = 1;mi = 12;di = 1;
%% parameters to be specified by input
s01  = [0.05,0.04,0.03,0.02]; % fitness difference for neighbor bins
n01  = 10.^(-4)*s0^2;   % host population
m01  = s0*0.005*[0.01,0.02,0.05,0.1,0.2,0.5,(1:12),(14:2:24)]; % mutation rate
s0   = s01(di); %
n0   = n01(1);  %
m0   = m01(mi);  %
%% parameters
ntot = 100; % total number of experiments
tcros= max(1,log(s0/m0))/s0^2;

h = 1.; %*floor(0.1/d0); %min(1,.01/d0); *floor(0.1/d0)  % generation
T = 2*round(tcros); %1000*round(1./s0); % maximum simulation
Tburn = 0; %round(10/s0);  % burn time roughly one sweeping time
totStep = floor(T/h)+1;
bnStep  = floor(Tburn/h);

Nrec = T;
hT   = floor(totStep/Nrec);

p0   = 1000;
nmax = 30;
d0   = s0*log(2);
fitx = d0*(-nmax:nmax);

s = rng('shuffle');
%% variables to record the data
%
totI  = zeros(2,Nrec);  % total population
sm = zeros(2,Nrec);     % nose fitness
ms = zeros(2,Nrec);     % mean fitness
vs = zeros(2,Nrec);     % variance of fitness
s3 = zeros(2,Nrec);     % mean fitness
s4 = zeros(2,Nrec);     % variance of fitness
mx = zeros(2,Nrec);     % mean fitness
%}
xn   = zeros(1,ntot);
mphi = zeros(1,ntot);

extmp = 0;
et    = [];
bt    = [];
scount= 0; % successful
ncount= 0;
%% evolve the system
tic;
for ni = 1:ntot
    flag = 1;
    %% initiate the problem
    pop  = zeros(2,2*nmax+1);  % two populations
    mfit = zeros(2,1);
    pop(1,:) = 0.1*d0^2/n0/sqrt(2*pi)*exp(-fitx.^2/2/d0^2);
    pop(2,:) = 0.1*d0^2/n0/sqrt(2*pi)*exp(-fitx.^2/2/d0^2);
    fl   = pop>p0;
    fs   = pop>0&pop<=p0;
    pop(fl) = pop(fl); %+sqrt(pop(fl)).*randn(1,sum(fl));
    pop(fs) = poissrnd(pop(fs));
    pop(pop<=0) = 0;
    for i = 1:2
        mfit(i) = (sum(fitx.*pop(1,:))+sum(fitx.*pop(2,:)))/sum(sum(pop));
    end
    dp = zeros(1,2*nmax);

    % burn time two populations fully interacting with each other 
    for t = 1:bnStep
        % forward Eular for two populations
        for i = 1:2
            pop(i,:) = pop(i,:).*(1+h*(fitx-mfit(i)));  %reproduction
            mp = m0*h*pop(i,1:end-1);
            p1 = mp>p0;
            dp(p1) = mp(p1)+sqrt(mp(p1)).*randn(1,sum(p1));
            dp(~p1)= poissrnd(mp(~p1));
            pop(i,1:end-1) = pop(i,1:end-1)-dp;
            pop(i,2:end) = pop(i,2:end)+dp;
            fl   = pop>p0;
            fs   = pop>0&pop<=p0;
            pop(fl) = pop(fl)+sqrt(pop(fl)*h).*randn(1,sum(fl));
            pop(fs) = poissrnd(pop(fs));
            pop(i,pop(i,:)<1) = 0;
        end
        
        % adapt the mean fitness for the two populations
        for i = 1:2
            mfit(i) = mfit(i)+h*((sum(pop(i,:))+sum(pop(mod(i,2)+1,:)))*n0);
        end
        if max(mfit)>d0/2
            ii = round(max(mfit)/d0);
            for i = 1:2
                pop(i,1:end-ii) = pop(i,ii+1:end);
                pop(i,end-ii+1:end) = 0;
            end
            mfit = mfit-d0*ii;
        end
        
        % one branch dies during burning
        if any(sum(pop,2)==0)
            flag = 0;
            break;
        end
    end
    
    if flag  % if two populations both exist, their interaction starts to fade away
        ncount = ncount+1;
        for t = 2:totStep
            % forward Euler
            %
            for i = 1:2
                pop(i,:) = pop(i,:).*(1+h*(fitx-mfit(i)));  %reproduction
                mp = m0*h*pop(i,1:end-1);
                p1 = mp>p0;
                dp(p1) = mp(p1)+sqrt(mp(p1)).*randn(1,sum(p1));
                dp(~p1)= poissrnd(mp(~p1));
                pop(i,1:end-1) = pop(i,1:end-1)-dp;
                pop(i,2:end) = pop(i,2:end)+dp;
                fl   = pop>p0;
                fs   = pop>0&pop<=p0;
                pop(fl) = pop(fl)+sqrt(pop(fl)*h).*randn(1,sum(fl));
                pop(fs) = poissrnd(pop(fs));
                pop(i,pop(i,:)<1) = 0;
            end
            
            % adapt the mean fitness, drift velocity by the other
            % population decays in time
            for i = 1:2
                mfit(i) = mfit(i)+h*(sum(pop(i,:))+exp(-2*h*t/tcros)*sum(pop(mod(i,2)+1,:)))*n0;
            end
            if max(mfit)>d0/2
                ii = round(max(mfit)/d0);
                for i = 1:2
                    pop(i,1:end-ii) = pop(i,ii+1:end);
                    pop(i,end-ii+1:end) = 0;
                end
                mfit = mfit-d0*ii;
            end
            
            % record
            %
            if mod(t-1,hT)==0
                tt = (t-1)/hT;
                for i = 1:2
                    totI(i,tt) = sum(pop(i,:));
                    %if mod(t-1,1000)==0
                    %    plot(fitx,pop);
                    %end
                    if totI(i,tt)>0
                        ms(i,tt)   = sum((fitx-mfit(i)).*pop(i,:))/totI(i,tt);
                        sm(i,tt)   = max(fitx(pop(i,:)>0))-mfit(i);
                        vs(i,tt)   = sum((fitx-mfit(i)).^2.*pop(i,:))/totI(i,tt)-ms(i,tt)^2;
                        s3(i,tt)   = sum((fitx-mfit(i)-ms(i,tt)).^3.*pop(i,:))/totI(i,tt);
                        s4(i,tt)   = sum((fitx-mfit(i)-ms(i,tt)).^4.*pop(i,:))/totI(i,tt);
                        mx(i,tt)   = mfit(i);
                    end
                end
            end
            %}
            
            % successful speciation, both pop survives to fully lose the
            % interaction
            if 2*h*t>tcros
                scount = scount+1;
                break;
            end
            
            % one branch dies during burning
            if any(sum(pop,2)==0)
                brt = (t-1)*h;
                bt = [bt,brt];
                % extinction both branches die
                if all(sum(pop,2))==0
                    ext   = (t-1)*h;
                    if (t-1)*h>T/10
                        et    = [et,ext-extmp];
                    end
                end
                break;
            end
            
        end
        if tt>1
            xn(ni) = mean(max(sm(:,1:tt)));
            mphi(ni) = exp(mean(log(sum(totI(:,1:tt)))));
        end
    end
end
toc;

%% save data 
wname = 'Epidemics';   % to record the pop dyn
ename = 'Extinct';     % to record the extinctions
sname = 'Success';     % to record the number of successful branching
uname = sprintf('%.3f',0);
mname = sprintf('%.3f',m0*10000);
dname = sprintf('%.3f',d0);
rname = sprintf('%03d',replica);
nname = sprintf('%.2f',-log10(n0));
%pname = sprintf('%d',p0);
dtype = '.dat';
epiname = [wname,'_',uname,'_',mname,'_',dname,'_',nname,'_',rname,dtype]; %,'_',pname
etname  = [ename,'_',uname,'_',mname,'_',dname,'_',nname,'_',rname,dtype];
scname  = [sname,'_',uname,'_',mname,'_',dname,'_',nname,'_',rname,dtype];

dlmwrite(epiname,[xn;mphi]); %h*hT:h*hT:T;totI;ms;vs;sm;mx;s3;s4]);
dlmwrite(etname,bt);
dlmwrite(scname,[scount,ncount]);
