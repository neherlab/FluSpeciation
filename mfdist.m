function mfdist(replica,mi,ni)
% the code simulates the binned traveling-wave dynamics with drift velocity
% equals to the total fraction of infection
% the function takes in three variables
% replica labels different rounds of simulations with the same parameter
% mi gives which of the mutation rates is using
% ni gives which of the host population size is using
%% parameters to be specified by input
d0   = 0.01; % inverse of cross-immunity range, antigenic advance for each mutation
m01  = 10.^(-3)/4;  % mutation rate % (d0)*0.01*1.5; %
n01  = 10.^(-9)/16; % inverse of host population % (d0)^2*10.^(-3);
m0   = m01(mi);  %
n0   = n01(ni);  %
%% parameters
h = 1; % forward simulation time interval *floor(0.1/d0); %min(1,.01/d0);% generation
T = 2000000; % total simulation period % *floor(0.1/d0);
totStep = floor(T/h)+1;  % total number of steps

Nrec = 10000;  % number of recording steps
hT   = floor(totStep/Nrec);  % steps before each record

p0   = 1000;   % threshold for generating random numbers with poisson or gaussian
nmax = round(1/d0);  % maximal number of bins, with maximal fitness 1
fitx = d0*(-nmax:nmax);  % fitness at each bin

s = rng('shuffle'); % seed for random numbers
%% initiate the problem
pop  = zeros(1,2*nmax+1);
%{
pop  = 100/sqrt(2*pi)*exp(-fitx.^2/2/d0^2);
fl   = pop>p0;
fs   = pop>0&pop<=p0;
pop(fl) = pop(fl)+sqrt(pop(fl)).*randn(1,sum(fl));
pop(fs) = poissrnd(pop(fs));
pop(pop<=0) = 0;
%}
pop(nmax+1) = 100; % initiate with the bin close to 0
if sum(pop)>0
    mfit = -d0/2; % shift the mean fitness to start the initial pop growth
end

%% variables to record the data

totI  = zeros(1,Nrec);  % total population
sm = zeros(1,Nrec);     % nose fitness
ms = zeros(1,Nrec);     % mean fitness
vs = zeros(1,Nrec);     % variance of fitness
s3 = zeros(1,Nrec);     % mean fitness
s4 = zeros(1,Nrec);     % variance of fitness
mx = zeros(1,Nrec);     % mean fitness

ext = 0;     % time of extinction event
extmp = 0;   % last time of extinction
et    = [];  % record all extinction time intervals
dp = zeros(1,2*nmax);  % change in pop in each bins except the last bin
%% evolve the system
tic;
for t = 2:totStep
    % Forward Eular
    pop = pop.*(1+h*(fitx-mfit)); % reproduction
    pop(pop<=0) = 0;  % no negative
    %{
    flg = find(pop>0);
    dp = poissrnd(m0*h*pop(flg(end))); %+sqrt(m0*h*pop(1:end-1)).*randn(1,2*nmax);
    if flg(end)<2*nmax+1
    pop(flg(end))   = pop(flg(end))-dp;
    pop(flg(end)+1) = dp;  % mutation
    end
    %}
    % fluctuation in mutation
    mp = m0*h*pop(1:end-1);
    p1 = mp>p0;
    dp(p1) = mp(p1)+sqrt(mp(p1)).*randn(1,sum(p1));
    dp(~p1)= poissrnd(mp(~p1));
    pop(1:end-1) = pop(1:end-1)-dp;
    %pop(2:end) = pop(2:end)+dp;
    if mfit>0
        pop(nmax+2:end) = pop(nmax+2:end)+dp(nmax+1:end);
    else
        pop(nmax+1:end) = pop(nmax+1:end)+dp(nmax:end);
    end
    % fluctuation in birth death process
    fl   = pop>p0;
    fs   = pop>0&pop<=p0;
    pop(fl) = pop(fl)+sqrt(pop(fl)).*randn(1,sum(fl));
    pop(fs) = poissrnd(pop(fs));
    pop(pop<=0) = 0;
   
    % mean fitness shifts by total fraction of infection
    mfit = mfit+h*(sum(pop)*n0);
    % adjust the mean fitness with keeping the middle bin near zero fitness
    if mfit>d0/2
        ii = round(mfit/d0);
        pop(1:end-ii) = pop(ii+1:end);
        pop(end-ii+1:end) = 0;
        mfit = mfit-d0*ii;
    end
    
    % record
    if mod(t-1,hT)==0
        tt = (t-1)/hT;
        totI(tt) = sum(pop);
        if totI(tt)>0
            ms(tt)   = sum((fitx-mfit).*pop)/totI(tt);
            sm(tt)   = max(fitx(pop>0))-mfit;
            vs(tt)   = sum((fitx-mfit).^2.*pop)/totI(tt)-ms(tt)^2;
            s3(tt)   = sum((fitx-mfit-ms(tt)).^3.*pop)/totI(tt);
            s4(tt)   = sum((fitx-mfit-ms(tt)).^4.*pop)/totI(tt);
            mx(tt)   = mfit;
        end
    end
    
    % extinction
    if sum(pop)==0
        if ext == 0
            ext   = (t-1)*h;
        else
            ext   = (t-1)*h;
            et = [et,ext-extmp];
        end
        %break;
        % reseeding
        pop(nmax+1) = poissrnd(1/d0);
        mfit = -d0/2;
        extmp = ext;
    end

end
toc;
et = [et,(t-1)*h-extmp];
%% save data
%dirc  = './mfdist/';
wname = 'Epidemics';   % to record population dynamics
ename = 'Extinct';     % to record extinctions
uname = sprintf('%.3f',0);
mname = sprintf('%.3f',m0*10000);
dname = sprintf('%.3f',d0);
rname = sprintf('%03d',replica);
nname = sprintf('%d',-log10(n0));
pname = sprintf('%d',p0);
dtype = '.dat';
epiname = [dirc,wname,'_',uname,'_',mname,'_',dname,'_',nname,'_',pname,'_',rname,dtype]; %,'_',pname
etname  = [dirc,ename,'_',uname,'_',mname,'_',dname,'_',nname,'_',pname,'_',rname,dtype];

dlmwrite(epiname,[h*hT:h*hT:T;totI;ms;vs;sm;mx;s3;s4]);
dlmwrite(etname,et);