function FluEpiTreeNM(replica,mui,mi,di,ni)
% the code simulates the SIR type epidemics with antigentic mutation for the flu
% the function takes in five variables
% replica labels different rounds of simulations with the same parameter
% mui gives which of the mortality rates is using
% mi gives which of the mutation rates is using
% di gives which of the deleterious effect is using
% ni gives which of the host population sizes is using
%replica=1;mui=1;mi=12;di=4;ni=5;
%% parameters to be determined by input
mu01 = 0.;  % no mortal
d01  = 0.02;  % inverse of range of cross-immunity, antigentic advance in one mutation 
m01  = d01*0.01*(1:24); % mutation rates 
%[0.002:0.002:0.01,(1:7)*0.02]; %(4:15); %(10:2:32); %[0.02:0.01:0.08,0.1:0.02:0.18]; %(1:0.5:6.5); %[(0.002:.002:0.01),(0.015:0.005:0.045)]; %[(0.06:0.02:0.18),0.2:0.05:0.4];  %0.001*mu01;
b01  = [1,0.5,0.2,0.1,0.01];  % deleterious effect in transmissibility
n01  = 10.^(-8:-3)*d01^2;   % inverse of host population sizes
mu0  = mu01(mui);  % mortality
m0   = m01(mi);  % mutation rate
d0   = d01;  % inverse of cross-immunity range
b0   = d0*b01(di);  % deleterious change
n0   = n01(ni);  % inverse of host population
%replica = 0;
%% simulation parameters
h = min(0.5,0.01*round(0.1/d0)); % forward time inteval
T = 50000*round(0.1/d0); % total simulation period
totStep = floor(T/h)+1;  % total number of forward steps
lT = round(T/100);       % maximal number of strains to quit the simulation
nm = min(30,max(1,round(m0*d0./n0/200)));

nmax = 10;  % maximal number of decendent strains of a given strain (to reduce computation cost)
N = nm*400; % maximal number of alive genotype
Nrec = 50000; % number of steps to be recorded
hT   = floor(totStep/Nrec);  % step inteval for recording the data

sig  = 1;    % range of cross-immunity
beta0= 10;   % transmissibility 
db   = beta0*b0;  % deleterious effect in transmissibility
nu   = 5;   % recover rate from infection
mort = mu0; % mort rate could be strain dependent
mu   = m0;  % mutation rate per individual
I0   = n0;  % minimal fraction of infectious (inverse of effective population, epidemics starts from this value)

s = rng('shuffle');  % random seeding
%% initialize the problem
%nid   = zeros(1,N*N); % index on the alive
prts  = [];  % parents of genotype
%taus  = zeros(1,N);  % as a common ancestor
birth = [];  % birth time of genotype
death = [];  %zeros(lT*N,1);  % death time of genotype
noff  = [];  % no of alive offspring lineages
droot = [];  % distance to root
Imax  = [];  % maximum I of strain
talp  = [];  % time to the latest offsprings
beta  = [];  % transmissibility

xx    = [];  % fitness at birth
dx    = [];  % antigenic gain from mutation
phi   = [];  % local cross-immunity at born

Kij  = (eye(N));  % antigenic coupling
Spop = zeros(1,N); % population fraction of susceptible
Ipop = zeros(1,N); % population of infected
nn   = zeros(1,N); % number of offsprings
gid  = zeros(1,N); % id of genotype

% initialize the strain
gn   = 1;     % label the first seeding strain
iniI = d0^2;  % initial fraction of infectious
iniq = log(beta0/nu)-d0;   % initial immunity (close to the steady state to avoid fast extinction)
beta(gn) = beta0;  % initial transmissibility
prts(gn) = 0;   % no parents
birth(gn) = 0;  % born at t=0
death(gn) = 0;  % not yet dead
droot(gn) = 0;  % is the root
xx(gn)    = beta(gn)/nu*exp(-iniq);  % fitness at birth
dx(gn)    = 0;  % fitness gain at birth
phi(gn)   = (mort/nu*(log(beta(gn)/nu)));  % drift at birth
Imax(gn)  = iniI;  % maximal infection for strains
noff(gn) = 1;   % no of alive offspring branches
%id = noff>1;   % index of the strain with more than one branch
for i = 1:gn
    Ipop(i) = iniI;  
    Spop(i) = 1;  % initial Spop will be updated later.
    gid(i)  = i;
end
alive = find(Ipop>0);  % alive strains

%% variables to record the data

mdist = zeros(1,Nrec);  % shortest route to the root

totR  = zeros(1,Nrec);  % total number of resist inividuals at each time
totI  = zeros(1,Nrec);  % tot number of infectious
nstr = zeros(1,Nrec);   % number of alive strains
ms = zeros(1,Nrec);   % mean number of susceptible individuals
vs = zeros(1,Nrec);   % variance
sm = zeros(1,Nrec);   % maximum
%s3 = zeros(1,Nrec);   % 3rd moment
%s4 = zeros(1,Nrec);   % 4th moment
mb = zeros(1,Nrec);   % mean beta
vb = zeros(1,Nrec);   % variance of beta
mf = zeros(1,Nrec);   % mean number of cross-immunity
vf = zeros(1,Nrec);   % variance
md = zeros(1,Nrec);   % mean distance to the root
vd = zeros(1,Nrec);   % variance
mr = zeros(1,Nrec);   % mean cross-immunity distance
vr = zeros(1,Nrec);   % variance
ts = zeros(1,Nrec);   % total time distance to the common root
tp = zeros(1,Nrec);   % time increment to the last common root

mlca = zeros(1,Nrec);   % time to the last common root
dlca = zeros(1,Nrec);   % generation to the last common root

% record extinction time and waiting time for successive branching
et    = []; % extinctions
wt    = []; % wait time
et_br = []; % extinction of branches
ngbr  = []; % number of generations in branches
naa   = []; % number of alive strains
net   = 0;
nwt   = 0;

% record the initial status

Sp = 1 - sum(Ipop);  % total number of susceptibles
for j = alive
    Kij((1:gn),(1:gn)) = 1;
    Spop(j) = Sp*exp(-sum(iniq*Kij((1:gn),(1:gn))));
end
nalv = 1;   % number of alive strains

%flag  = 0;
ext   = 0;  % extinction time
extmp = 0;  % last extinction time
tautp0= 0;  % if extinct, corrected by the time to the reseeding strain
taus = 0;   % time to the MRCA
taup = 0;   % generations to MRCA
%% evolve the system
tic;
for t = 2:totStep
    % if the number of coexisting strains or total number of emerged
    % strains larger than the recording threshold
    if length(alive)>=N || gn>=lT*N
        break;
    end
    
    Itemp = Ipop;
    Stmp0 = Spop;
    Stemp = Spop;
    tautp = taus;  % time to the MRCA
    taupp = taup;
    %idtmp = id;
    
    % Midpoint method for forward SIR simulation
    Itemp(alive) = Ipop(alive)+h/2*(beta(gid(alive)).*Spop(alive)-nu-mort).*Ipop(alive);
    dItmp = - sum(Itemp(alive)-Ipop(alive));
    Stemp(alive) = Spop(alive)+dItmp/Sp*Spop(alive)-h/2*(nu*Ipop(alive)*Kij(alive,alive)-mort*log(Spop(alive))).*Spop(alive);
    Sptmp = Sp+dItmp;
    
    dI = h*(beta(gid(alive)).*Stemp(alive)-nu-mort).*Itemp(alive);
    Ipop(alive) = Ipop(alive)+dI;
    Spop(alive) = Spop(alive)-sum(dI)/Sptmp*Stemp(alive)-h*(nu*Itemp(alive)*Kij(alive,alive)-mort*log(Stemp(alive))).*Stemp(alive);
    Sp = 1 - sum(Ipop);
    
    % extinction
    % remove the strains from alive set
    id = find(Ipop(alive)<n0);
    Ipop(alive(id)) = 0;    
    nn(alive(id)) = 0;
    flg  = death(gid(alive))==0 & Ipop(alive)<n0;
    death(gid(alive(flg)))= (t-1)*h;
    for ii = id
        i = gid(alive(ii));
        noff(i) = noff(i)-1;
        if noff(i)<=0
            talp(i) = (t-1)*h;
            noff(i) = 0;
        end
        while prts(i)>0 && noff(i)==0
            i = prts(i);
            noff(i) = noff(i)-1;
            if noff(i)<=0
                talp(i) = (t-1)*h;
                noff(i) = 0;
            end
        end
    end
    alive(id) = [];
    
    if isempty(alive)  % no strain left all extinct
        % not record the first extinction
        if ext == 0
            ext   = (t-1)*h;
        else
            ext   = (t-1)*h;
            net   = net+1;
            et(net)= ext-extmp;
        end
        % reseeding with a previous strain
        [~,i] = max(Spop);
        Ipop(i) = 10*I0;
        %totI(t) = 10*I0;
        alive = i;
        Spop(i) = Spop(i)^(exp(-d0/sig));
        noff(gid(i)) = noff(gid(i))+1;
        extmp = ext;
        tautp0 = (t-1)*h-birth(gid(i));
    end
    
    % update the time to the MRCA
    id = find(noff>1);
    if ~isempty(id)
        taus = max((t-1)*h-birth(id)-tautp0);
        taup = max(droot(gid(alive)))-min(droot(id));
    else
        taus = max((t-1)*h-birth(gid(alive))-tautp0);
        taup = max(droot(gid(alive)))-min(droot(gid(alive)));
    end
    
    % if time to MRCA jumps
    if taus < tautp  % successive branching point
        if abs(taus)>0.9*h %any(id)
            taus = taus+tautp0;
            tautp0 = 0;
        end
        nwt = nwt+1;
        wt(nwt) = tautp-taus;  % jump in tMRCA, time between successive MRCA
        et_br(nwt) = tautp;  % depth in time of the dead branch
        ngbr(nwt)  = taupp;  % depth in generations of the dead branch
        naa(nwt)   = nalv;   % number of surviving strains
    end
    
    % add time interval to all strains with alive offsprings
    % mutation (stochastic)
    for j=alive
        if rand(1)<mu*Ipop(j)/I0*h  && nn(j)<nmax  % if mutation happens by chance
            dd    = d0; % distance to its root (can be draw from some distribution)
            if rand(1)<(1-exp(1-beta0/beta(j)))
                beta1 = beta0;  % compensatory mutation to transmissibility
            else
                beta1 = beta(j)-db; %deleterious effect on transmissibility of an antigenic mutation *log(rand(1));
            end
            SnewP = Sp*(Spop(j)/Sp)^(exp(-dd/sig)); % number of susceptibles to the new strain
            if beta1*SnewP>nu  % the new strain will survive
                % adding new strain and its potential genotype
                if length(alive)>=N || gn>=lT*N
                    break;
                end
                % update tree
                gn = gn+1; % a new genotype added
                nal = setdiff(1:N,alive); % not used ids
                if j>1
                    nn(j)     = nn(j)+1; 
                    % number of offsprings considered is smaller than a fixed number
                end
                prts(gn)  = gid(j);  % parents
                birth(gn) = (t-1)*h;  % birth time
                death(gn) = 0;  % death
                talp(gn)  = 0;  % time to the last alive
                beta(gn)  = beta1;
                droot(gn) = droot(gid(j))+1;  % depth to root
                xx(gn)    = beta1*(SnewP)/nu;   % fitness at birth
                dx(gn)    = beta1*(SnewP-Spop(j))/nu;  % fitness gain at birth
                phi(gn)   = (beta1/nu*(Stmp0(j)-Spop(j))/h/nu+mort/nu*(log(beta1/nu)-beta1/nu*Stmp0(j)+1));  % drift at birth
                Imax(gn)  = I0;
                noff(gid(j)) = noff(gid(j))+1;  % add an offspring
                noff(gn)  = 1; % newly born strain
                %ancs{gn} = [ancs{j},gn];
                
                % update epidemic equations for the new strain
                j1 = nal(1);
                Kij(j1,alive)= Kij(j,alive)*exp(-dd/sig);
                Kij(alive,j1)= Kij(alive,j)*exp(-dd/sig); % update the distance in all genotype
                Ipop(j1) = I0;
                Spop(j1) = SnewP;  % initial Spop will be updated later.
                alive = union(alive,j1);  % alive strains% immunity types
                gid(j1) = gn;
            end
            Ipop(j) = Ipop(j)-I0; % failed mutation
        end
        if Ipop(j)>Imax(gid(j)) % record the maximal infection for a strain.
            Imax(gid(j)) = Ipop(j);
        end
    end
    
    nalv = length(alive); 
    
%
    % extinguish the branches other than the most populated one to root if
    % survives longer than the cross-immunity range.
    if taup>sig/d0 || length(alive)>=N
        id = find(noff>1);
        if length(id)>1
            [~,ord] = min(birth(id));
            r0 = id(ord);  % current root
            brchs = {};  % strains in brch
            balv  = {};  % alive in brch
            godie = [];
            nr    = 0;
            for i = alive
                j = gid(i);
                pp = j;
                flg = 0;
                if nr>0
                    for ri=1:nr
                        if ismember(j,brchs{ri})
                            flg = ri;
                        end
                    end
                end
                while j>r0 && flg==0
                    j = prts(j);
                    pp(end+1) = j;
                    if nr>0
                        for ri=1:nr
                            if ismember(j,brchs{ri})
                                flg = ri;
                            end
                        end
                    end
                end
                if flg==0
                    dp = setdiff(pp,r0);
                    if ~isempty(dp)
                        nr = nr+1;
                        brchs{nr} = dp;
                        balv{nr} = i;
                    else
                        godie = i;
                    end
                else
                    brchs{flg} = union(brchs{flg},pp);
                    balv{flg} = [balv{flg},i];
                end
            end  % find different branches
            rr = zeros(1,nr);
            bI = zeros(1,nr);
            for ri = 1:nr
                bid = intersect(id,brchs{ri});
                if isempty(bid)  % no branch along the lineage
                    rr(ri) = max(brchs{ri});
                else
                    rr(ri) = min(bid);
                end
                bI(ri) = sum(Ipop(balv{ri}));
            end
            [~,mid] = max(bI);
            p0 = [];
            
            for ri = 1:nr
                if ri~=mid
                    p0 = union(p0,brchs{ri});
                    godie = union(godie,balv{ri});
                end
            end
            Ipop(godie) = 0;
            nn(godie)   = 0;
            death(gid(godie))= (t-1)*h;
            noff(p0) = 0;  % no offsprings in the path
            talp(p0) = (t-1)*h;
            noff(r0) = 1;
            alive = setdiff(alive,godie);
        end
    end
%
    
    % record data
    if mod(t-1,hT)==0
        tt = (t-1)/hT;
        % infections
        totI(tt) = sum(Ipop(alive));
        nstr(tt) = length(alive);
        % susceptible
        ms(tt)   = sum(beta(gid(alive))/nu.*Spop(alive).*Ipop(alive))/totI(tt);
        sm(tt)   = max(beta(gid(alive))/nu.*Spop(alive));
        vs(tt)   = sum((beta(gid(alive))/nu.*Spop(alive)).^2.*Ipop(alive))/totI(tt)-ms(tt)^2;
        %s3(tt)   = sum((beta(gid(alive))/nu.*Spop(alive)-ms(tt)).^3.*Ipop(alive))/totI(tt);
        %s4(tt)   = sum((beta(gid(alive))/nu.*Spop(alive)-ms(tt)).^4.*Ipop(alive))/totI(tt);
        mb(tt)   = sum(beta(gid(alive))/nu.*Ipop(alive))/totI(tt);
        vb(tt)   = sum((beta(gid(alive))/nu-mb(tt)).^2.*Ipop(alive))/totI(tt);
        % drift
        fftmp    = (beta(gid(alive))./nu.*(Stmp0(alive)-Spop(alive))/h/nu+mort/nu.*(log(beta(gid(alive))./nu)-beta(gid(alive))./nu.*Stmp0(alive)+1));
        mf(tt)   = sum(fftmp.*Ipop(alive).*(fftmp>0))/totI(tt);
        vf(tt)   = sum((fftmp-mf(tt)).^2.*Ipop(alive).*(fftmp>0))/totI(tt);
        % antigenic distance
        md(tt)   = sum(droot(gid(alive)).*Ipop(alive))/totI(tt);
        vd(tt)   = sum((droot(gid(alive))-md(tt)).^2.*Ipop(alive))/totI(tt);
        if length(alive)>1
            uptr     = triu(Kij(alive,alive),1);
            ddd      = -sig*log(reshape(uptr(uptr>0),1,[]));
            mr(tt)   = mean(ddd);  % mean distance between strains
            vr(tt)   = var(ddd);
        end
        mdist(tt)= min(droot(gid(alive)));
        % common ancestor
        ts(tt)   = taus;
        tp(tt)   = taup;
    end

end
talp(talp==0) = (t-1)*h;
net   = net+1;
et(net)= (t-1)*h-extmp;
toc

%% save data
%dirc  = './';
wname = 'Epidemics';  % to record epidemic data
tname = 'WaitTime';   % to record branching data
ename = 'Extinct';    % to record extinction data
gname = 'Genotype';   % to record strain data
uname = sprintf('%.3f',100*mu0);    % mortality rate
mname = sprintf('%.3f',m0*10000);   % mutation rate
dname = sprintf('%.3f',d0);         % antigenic advance
nname = sprintf('%.2f',-log10(n0)); % population size
sname = sprintf('%.2f',(d0/b0));    % deleterious effect in transmissibility
rname = sprintf('%03d',replica);    % replica of the same set of parameters
dtype = '.dat';
epiname = [wname,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];
wtname  = [tname,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];
etname  = [ename,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];
gnname  = [gname,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];

dlmwrite(epiname,[h*hT:h*hT:T;totI;totR;ms;vs;sm;nstr;mb;vb;mf;vf;ts;tp;mlca;dlca;mdist;md;vd]);
dlmwrite(wtname,[wt;et_br;ngbr;naa]);
dlmwrite(etname,et);
dlmwrite(gnname,[prts;birth;death;beta;talp;droot;dx;xx;Imax;phi]);
