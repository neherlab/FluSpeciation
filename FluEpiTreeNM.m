%function FluEpiTreeNM(replica,mui,mi,di,ni)
% the code simulates the SIR type epidemics with antigentic mutation for the flu
% the function takes in five variables
% replica labels different rounds of simulations with the same parameter
% mui gives which of the mortality rates is using
% mi gives which of the mutation rates is using
% di gives which of the deleterious effect is using
% ni gives which of the host population sizes is using
replica=1;mui=1;mi=24;di=7;ni=4;
% no mortal
mu01 = 0.; 
d01  = 0.02;
m01  = d01*0.01*(1:24); %[0.002:0.002:0.01,(1:7)*0.02]; %(4:15); %(10:2:32); %[0.02:0.01:0.08,0.1:0.02:0.18]; %(1:0.5:6.5); %[(0.002:.002:0.01),(0.015:0.005:0.045)]; %[(0.06:0.02:0.18),0.2:0.05:0.4];  %0.001*mu01;
b01  = [1,0.5,0.3,0.2,0.1,0.01,0];
mu0  = mu01(mui);  %(mui)
m0   = m01(mi);  %(mi)
d0   = d01;  %(di)
b0   = d0*b01(di);
n01  = 10.^(-8:-3)*d01^2;
n0   = n01(ni);
%replica = 0;
%% parameters
h = min(0.5,0.01*round(0.1/d0)); %0.02;
T = 200*round(0.1/d0); %max(5000,10000*round(n0/(d0^2*m0)))
totStep = floor(T/h)+1;
lT = round(T/100);
nm = min(30,max(1,round(m0*d0./n0/200)));

nmax = 10;
N = nm*400; % maximum number of alive genotype
Nrec = 5000;
hT   = floor(totStep/Nrec);

sig  = 1;   % range of cross-immunity
beta0= 10;   % transimission rate
db   = beta0*b0;
nu   = 5;   % recover rate from infection
mort = mu0;  % mort rate could be strain dependent
mu   = m0; % mutation rate per individual
I0   = n0;  % minimal fraction of infectious (inverse of effective population, epidemics starts from this value)

s = rng('shuffle');
%% initiate the problem
%nid   = zeros(1,N*N); % index on the alive
prts  = [];  % parents of genotype
%taus  = zeros(1,N);  % as a common ancestor
birth = [];  % birth time of genotype
death = []; %zeros(lT*N,1);  % death time of genotype
noff  = [];  % no of alive offspring lineages
droot = [];  % distance to root
Imax  = [];  % maximum I of strain
talp  = [];  % time to the latest offsprings
beta  = [];

xx    = [];  % fitness at birth
dx    = [];  % antigenic gain from mutation
phi   = [];  % local cross-immunity at born

Kij  = (eye(N));  % antigenic coupling
Spop = zeros(1,N); % population fraction of susceptible
Ipop = zeros(1,N); % population of infected
nn   = zeros(1,N); % number of offsprings
gid  = zeros(1,N); % id of genotype

%ancs = cell(N);

% initialize the strain
gn   = 1;   % label the number of appeared genotype
iniI = d0^2;  % initial fraction of infectious
iniq = log(beta0/nu)-d0;   % initial immune
beta(gn) = beta0;
prts(gn) = 0;
birth(gn) = 0;
death(gn) = 0;
droot(gn) = 0;  % depth to root
xx(gn)    = 1;   % fitness at birth
dx(gn)    = 0;  % fitness gain at birth
phi(gn)   = (mort/nu*(log(beta(gn)/nu)));  % drift at birth
Imax(gn)  = iniI;
noff(gn) = 1;  % no of alive offspring branches
%id = noff>1;    % index of the strain with more than one branch
for i = 1:gn
    Ipop(i) = iniI;
    Spop(i) = 1;  % initial Spop will be updated later.
    %nid(i)  = i;
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
vn = zeros(1,Nrec);   % velocity of nose

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

% record the initial 

Sp = 1 - sum(Ipop);  % total number of susceptible individuals
for j = alive
    Kij((1:gn),(1:gn)) = 1;
    Spop(j) = Sp*exp(-sum(iniq*Kij((1:gn),(1:gn))));
end
%ancs{1} = 1;
nalv = 1;

%flag  = 0;
ext   = 0;  % extinction time
extmp = 0;  % last extinction time
tautp0= 0;
taus = 0;
taup = 0;
dtmp = 0;
%% evolve the system
tic;
for t = 2:totStep
    
    if length(alive)>=N || gn>=lT*N
        break;
    end
    
    Itemp = Ipop;
    Stmp0 = Spop;
    Stemp = Spop;
    tautp = taus;  % time to the last common ancestor
    taupp = taup;
    %idtmp = id;
    % Midpoint
    %
    Itemp(alive) = Ipop(alive)+h/2*(beta(gid(alive)).*Spop(alive)-nu-mort).*Ipop(alive);
    dItmp = - sum(Itemp(alive)-Ipop(alive));
    Stemp(alive) = Spop(alive)+dItmp/Sp*Spop(alive)-h/2*(nu*Ipop(alive)*Kij(alive,alive)-mort*log(Spop(alive))).*Spop(alive);
    Sptmp = Sp+dItmp;
    
    dI = h*(beta(gid(alive)).*Stemp(alive)-nu-mort).*Itemp(alive);
    Ipop(alive) = Ipop(alive)+dI;
    Spop(alive) = Spop(alive)-sum(dI)/Sptmp*Stemp(alive)-h*(nu*Itemp(alive)*Kij(alive,alive)-mort*log(Stemp(alive))).*Stemp(alive);
    Sp = 1 - sum(Ipop);
    
    % extinction
    id = find(Ipop(alive)<n0);
    Ipop(alive(id)) = 0;    
    nn(alive(id)) = 0;
    flg  = death(gid(alive))==0 & Ipop(alive)<n0;
    death(gid(alive(flg)))= (t-1)*h;
    %if ~isempty(id)
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
    %end
    
    if isempty(alive)  % no strain left all extinct
        % not record the first extinction
        if ext == 0
            ext   = (t-1)*h;
        else
            ext   = (t-1)*h;
            net   = net+1;
            et(net)= ext-extmp;
        end
        % reservior
        [~,i] = max(Spop);
        Ipop(i) = 10*I0;
        %totI(t) = 10*I0;
        alive = i;
        Spop(i) = Spop(i)^(exp(-d0/sig));
        noff(gid(i)) = noff(gid(i))+1;
        extmp = ext;
        tautp0 = (t-1)*h-birth(gid(i));
    end
    
    % tau
    id = find(noff>1);
    if ~isempty(id)
        taus = max((t-1)*h-birth(id)-tautp0);
        taup = max(droot(gid(alive)))-min(droot(id));
    else
        taus = max((t-1)*h-birth(gid(alive))-tautp0);
        taup = max(droot(gid(alive)))-min(droot(gid(alive)));
    end
    
    
    if taus < tautp  % successive branching point
        if abs(taus)>0.9*h %any(id)
            taus = taus+tautp0;
            tautp0 = 0;
        end
        nwt = nwt+1;
        wt(nwt) = tautp-taus;
        et_br(nwt) = tautp;
        ngbr(nwt)  = taupp;
        naa(nwt)   = nalv; 
        %tautp0 = 0;
    end
    
    % add time interval to all strains with alive offsprings
    % mutation (stochastic)
    xn = max(beta(gid(alive))/nu.*Spop(alive))-1;
    vx = var(beta(gid(alive))/nu.*Spop(alive));
    for j=alive
        if rand(1)<mu*Ipop(j)/I0*h  && nn(j)<nmax
            dd    = d0; % distance to its root (can be draw from some distribution)
            if rand(1)<(1-exp(1-beta0/beta(j)))
                beta1 = beta0;
            else
                beta1 = beta(j)-db; %+db*log(rand(1));
            end
            SnewP = Sp*(Spop(j)/Sp)^(exp(-dd/sig));
            x1 = beta1/nu.*SnewP-1;
            if x1>0 && rand(1)<min(1,exp((x1^2-xn^2)/vx/2)) %establish beta1*SnewP>nu  % will survive
                % adding new genotype and its potential genotype
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
                Imax(gn)  = I0/x1;  % start from 1/x I0
                noff(gid(j)) = noff(gid(j))+1;  % add an offspring
                noff(gn)  = 1; % newly born strain
                %ancs{gn} = [ancs{j},gn];
                
                % update epidemics
                j1 = nal(1);
                Kij(j1,alive)= Kij(j,alive)*exp(-dd/sig);
                Kij(alive,j1)= Kij(alive,j)*exp(-dd/sig); % update the distance in all genotype
                Ipop(j1) = I0/x1;  % start from 1/x I0
                Spop(j1) = SnewP;  % initial Spop will be updated later.
                alive = union(alive,j1);  % alive strains% immunity types
                gid(j1) = gn;
            end
            Ipop(j) = Ipop(j)-I0;
            %totI(t) = totI(t)-I0; % failed mutation
        end
        if Ipop(j)>Imax(gid(j))
            Imax(gid(j)) = Ipop(j);
        end
    end
    
    nalv = length(alive); 
    
%
    % extinguish the branches other than the most populated one to root
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
    
    % record
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
        vn(tt)   = (max(droot)-dtmp)/h/hT;
        dtmp  = max(droot);
        % distance generation to the last common ancestor
        %{
        if ~isempty(alive)
            if isempty(id0) || ~isempty(setdiff(id0,alive))
                aset = ancs{alive(1)};
                for ii = alive
                    aset = intersect(aset,ancs{ii});
                end
            else
                for ii = setdiff(alive,id0)
                    aset = intersect(aset,ancs{ii});
                end
            end
            if isempty(mast) || mast~=max(aset)
                mast = max(aset);
                ml0  = length(ancs{alive(1)})-find(ancs{alive(1)}==mast);
                mlca(tt) = ml0;
                for ii = alive
                    ml0  = length(ancs{ii})-find(ancs{ii}==mast);
                    if ml0>mlca(tt)
                        mlca(tt) = ml0;
                    end
                end
            else
                mlca(tt) = mlca(tt);
                for ii = setdiff(alive,id0)
                    ml0  = length(ancs{ii})-find(ancs{ii}==mast);
                    if ml0>mlca(tt)
                        mlca(tt) = ml0;
                    end
                end
            end
            id0 = alive;
            dlca(tt) = (t-1)*h-birth(mast);
        else
            dlca(tt) = dlca(tt-1);
            mlca(tt) = mlca(tt-1);
        end
        %}
    end

end
talp(talp==0) = (t-1)*h;
net   = net+1;
et(net)= (t-1)*h-extmp;
toc

%% save data
%dirc  = './';
wname = 'Epidemics';
tname = 'WaitTime';
ename = 'Extinct';
gname = 'Genotype';
uname = sprintf('%.3f',100*mu0);
mname = sprintf('%.3f',m0*10000);
dname = sprintf('%.3f',d0);
nname = sprintf('%.2f',-log10(n0));
rname = sprintf('%03d',replica);
sname = sprintf('%.2f',(d0/b0));
dtype = '.dat';
epiname = [wname,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];
wtname  = [tname,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];
etname  = [ename,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];
gnname  = [gname,'_',uname,'_',mname,'_',dname,'_',sname,'_',nname,'_',rname,dtype];

dlmwrite(epiname,[h*hT:h*hT:T;totI;totR;ms;vs;sm;nstr;mb;vb;mf;vf;ts;tp;mlca;dlca;mdist;md;vd]);
dlmwrite(wtname,[wt;et_br;ngbr;naa]);
dlmwrite(etname,et);
dlmwrite(gnname,[prts;birth;death;beta;talp;droot;dx;xx;Imax;phi]);
