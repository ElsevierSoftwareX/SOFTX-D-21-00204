function [Freq,R1Th,R1_HM_best,R1_HH_best,tL,tB,tD,h,x_opt,...
N_opt,T1T2_best,Resid,R1b_HM_best, R1L_HM_best, R1_constantOffset_best, ...
R1bb_HH_Intra_best, ...
R1bb_HH_best,R1LL_HH_best,R1bL_HH_best,R1Lb_HH_best,...
R1b_HH_best,R1L_HH_best, xN_count_min_chi2,R1Xp] = Fit_3TM(R1Xp_in,...
FreqXp,ValXpFreq,ValHeaderSdf,...
HHTau,HMTau,IntPts,B_ind,L_ind,D_ind,Spin,R1_Offset,...
Npara,Npara_start,dNpara,Nx,x_start,dx,HH_triggerOn,...
HM_triggerOn,T1T2_freq,R1bb_Intra,chi_style,NbWdens,...
No_freq_start,No_freq_end,TauB_imp,TauD_imp,TauL_imp,Nv_imp,x_imp,Ns,T1T2Ratio)

%Set T1T2Ratio if default values 

            if isempty(T1T2Ratio)==1
                T1T2min=0
                T1T2max=100
                T1T2Ratio=[T1T2min,T1T2max]
            end
            
%Set experimental data arrays 
   
NFreq = size(FreqXp,1);       % the number of expt frequencies
Freq = zeros(NFreq,1);                     % experimental frequencies
R1Xp = zeros(NFreq,1);                     % experimental R1 data
No_freq_end=NFreq-No_freq_end;

% This sequence reverses order so that lowest frequency = index 1
    if ValXpFreq=="Descending"
        for i=1:NFreq                             
            Freq(i,1)=FreqXp(NFreq-i+1,1); 
            R1Xp(i,1)=R1Xp_in(NFreq-i+1,1); 
        end
    elseif ValHeaderSdf=="On"
        for i=1:NFreq                             
            Freq(i,1)=FreqXp(NFreq-i+1,1); 
            R1Xp(i,1)=R1Xp_in(NFreq-i+1,1); 
        end
    elseif ValXpFreq=="Ascending"
        for i=1:NFreq                             
            Freq(i,1)=FreqXp(i,1); 
            R1Xp(i,1)=R1Xp_in(i,1); 
        end
    end

    if HH_triggerOn==1
        if HM_triggerOn==0
            Npara=1
        else Npara=Npara
        end
    end
    
% Index of TauL
tau_L_rd = dlmread('tau_L.dat');
L_ind;

if isempty(L_ind)==1
    L_ind='[100000 0;10000 0]'
end  

strL = convertCharsToStrings(L_ind);
TestL = cellfun(@str2num,strL,'UniformOutput',false);
L=cell2mat(TestL);
L_min=min(L(L > 1));
L_min=round(L_min,2);
L_max=max(L(L > 1));
L_max=round(L_max,2);

if isempty(TauL_imp)==1
    TauL_imp=0
end

if TauL_imp>0
    L_imp=[]
    for i=1:113
            L_imp(end+1)=abs(TauL_imp-tau_L_rd(i,1));
    end
    L_i=min(L_imp)
    il_imp=find(L_imp==L_i)
        iL_start=il_imp
        iL_end=il_imp
else  
    Lind=[L_min,L_max];
    iL_start=(16*(log10(Lind(1,1))-1))+1;
    iL_end=(16*(log10(Lind(1,2))-1))+1;
    iL_start=round(iL_start);
    iL_end=round(iL_end); 
end

NL=iL_end - iL_start + 1;

% Index of TauB

if isempty(B_ind)==1
    B_ind='[100 0;4 0]'
end 

B_ind;
strB = convertCharsToStrings(B_ind);
TestB = cellfun(@str2num,strB,'UniformOutput',false);
B=cell2mat(TestB);
B_min=min(B(B > 1));
B_min=round(B_min);
B_max=max(B(B > 1));
B_max=round(B_max);
Bind=[B_min,B_max];

tau_b_rd = dlmread('tau_b.dat');

if isempty(TauB_imp)==1
    TauB_imp=0
end
    
if TauB_imp>0
    ib_end=22
    ib_start=1
else    
    B_start=[];
        for i=1:22
            B_start(end+1)=abs(Bind(1,1)-tau_b_rd(i,1));
        end
    B_s=min(B_start);
    ib_start = find(B_start==B_s);
    
    B_end=[];
        for i=1:22
            B_end(end+1)=abs(Bind(1,2)-tau_b_rd(i,1));
        end
    B_e=min(B_end);
    ib_end = find(B_end==B_e); 
end

Nb=ib_end - ib_start + 1; 

% Index of TauD
tau_d_in_rd = dlmread('tau_d.dat');
D_ind;

if isempty(D_ind)==1
    D_ind='[1 0;-1 0]'
end 

strD = convertCharsToStrings(D_ind);
TestD = cellfun(@str2num,strD,'UniformOutput',false);
D=cell2mat(TestD);
D_min=D(2,1);
D_max=D(1,1);

if isempty(TauD_imp)==1
    TauD_imp=0
end

if TauD_imp>0
    D_imp=[];
    for i=49*(il_imp-1)+1:49*(il_imp)+1
            D_imp(end+1)=abs(TauD_imp-tau_d_in_rd(i,1));
    end
    D_i=min(D_imp)
    iD_imp=find(D_imp==D_i)
    display(tau_d_in_rd(67*40,1))
    D=iD_imp
    id_start=iD_imp-25
    id_end=iD_imp-25
    
else
Dind=[D_min,D_max];
minD=min(Dind)
maxD=max(Dind)
% id_start=(16^(minD));
% id_end=(16^(maxD));
 id_start=(16*(minD));
 id_end=(16*(maxD));
id_start=round(id_start);
id_end=round(id_end);

end
Nd=abs(id_end - id_start + 1);

if Nd>48
    Nd=48
else Nd=Nd
end

% Input imposed regression values

if isempty(Nv_imp)==1
  Nv_imp=0  
end

if isempty(x_imp)==1
  x_imp=0
end

if Nv_imp>0
    Npara=1
    Npara_start=Nv_imp
end

if x_imp>0
    Nx=1
    x_start=x_imp
end 

% Input of Tau Numbers in the case of HM and HH interaction

    if isempty(HMTau)==1
        HMTau={'22','49','113'};
        HMTau = convertCharsToStrings(HMTau);
        HMTau = cellfun(@str2num,HMTau,'UniformOutput',false);
        HMTau=cell2mat(HMTau);
    else HMTau=HMTau
    end
    
    if isempty(HHTau)==1
        HHTau={'22','33','65'};
        HHTau = convertCharsToStrings(HHTau);
        HHTau = cellfun(@str2num,HHTau,'UniformOutput',false);
        HHTau=cell2mat(HHTau)
    else HHTau=HHTau
    end
   
    if isempty(IntPts)==1
        IntPts={'201'};
        IntPts = convertCharsToStrings(IntPts);
        IntPts = cellfun(@str2num,IntPts,'UniformOutput',false);
        IntPts=cell2mat(IntPts)
    else IntPts=IntPts
    end
    
    if isempty(Spin)==1
        Spin={'2.5'};
        Spin = convertCharsToStrings(Spin);
        Spin = cellfun(@str2num,Spin,'UniformOutput',false);
        Spin=cell2mat(Spin)
    else Spin=Spin
    end
    
    if isempty(T1T2_freq)==1
        T1T2_freq={'20'};
        T1T2_freq = convertCharsToStrings(T1T2_freq);
        T1T2_freq = cellfun(@str2num,T1T2_freq,'UniformOutput',false);
        T1T2_freq=cell2mat(T1T2_freq)
    else T1T2_freq=T1T2_freq
    end
    
    if isempty(NbWdens)==1
    NbWdens={'66.6'};
    NbWdens = convertCharsToStrings(NbWdens);
    NbWdens = cellfun(@str2num,NbWdens,'UniformOutput',false);
    NbWdens=cell2mat(NbWdens)
    else NbWdens=NbWdens
    end
 
Nd_half_range = (HMTau(1,2)+1)/2;

Nf=IntPts(1,1);

NLW=HHTau(1,3);NdW=HHTau(1,2); NbW=HHTau(1,1); NfW=IntPts;
NL_rd=HMTau(1,3) ;Nb_rd=HMTau(1,1) ;Nd_rd=HMTau(1,2);
NLW_rd=NLW; NdW_rd=NdW; NbW_rd=NbW; NfW_rd=NfW;

% Read in and organise the time arrays
 tau_L_rd = zeros(NL_rd,1);
 tau_d_rd = zeros(NL_rd,Nd_rd);          tau_d_in_rd = zeros(NL_rd*Nd_rd,1);
 tau_b_rd = zeros(Nb_rd,1);
 
 tau_L_rd = dlmread('tau_L.dat');
 tau_d_in_rd = dlmread('tau_d.dat');
 tau_b_rd = dlmread('tau_b.dat');


 tau_L = zeros(NL,1);
 tau_d = zeros(NL,Nd);          tau_d_in = zeros(NL*Nd,1);
 tau_b = zeros(Nb,1);
 
 % tau_L
    for iL = 1:NL
       tau_L(iL,1) = tau_L_rd(iL_start-1+iL,1);
%         tau_L(iL,1) = tau_L_rd(iL_start+iL,1);
    end
 % tau_d
    tcount = 0;
    for iL=1:NL_rd    % reorganise tau_d read input into array
        for id=1:Nd_rd
            tcount = tcount+1;
            tau_d_in(iL,id) = tau_d_in_rd(tcount);
        end
    end

    for iL=1:NL                    % reorganise tau_d read input into array
        for id=1:Nd
             tau_d(iL,id) = tau_d_in(iL_start+iL,id_start+Nd_half_range+id-1);
        end
    end

    for ib = 1:Nb
       tau_b(ib,1) = tau_b_rd(ib_start-1+ib,1);
    end
%--------------------------------------------------------------------------
% Arrays for the R1(w) and R2(w) pre-calculated data sets HM HM
 R1_L_rd = zeros(NL*Nd*Nf,1);        R2_L_rd = zeros(NL*Nd*Nf,1);
 R1_b_rd = zeros(Nb*Nf,1);           R2_b_rd = zeros(Nb*Nf,1);
 
 R1_L = zeros(NL,Nd,Nb,Nf);         R2_L = zeros(NL,Nd,Nb,Nf);     
 R1_b = zeros(NL,Nd,Nb,Nf);         R2_b = zeros(NL,Nd,Nb,Nf);
 
 R1_L_T1T2 = zeros(NL,Nd,Nb);       R2_L_T1T2 = zeros(NL,Nd,Nb);
 R1_b_T1T2 = zeros(NL,Nd,Nb);       R2_b_T1T2 = zeros(NL,Nd,Nb);
 
 R1_constantOffset_X = zeros(NL,Nd,Nb,NFreq);
 R1_constantOffset_X_T1T2 = zeros(NL,Nd,Nb);
 
 Freq_th = zeros(Nf,1);         lg_Freq_th = zeros(Nf,1);
 
 % read HM data
 fprintf('\n Reading HM input relaxation rate data sets ... \n');
 
 Freq_th=dlmread('Freq_th.dat');   % all frequencies for theory R1

 R1_L_rd = dlmread('R1_L.dat');
 R1_b_rd = dlmread('R1_b.dat');
 R2_L_rd = dlmread('R2_L.dat');
 R2_b_rd = dlmread('R2_b.dat');
 
 % pre-calculated data sets do not include this spin factor
 R1_L_rd = R1_L_rd*Spin*(Spin+1.0);       
 R1_b_rd = R1_b_rd*Spin*(Spin+1.0);
 R2_L_rd = R2_L_rd*Spin*(Spin+1.0);
 R2_b_rd = R2_b_rd*Spin*(Spin+1.0);
 
 R1_L_ind = zeros(NL,Nd,Nb,NFreq);      R1_b_ind = zeros(NL,Nd,Nb,NFreq);
 
 fprintf('... input relaxation rate data sets read successfully \n\n');
 % end HM end HM end HM end HM end HM end HM end HM end HM end HM 
 
 
 % WATER WATER WATER WATER WATER WATER WATER WATER WATER WATER 
 % Arrays for the R1(w) and R2(w) pre-calculated data sets HH WATER
 
 %%%  set the arrays ==============================================
 % These must be set to the same array dimensions as for HM
 R1W_LL = zeros(NL,Nd,Nb,Nf);   R1W_bb = zeros(NL,Nd,Nb,Nf);
 R1W_Lb = zeros(NL,Nd,Nb,Nf);
 R2W_LL = zeros(NL,Nd,Nb,Nf);   R2W_bb = zeros(NL,Nd,Nb,Nf);
 R2W_Lb = zeros(NL,Nd,Nb,Nf);
  
 % these arrays for reading data from file
 % These must be set to the WATER array dimensions
 R1W_LL_rd = zeros(NLW*NdW*NfW,1);    R2W_LL_rd = zeros(NLW*NdW*NfW,1);
 R1W_Lb_rd = zeros(NLW*NdW*NbW*NfW,1);R2W_Lb_rd = zeros(NLW*NdW*NbW*NfW,1);
 R1W_bb_rd = zeros(NbW*NfW,1);        R2W_bb_rd = zeros(NbW*NfW,1);
 
 % These must be set to the same array dimensions as for HM
 R1W_L = zeros(NL,Nd,Nb,Nf);     
 R2W_L = zeros(NL,Nd,Nb,Nf);     
 R1W_b = zeros(NL,Nd,Nb,Nf);        
 R2W_b = zeros(NL,Nd,Nb,Nf);
 
 R1W_LL_T1T2 = zeros(NL,Nd,Nb);  R2W_LL_T1T2 = zeros(NL,Nd,Nb);
 R1W_bb_T1T2 = zeros(NL,Nd,Nb);  R2W_bb_T1T2 = zeros(NL,Nd,Nb);
 R1W_Lb_T1T2 = zeros(NL,Nd,Nb);  R2W_Lb_T1T2 = zeros(NL,Nd,Nb);
 
 % Read in R1W(w) and R2W(w) pre-calculated WATER data sets
 fprintf('\n Reading HH input relaxation rate data sets ... \n');
 
 R1W_LL_rd = dlmread('R1_LL_water.dat');
 R1W_bb_rd = dlmread('R1_bb_water.dat');
 R1W_Lb_rd = dlmread('R1_Lb_water.dat');
 R2W_LL_rd = dlmread('R2_LL_water.dat');
 R2W_bb_rd = dlmread('R2_bb_water.dat');
 R2W_Lb_rd = dlmread('R2_Lb_water.dat');
 
 R1W_LL_ind = zeros(NL,Nd,Nb,NFreq);    R1W_bb_ind = zeros(NL,Nd,Nb,NFreq);
 R1W_Lb_ind = zeros(NL,Nd,Nb,NFreq);      
  
 fprintf('... input relaxation rate data sets read successfully \n\n');
 % end WATER end WATER end WATER end WATER end WATER end WATER end WATER  
 % =======================================================================
 
 R1_constantOffset_ind = zeros(NL,Nd,Nb,NFreq); 
 
 R1_constantOffset = zeros(NFreq,1); 
 R1_constantOffset = R1_Offset;

 %%%  set the arrays =====================================================

 R1 = zeros(NL,Nd,Nb,NFreq); 
 
 R1b_HM = zeros(NL,Nd,Nb,NFreq);        R1L_HM = zeros(NL,Nd,Nb,NFreq);
 
 R1bb_HH = zeros(NL,Nd,Nb,NFreq);       R1LL_HH = zeros(NL,Nd,Nb,NFreq);
 R1Lb_HH = zeros(NL,Nd,Nb,NFreq);       R1bL_HH = zeros(NL,Nd,Nb,NFreq);
 R1b_HH = zeros(NL,Nd,Nb,NFreq);        R1L_HH = zeros(NL,Nd,Nb,NFreq);
 
 R1b_HM_best = zeros(NFreq,1);          R1L_HM_best = zeros(NFreq,1); 
 R1_HM_best = zeros(NFreq,1);
 
 R1_constantOffset_best = zeros(NFreq,1);
 R1bb_HH_Intra_best = zeros(NFreq,1);

 R1bb_HH_best = zeros(NFreq,1);         R1LL_HH_best = zeros(NFreq,1);
 R1bL_HH_best = zeros(NFreq,1);         R1Lb_HH_best = zeros(NFreq,1);
 R1b_HH_best = zeros(NFreq,1);          R1L_HH_best = zeros(NFreq,1);
 R1_HH_best = zeros(NFreq,1);
 
 R1_off = zeros(NFreq,1); 
 R1Xp_ind = zeros(NL,Nd,Nb,NFreq);
 R1Th = zeros(NFreq,1); 
 R1_best = zeros(NFreq,1);       
 Resid = zeros(NFreq,1);                % residuals       
 
 chi2_f = zeros(NL,Nd,Nb,NFreq);        
 chi2 = zeros(NL,Nd,Nb);                      

 Freq_ind = zeros(NFreq,1);     lg_Freq = zeros(NFreq,1);
 
 NxN = Npara*Nx;
 xx = zeros(NxN,1);
 NN = zeros(NxN,1);
 min_chi2 = zeros(NxN,1);
 
% Set the trial s/v ratio values for fit
 x = zeros(Nx,1);                 
 for i=1:Nx; x(i,1)=x_start + i*dx; end       

% Set the trial Nlayer values for fit
 N = zeros(Npara,1);      
 for i=1:Npara; N(i,1)=Npara_start+((i-1)*dNpara); end
 
 % Find the theory frequency indexes that best correspond to expt f
 for i=1:NFreq
      
    delta_f(:,1) = abs(log10(Freq_th(:,1)) - log10(Freq(i,1)));
    minimum = min(delta_f);
  
    Freq_ind(i,1) = find(delta_f == minimum);
    
 end
 
 lg_Freq = log10(Freq);
 lg_Freq_th = log10(Freq_th);
  
% Here the index in the list of theory frequencies that corresponds to T1T2_freq is found
  delta_T1T2(:,1) = abs(Freq_th(:,1) - T1T2_freq);
  minimum = min(delta_T1T2);
  T1T2_index = find(delta_T1T2 == minimum);
  
  % reorganise read input into 4D arrays: HM rates
RcountL = 0;
for iL=1:NL                             
    for id=1:Nd
        idn = id_start+Nd_half_range+id;
        iLn = iL_start+iL-1;
        RcountL = (iLn-1)*Nd_rd + idn;
        
        for ib=1:Nb
            
            ibn = ib_start-1+ib;
            R1Xp_ind(iL,id,ib,1:NFreq) = R1Xp(1:NFreq,1);
            
            for ifr=1:Nf
            Rindexb = (ibn-1)*Nf + ifr;
            RindexL = (RcountL-1)*Nf + ifr;

            R1_L(iL,id,ib,ifr) = R1_L_rd(RindexL,1);
            R1_b(iL,id,ib,ifr) = R1_b_rd(Rindexb,1);
            R2_L(iL,id,ib,ifr) = R2_L_rd(RindexL,1);
            R2_b(iL,id,ib,ifr) = R2_b_rd(Rindexb,1);
            end
            
            R1_L_T1T2(iL,id,ib) = R1_L(iL,id,ib,T1T2_index);
            R1_b_T1T2(iL,id,ib) = R1_b(iL,id,ib,T1T2_index);
            R2_L_T1T2(iL,id,ib) = R2_L(iL,id,ib,T1T2_index);
            R2_b_T1T2(iL,id,ib) = R2_b(iL,id,ib,T1T2_index);

        end
    end
end

% reorganise read input into 4D arrays: HH WATER rates
RcountL = 0;
Nd_W2HM_shift = ((Nd_rd-NdW_rd)/2)+1;
NL_W2HM_shift = NL_rd-NLW_rd+1;

for iL=1:NL                             
    for id=1:Nd
        
        idn = id_start+Nd_half_range+id;
        iLn = iL_start-NL_W2HM_shift+iL;
        RcountL = (iLn-1)*NdW_rd + idn;
              
        for ib=1:Nb
            
            ibn = ib_start-1+ib;
            
            for ifr=1:Nf
            Rindexb = (ibn-1)*Nf + ifr;
            RindexL = (RcountL-1)*Nf + ifr;

            R1W_LL(iL,id,ib,ifr) = R1W_LL_rd(RindexL,1);
            R1W_bb(iL,id,ib,ifr) = R1W_bb_rd(Rindexb,1);
            R1W_Lb(iL,id,ib,ifr) = R1W_Lb_rd(Rindexb,1);
            R2W_LL(iL,id,ib,ifr) = R2W_LL_rd(RindexL,1);
            R2W_bb(iL,id,ib,ifr) = R2W_bb_rd(Rindexb,1);
            R2W_Lb(iL,id,ib,ifr) = R2W_Lb_rd(Rindexb,1);
            end
            
            R1W_LL_T1T2(iL,id,ib) = R1W_LL(iL,id,ib,T1T2_index);
            R1W_bb_T1T2(iL,id,ib) = R1W_bb(iL,id,ib,T1T2_index);
            R1W_Lb_T1T2(iL,id,ib) = R1W_Lb(iL,id,ib,T1T2_index);
            R2W_LL_T1T2(iL,id,ib) = R2W_LL(iL,id,ib,T1T2_index);
            R2W_bb_T1T2(iL,id,ib) = R2W_bb(iL,id,ib,T1T2_index);
            R2W_Lb_T1T2(iL,id,ib) = R2W_Lb(iL,id,ib,T1T2_index);

        end
    end
end

% Here we match the theory log10 frequencies to expt frequencies for fitting
% and use quadratic scaling to obtain the best R1_L and R1_b at the 
% expt frequency interpolated between the theory log10 frequencies

%Index of TauB_imp
    if TauB_imp==0 
        B=1:Nb
    else
        B_imp=[];
        for i=1:22
                B_imp(end+1)=abs(TauB_imp-tau_b_rd(i,1));
        end
        B_i=min(B_imp)
        ib_imp=find(B_imp==B_i)
        B=ib_imp
    end
    
% Index of TauL_imp
if TauL_imp==0
    L=1:NL
else L_imp=[];
    for i=1:113
            L_imp(end+1)=abs(TauL_imp-tau_L_rd(i,1));
    end
    L_i=min(L_imp)
    il_imp=find(L_imp==L_i)
        L=il_imp-iL_start+1
end

% Index of TauD_imp
    if TauD_imp==0 
        D=1:Nd
    else D_imp=[];
    for i=49*(il_imp-1)+1:49*(il_imp)+1
            D_imp(end+1)=abs(TauD_imp-tau_d_in_rd(i,1));
    end
    D_i=min(D_imp)
    iD_imp=find(D_imp==D_i)
    display(tau_d_in_rd(67*40,1))
    D=iD_imp-id_start-24
    end

d_log_freq_th = lg_Freq_th(2,1) - lg_Freq_th(1,1);

% Loop for creating the fitting arrays

for iL=L                 
    for id=D
        for ib=B  
            for i=1:NFreq
              
              ind = Freq_ind(i,1);
              
              xlgfXp = (lg_Freq(i,1) - lg_Freq_th(ind-1,1))/d_log_freq_th;
              factor0 = (xlgfXp-1)*(xlgfXp-2)/2;
              factor1 = -xlgfXp*(xlgfXp-2);
              factor2 = xlgfXp*(xlgfXp-1)/2;
              
              R1_b_ind(iL,id,ib,i) =  factor0*R1_b(iL,id,ib,ind-1) + ...
                                      factor1*R1_b(iL,id,ib,ind) + ...
                                      factor2*R1_b(iL,id,ib,ind+1);
                
              R1_L_ind(iL,id,ib,i) =  factor0*R1_L(iL,id,ib,ind-1) + ...
                                      factor1*R1_L(iL,id,ib,ind) + ...
                                      factor2*R1_L(iL,id,ib,ind+1);
                                  
              R1_constantOffset_ind(iL,id,ib,:) = R1_constantOffset(:,1);
                                  
              R1W_bb_ind(iL,id,ib,i) =  factor0*R1W_bb(iL,id,ib,ind-1) + ...
                                        factor1*R1W_bb(iL,id,ib,ind) + ...
                                        factor2*R1W_bb(iL,id,ib,ind+1);
                
              R1W_LL_ind(iL,id,ib,i) =  factor0*R1W_LL(iL,id,ib,ind-1) + ...
                                        factor1*R1W_LL(iL,id,ib,ind) + ...
                                        factor2*R1W_LL(iL,id,ib,ind+1);
                                    
              R1W_Lb_ind(iL,id,ib,i) =  factor0*R1W_Lb(iL,id,ib,ind-1) + ...
                                        factor1*R1W_Lb(iL,id,ib,ind) + ...
                                        factor2*R1W_Lb(iL,id,ib,ind+1);
                                    
            end    
        end
    end
end

fprintf('... fitting arrays complete \n\n');
fprintf('Starting fit ... \n');

%%% Here the fitting starts - find best x and N for each tau

% The sole purpose of the initial scoping is to locate the global min chi^2
% fit parameter for each (x,N) pairing for each tau_L, tau_d and tau_b  


xN_count = 0;
for ix=1:Nx                 % x loop
    
    for iN=1:Npara          % N loop
        
    xN_count = xN_count + 1;
    
    xx(xN_count,1) = x(ix,1);
    NN(xN_count,1) = N(iN,1); 
    
    % Calculate R1
    R1b_HM = HM_triggerOn*( (1-x(ix,1))*N(iN,1)*R1_b_ind );
    R1L_HM = HM_triggerOn*( N(iN,1)*x(ix,1)*R1_L_ind );
    R1_HM = R1b_HM + R1L_HM;
    
    R1_constantOffset_X = R1_constantOffset_ind;
    
    R1bb_HH_Intra = R1bb_Intra;
    
    R1bb_HH = HH_triggerOn*( (1-x(ix,1))*NbWdens*R1W_bb_ind + R1bb_HH_Intra);
    R1LL_HH = HH_triggerOn*( x(ix,1)*Ns*R1W_LL_ind );    
    R1bL_HH = HH_triggerOn*( (1-x(ix,1))*NbWdens*R1W_Lb_ind );
    R1Lb_HH = HH_triggerOn*( x(ix,1)*Ns*R1W_Lb_ind ); 
    
    R1b_HH = R1bb_HH + R1bL_HH;
    R1L_HH = R1LL_HH + R1Lb_HH;
    R1_HH = R1b_HH + R1L_HH;
    
    R1 = R1_HM + R1_constantOffset_X + R1_HH;
    
                R1b_HM_T1T2 = HM_triggerOn*( (1-x(ix,1))*N(iN,1)*R1_b_T1T2 );
                R1L_HM_T1T2 = HM_triggerOn*( N(iN,1)*x(ix,1)*R1_L_T1T2 );
                R1_HM_T1T2 = R1b_HM_T1T2 + R1L_HM_T1T2;
                
                R1bb_HH_T1T2 = HH_triggerOn*( (1-x(ix,1))*NbWdens*R1W_bb_T1T2 + R1bb_HH_Intra);
                R1LL_HH_T1T2 = HH_triggerOn*( x(ix,1)*Ns*R1W_LL_T1T2 );    
                R1bL_HH_T1T2 = HH_triggerOn*( (1-x(ix,1))*NbWdens*R1W_Lb_T1T2 );
                R1Lb_HH_T1T2 = HH_triggerOn*( x(ix,1)*Ns*R1W_Lb_T1T2 );
                
                R1b_HH_T1T2 = R1bb_HH_T1T2 + R1bL_HH_T1T2;
                R1L_HH_T1T2 = R1LL_HH_T1T2 + R1Lb_HH_T1T2;
                R1_HH_T1T2 = R1b_HH_T1T2 + R1L_HH_T1T2;
    
                R1_T1T2 = R1_HM_T1T2 + R1_constantOffset_X_T1T2 + R1_HH_T1T2;
                
                R2b_HM = HM_triggerOn*( (1-x(ix,1))*N(iN,1)*R2_b_T1T2 );
                R2L_HM = HM_triggerOn*( N(iN,1)*x(ix,1)*R2_L_T1T2 );
                R2_HM = R2b_HM + R2L_HM;

                R2bb_HH = HH_triggerOn*( (1-x(ix,1)*NbWdens*R2W_bb_T1T2 + R1bb_HH_Intra));
                R2LL_HH = HH_triggerOn*( x(ix,1)*Ns*R2W_LL_T1T2 );    
                R2bL_HH = HH_triggerOn*( (1-x(ix,1)*NbWdens*R2W_Lb_T1T2 ));
                R2Lb_HH = HH_triggerOn*( x(ix,1)*Ns*R2W_Lb_T1T2 ); 

                R2b_HH = R2bb_HH + R2bL_HH;
                R2L_HH = R2LL_HH + R2Lb_HH;
                R2_HH  = R2b_HH  + R2L_HH;

                R2_T1T2 =  R2_HM + R1_constantOffset_X_T1T2 ...
                                + R2_HH;

                T1_T2 = R2_T1T2./R1_T1T2;
       
    if strcmp(chi_style,'log')==1
        
        for i=No_freq_start:(NFreq-No_freq_end)
        chi2_f(:,:,:,i) = (log10(R1(:,:,:,i)) - log10(R1Xp_ind(:,:,:,i))).^2;          % log chi^2
        end
        
    elseif strcmp(chi_style,'lin')==1
        
        for i=No_freq_start:(NFreq-No_freq_end)
        chi2_f(:,:,:,i) = (R1(:,:,:,i) - R1Xp_ind(:,:,:,i)).^2;                     % lin chi^2
        end
        
     elseif strcmp(chi_style,'wtlin')==1
        
        for i=No_freq_start:(NFreq-No_freq_end)
        chi2_f(:,:,:,i) = ((R1(:,:,:,i) - R1Xp_ind(:,:,:,i)).^2)/R1Xp_err(i,1)^2;   % weighted lin chi^2                   % lin chi^2
        end
        
    end


    chi2 = sum(chi2_f,4);                               % sum chi^2 for each expt frequency
   
        for iL=1:numel(L)                 
      for id=1:numel(D)
        for ib=1:numel(B) 

         if T1_T2(iL,id,ib)<T1T2Ratio(1,1) || T1T2Ratio(1,2)<T1_T2(iL,id,ib)
         chi2(iL,id,ib)=1E8;  
             end
         end
      end
    end
    
   min_chi2(xN_count,1) = min(min(min(chi2)));      % global minimum chi^2 value for this x and Nlayer
    end
end

% This is the minimum chi^2 from xN_count list to find x and N

xN_count_min_chi2 = min(min_chi2(min_chi2>0));          % this is the minimum chi2
% xN_count_min_chi2 = min(min_chi2(min_chi2>1));          % this is the minimum chi2
i_xn_opt = find(min_chi2 == xN_count_min_chi2);         % find the index this refers to
x_opt = xx(i_xn_opt,1);                                 % This is the corresponding x at global min
N_opt = NN(i_xn_opt,1);                                 % This is the corresponding Nlayer at global min

% recalculate R1 etc for best choice of N and x for all times
    R1b_HM = HM_triggerOn*( (1-x_opt)*N_opt*R1_b_ind );
    R1L_HM = HM_triggerOn*( N_opt*x_opt*R1_L_ind );
    R1_HM = R1b_HM + R1L_HM;
    
    R1bb_HH = HH_triggerOn*( (1-x_opt)*NbWdens*R1W_bb_ind + R1bb_HH_Intra);
    R1LL_HH = HH_triggerOn*( x_opt*Ns*R1W_LL_ind );    
    R1bL_HH = HH_triggerOn*( (1-x_opt)*NbWdens*R1W_Lb_ind );
    R1Lb_HH = HH_triggerOn*( x_opt*Ns*R1W_Lb_ind ); 
    
    R1b_HH = R1bb_HH + R1bL_HH;
    R1L_HH = R1LL_HH + R1Lb_HH;
    R1_HH = R1b_HH + R1L_HH;
    
    R1 = R1_HM + R1_constantOffset_X + R1_HH;

%-----------------------------------------
    % Now recalculate chi_2
    
     if strcmp(chi_style,'log')==1
        for i=No_freq_start:(NFreq-No_freq_end)
        chi2_f(:,:,:,i) = (log10(R1(:,:,:,i)) - log10(R1Xp_ind(:,:,:,i))).^2;          % log chi^2
        end
     elseif strcmp(chi_style,'lin')==1
        for i=No_freq_start:(NFreq-No_freq_end)
        chi2_f(:,:,:,i) = (R1(:,:,:,i) - R1Xp_ind(:,:,:,i)).^2;                     % lin chi^2
        end
     elseif strcmp(chi_style,'wtlin')==1
        for i=No_freq_start:(NFreq-No_freq_end)
        chi2_f(:,:,:,i) = ((R1(:,:,:,i) - R1Xp_ind(:,:,:,i)).^2)/R1Xp_err(i,1)^2;   % weighted lin chi^2                   % lin chi^2
        end
     end
    
    chi2 = sum(chi2_f,4);                               % sum chi^2 for each expt frequency
    
    i_opt = find(chi2 == xN_count_min_chi2);          % 

% find time indexes for best fit
    N1 = (i_opt - mod(i_opt-1,Nd*NL) - 1)/(NL*Nd);
    tb_best_index = N1+1;
    N2 = i_opt - N1*Nd*NL;
    N3 = (N2 - mod(N2-1,NL) - 1)/NL;
    td_best_index = N3+1;
    tL_best_index = N2 - N3*NL;
    
    % recalculate R1 etc for best choice of N, x and time constants
    %-----------------------------------------
    R1_best(:,1) = R1(tL_best_index,td_best_index,tb_best_index,:);
    R1Th(:,1) = R1_best(:,1);

    R1b_HM_best(:,1) = R1b_HM(tL_best_index,td_best_index,tb_best_index,:);
    R1L_HM_best(:,1) = R1L_HM(tL_best_index,td_best_index,tb_best_index,:);
    R1_HM_best(:,1) = R1_HM(tL_best_index,td_best_index,tb_best_index,:);

    R1_constantOffset_best(:,1) = R1_constantOffset_X(tL_best_index,td_best_index,tb_best_index,:);

    R1bb_HH_Intra_best(:,1) = R1bb_HH_Intra;

    R1bb_HH_best(:,1) = R1bb_HH(tL_best_index,td_best_index,tb_best_index,:);
    R1LL_HH_best(:,1) = R1LL_HH(tL_best_index,td_best_index,tb_best_index,:);
    R1bL_HH_best(:,1) = R1bL_HH(tL_best_index,td_best_index,tb_best_index,:);
    R1Lb_HH_best(:,1) = R1Lb_HH(tL_best_index,td_best_index,tb_best_index,:);
    R1b_HH_best(:,1) = R1b_HH(tL_best_index,td_best_index,tb_best_index,:);
    R1L_HH_best(:,1) = R1L_HH(tL_best_index,td_best_index,tb_best_index,:);
    R1_HH_best(:,1) = R1_HH(tL_best_index,td_best_index,tb_best_index,:);
                   
% recalculate R1 for this choice N_opt and x_opt at T1T2 frequency

    iL=tL_best_index;
    id=td_best_index;
    ib=tb_best_index;   
    
    R1b_HM_T1T2 = HM_triggerOn*( (1-x_opt)*N_opt*R1_b_T1T2 );
    R1L_HM_T1T2 = HM_triggerOn*( N_opt*x_opt*R1_L_T1T2 );
    R1_HM_T1T2 = R1b_HM_T1T2 + R1L_HM_T1T2;
                
    R1bb_HH_T1T2 = HH_triggerOn*( (1-x_opt)*NbWdens*R1W_bb_T1T2 + R1bb_HH_Intra);
    R1LL_HH_T1T2 = HH_triggerOn*( x_opt*Ns*R1W_LL_T1T2 );    
    R1bL_HH_T1T2 = HH_triggerOn*( (1-x_opt)*NbWdens*R1W_Lb_T1T2 );
    R1Lb_HH_T1T2 = HH_triggerOn*( x_opt*Ns*R1W_Lb_T1T2 ); 
                
    R1b_HH_T1T2 = R1bb_HH_T1T2 + R1bL_HH_T1T2;
    R1L_HH_T1T2 = R1LL_HH_T1T2 + R1Lb_HH_T1T2;
    R1_HH_T1T2 = R1b_HH_T1T2 + R1L_HH_T1T2;
                
    R1T1T2_best = R1_HM_T1T2(iL,id,ib) + R1_constantOffset_X_T1T2(iL,id,ib) ...
                    + R1_HH_T1T2(iL,id,ib);
                                     
   % Calculate R2 for this choice N_opt and x_opt at T1T2 frequency
    R2b_HM_T1T2 = HM_triggerOn*( (1-x_opt)*N_opt*R2_b_T1T2 );
    R2L_HM_T1T2 = HM_triggerOn*( N_opt*x_opt*R2_L_T1T2 );
    R2_HM_T1T2 = R2b_HM_T1T2 + R2L_HM_T1T2;
                
    R2bb_HH_T1T2 = HH_triggerOn*( (1-x_opt)*NbWdens*R2W_bb_T1T2 + R1bb_HH_Intra);
    R2LL_HH_T1T2 = HH_triggerOn*( x_opt*Ns*R2W_LL_T1T2 );    
    R2bL_HH_T1T2 = HH_triggerOn*( (1-x_opt)*NbWdens*R2W_Lb_T1T2 );
    R2Lb_HH_T1T2 = HH_triggerOn*( x_opt*Ns*R2W_Lb_T1T2 ); 
                
    R2b_HH_T1T2 = R2bb_HH_T1T2 + R2bL_HH_T1T2;
    R2L_HH_T1T2 = R2LL_HH_T1T2 + R2Lb_HH_T1T2;
    R2_HH_T1T2  = R2b_HH_T1T2  + R2L_HH_T1T2;
                
    R2T1T2_best =  R2_HM_T1T2(iL,id,ib) + R1_constantOffset_X_T1T2(iL,id,ib) ...
                    + R2_HH_T1T2(iL,id,ib);
                
    T1T2_best = R2T1T2_best/R1T1T2_best;
    
h = 0.54/x_opt/1000;    % "planar pore equivalent" thickness (0.54 nm surface converted to um)


tL = tau_L(tL_best_index)/1000000;      % in um
tD = tau_d(tL_best_index,td_best_index)/1000000;      % in um

if TauB_imp>0
    tB = tau_b_rd(tb_best_index);      % in ps 
else
    tB = tau_b(tb_best_index);     % in ps 
end

Resid(:,1) = R1Th(:,1)-R1Xp(:,1);

fprintf('... fitting complete \n\n');
fprintf('Writing output ... \n');

fprintf('\n FIT COMPLETE!!\n');


end