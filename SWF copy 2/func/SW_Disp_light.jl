"SW_Disp: Surface Wave dispersion for Julia
This file is a small package to perform forward modelling of Surface Wave Dispersion.
It is based on the codes by Haney and Tsai (2015,2017,2019) as well as on many comments
Prof. Robert Herrmann and his Computer Programs in Seismology (version 3.30).
Use with care.
This package was developed by Mariano Arnaiz (marianoarnaiz@gmail.com)
at the Universidad Complutense de Madrid
in spring 2021 as part of the WINTER_C project with Dr. Javier Fullea."

"NOTE TO USER: THIS IS THE LIGHT VERSION OF THE CODE, SUITABLE FOR INVERSION BUT NOT FOR TESTING OF FORWARD MODELLING"
## Begin by loading the modules used
# All this are necessary to run SW_DISP properly
# You need to istall "brew install arpack"!!!!!
using PolynomialRoots, LinearAlgebra, Arpack, SparseArrays, Interpolations

## CALLING FUNCTION
"This function only CALLS ALL THE OTHER FUNCTIONS. It is intended to ease the input/output of the data"
##
function SW_Disp_light(Model::Array{Float64,2},T0::Int64,Tstep::Int64,TF::Int64,Nmodes::Int64)

# Print Message for user
# Prepare Inputs
@fastmath @inbounds fks,Mesh,vpv,vsv,vsv_v,vsv_h,rhovR,rhovL,qsv,qR,Nsolid,hsolid,Nfluid,hfluid,vpfluid,vpfv,rhofv,T,fnum = fix_SWF_input(T0,Tstep,TF,Model);
# Compute de Dispersion for R Wave
@fastmath @inbounds cR, UR= Rayleigh_Forwardsp(Nsolid,vsv_v,vpv,rhovR,fks,fnum,hsolid,Nmodes,Nfluid,vpfv,rhofv,hfluid,qsv);
# Compute de Dispersion for L Wave
@fastmath @inbounds cL, UL = Love_Forwardsp(Nsolid,vsv_h,rhovL,fks,fnum,hsolid,Nmodes,qsv);
# Organize the Outputs
@fastmath @inbounds SWV= fix_SWF_output(cR, UR, cL, UL, fks);
println("SWF_Disp_light > T(s) SWV[:, 1(cR) 2(UR) 3(cL) 4(UL), Mode]")
return T, SWV
end

## FUNTION FOR INPUT Data
"This function is intended to fix the inputs for both of the codes that require a bunch of inputs"
##
function fix_SWF_input(T0::Int64,Tstep::Int64,TF::Int64,Model::Array{Float64,2})
T=collect(TF:-Tstep:T0);
fks= @. 1/T; # vector of  frequencies  at which the velocities are measured (Hz) Change T (s) range
Step=1; #step for the fluid
#Depths
Mesh=Model[2:end-1,1];

# Properties
vpv=Model[2:end-1,2];
vsv=Model[2:end-1,3];
rhov=Model[2:end-1,4];

# Atenuaction Qs
qsv=Model[2:end-1,5];

# Radial Anysotropy
qR=Model[2:end-1,6]*0.01;
vsv_v=@. ((3-qR)*0.3333333333333333)*vsv; #Vertical velocity in realtion to Radial Anisitropy
vsv_h=@. (1+0.6666666666666666*qR)*vsv; #Horizontal velocity equal to total velocity


# FEM MODEL
# The model begins with the number of elements
# SOLID PART (the code can consider the water layer)
Nsolid = size(vpv,1); # number of elements in solid
hsolid =[diff(Mesh); Model[end,1]-Mesh[end,1]];
# Fluid Part
Nfluid = 0; # number of elements in fluid
hfluid = Step*ones(1,Nfluid); # grid spacing of mesh (meters)
vpfluid = 1.5; rhofluid = 1.03;
vpfv = vpfluid*ones(1,Nfluid);
rhofv = rhofluid*ones(1,Nfluid);

# Spherical Correction:
ER=6371;
HH=[0; hsolid]
RAD=ER .-[cumsum(HH[1:end-1]) cumsum(HH[2:end])]
CORRV=@. (2*ER)/(RAD[:,1]+RAD[:,2]);

#Apply to R wave data
CORRRHO_R=@. CORRV^-2.275;
vpv=@. vpv*CORRV
vsv_v=@. vsv_v*CORRV
rhovR=@. rhov*CORRRHO_R

#Apply to L wave data
CORRRHO_L=@. CORRV^-5;
vsv_h=@. vsv_h*CORRV
rhovL=@. rhov*CORRRHO_L

# frecuency
fnum=size(fks,1);

# other useful Components

return fks,Mesh,vpv,vsv,vsv_v,vsv_h,rhovR,rhovL,qsv,qR,Nsolid,hsolid,Nfluid,hfluid,vpfluid,vpfv,rhofv,T,fnum
end

## FUNCTION TO FIX Outputs
" This function cleans a bit the outut of the 2 main functions"
##
function fix_SWF_output(cR::Array{Float64,3}, UR::Array{Float64,3},  cL::Array{Float64,3}, UL::Array{Float64,3}, fks::Array{Float64,1})
# Organize the Velocity Dispersion Output
SWV=zeros(size(fks,1),4,Nmodes);
SWV[:,1,1:Nmodes]=cR[1:Nmodes,1,:]'#[cR[1,1,:] cR[2,1,:] cR[3,1,:]];  #safe Rayleigh phase velocity (m/s)
SWV[:,2,1:Nmodes]=UR[1:Nmodes,1,:]'#[UR[1,1,:] UR[2,1,:] UR[3,1,:]];  #safe Rayleigh group veloecity (m/s)
SWV[:,3,1:Nmodes]=cL[1:Nmodes,1,:]'#[cL[1,1,:] cL[2,1,:] cL[3,1,:]];  #safe Love phase velocity (m/s)
SWV[:,4,1:Nmodes]=UL[1:Nmodes,1,:]'#[UL[1,1,:] UL[2,1,:] UL[3,1,:]];  #safe Love group veloecity (m/s)
# Get R wave Vertical Eigenfunction

return SWV
end

## Compute Stoneley waves velocities (Explenation at the end)
" Compute Stoneley waves velocities"
##
function Stoneley(a,b,c,f,s)
# Compute all the polynomial coefficients
c16 = 256*(b^16) - (512*(b^18))/(a^2) + (256*(b^20))/(a^4);
c14 = -768*(b^14) + (1280*(b^16))/(a^2) - (512*(b^18))/(a^4) - (512*(b^16))/(c^2) + (1024*(b^18))/((a^2)*(c^2)) - (512*(b^20))/((a^4)*(c^2));
c12 = 832*(b^12) - (1024*(b^14))/(a^2) + (256*(b^16))/(a^4) + (256*(b^16))/(c^4) - (512*(b^18))/((a^2)*(c^4)) + (256*(b^20))/((a^4)*(c^4)) + (1536*(b^14))/(c^2) - (2560*(b^16))/((a^2)*(c^2)) + (1024*(b^18))/((a^4)*(c^2)) - (64*(b^12)*(f^2))/(s^2);
c10 = -416*(b^10) + (288*(b^12))/(a^2) - (768*(b^14))/(c^4) + (1280*(b^16))/((a^2)*(c^4)) - (512*(b^18))/((a^4)*(c^4)) - (1664*(b^12))/(c^2) + (2048*(b^14))/((a^2)*(c^2)) - (512*(b^16))/((a^4)*(c^2)) + (96*(b^10)*(f^2))/(s^2) + (96*(b^12)*(f^2))/((a^2)*(s^2)) + (64*(b^12)*(f^2))/((c^2)*(s^2));
c8 = 112*(b^8) - (32*(b^10))/(a^2) + (832*(b^12))/(c^4) - (1024*(b^14))/((a^2)*(c^4)) + (256*(b^16))/((a^4)*(c^4)) + (832*(b^10))/(c^2) - (576*(b^12))/((a^2)*(c^2)) - (48*(b^8)*(f^2))/(s^2) - (128*(b^10)*(f^2))/((a^2)*(s^2)) - (32*(b^12)*(f^2))/((a^4)*(s^2)) - (96*(b^10)*(f^2))/((c^2)*(s^2)) - (96*(b^12)*(f^2))/((a^2)*(c^2)*(s^2));
c6 = -16*(b^6) - (416*(b^10))/(c^4) + (288*(b^12))/((a^2)*(c^4)) - (224*(b^8))/(c^2) + (64*(b^10))/((a^2)*(c^2)) + (16*(b^6)*(f^2))/(s^2) + (48*(b^8)*(f^2))/((a^2)*(s^2)) + (32*(b^10)*(f^2))/((a^4)*(s^2)) + (48*(b^8)*(f^2))/((c^2)*(s^2)) + (128*(b^10)*(f^2))/((a^2)*(c^2)*(s^2)) + (32*(b^12)*(f^2))/((a^4)*(c^2)*(s^2));
c4 = (b^4) + (112*(b^8))/(c^4) - (32*(b^10))/((a^2)*(c^4)) + (32*(b^6))/(c^2) + ((b^4)*(f^4))/(s^4) - (2*(b^4)*(f^2))/(s^2) - (16*(b^6)*(f^2))/((a^2)*(s^2)) - (16*(b^6)*(f^2))/((c^2)*(s^2)) - (48*(b^8)*(f^2))/((a^2)*(c^2)*(s^2)) - (32*(b^10)*(f^2))/((a^4)*(c^2)*(s^2));
c2 = -((16*(b^6))/(c^4)) - (2*(b^4))/(c^2) - (2*(b^4)*(f^4))/((a^2)*(s^4)) + (2*(b^4)*(f^2))/((a^2)*(s^2)) + (2*(b^4)*(f^2))/((c^2)*(s^2)) + (16*(b^6)*(f^2))/((a^2)*(c^2)*(s^2));
c0 = (b^4)/(c^4) + ((b^4)*(f^4))/((a^4)*(s^4)) - (2*(b^4)*(f^2))/((a^2)*(c^2)*(s^2));
# Find the roots of the polynomial and tirn to velocity
Candidates=sqrt.(1 ./real(PolynomialRoots.roots([c0,c2,c4,c6,c8,c10,c12,c14,c16])));
# Delete all the values larger than the velocity in the Fluid
deleteat!(Candidates, Candidates .> c)
# test which root best satisfies stoneley equation
seq=zeros(length(Candidates));
for ii=1:length(Candidates)
    v = Candidates[ii];
    seq[ii] = sqrt((b/v)^2 - (b/c)^2)*((1 - 2*(b/v)^2)^2 - (4*(b/v)^2)*sqrt((b/v)^2 - (b/a)^2)*sqrt((b/v)^2 - 1)) + (f/s)*sqrt((b/v)^2 - (b/a)^2);
end
#
V_Stoneley=Candidates[argmin(seq)];
return V_Stoneley
end
    # calculate the fundamental
    # mode velocity of the guided wave for a model of a halfspace of water over
    # a halfspace of an elastic solid. This is called a Stoneley wave since its
    # velocity is less than the water velocity (i.e. it is trapped in both
    # directions, up and down). In contrast, a Scholte wave occurs for a finite
    # water depth when the guided wave velocity is greater than in water. The
    # Stoneley wave velocity is the solution of an eighth order polynomial.

## LOVE DISPERSION
" THIS IS THE MAIN LOVE WAVE DISPERSION FUNCTION"
##
function Love_Forwardsp(Nsolid::Int64,vsv::Array{Float64,1},rhov::Array{Float64,1},f::Array{Float64,1},fnum::Int64,hsolid::Array{Float64,1},Nmodes::Int64,qsv::Array{Float64,1})

# Some moduli
muv = @. rhov*vsv*vsv;
# Angular frecuency
# Calculate make angular frequency
ω = 2*pi*f;
## Initialize some matrixes
# initialize some local matrices
L1 = spzeros(2,2);
L3 = spzeros(2,2);
M1 = spzeros(2,2);
# initialize the global matrix
Ka1 = spzeros(Nsolid,Nsolid);
Ka3 = spzeros(Nsolid,Nsolid);
M = spzeros(Nsolid,Nsolid);
# For Solid Part of the Model
# Loop for the solid part of the Model
for ii=1:Nsolid
    # grab grid interval of current element
    h = hsolid[ii];
    #grab material properties of current element
    mu = muv[ii];
    # make elemental mass matrix
    M1 = spzeros(2,2);
    M1[1,1] = h*rhov[ii]*0.5;
    M1[2,2] = h*rhov[ii]*0.5;
    # make elemental stiffness matrices
    L1 = spzeros(2,2);
    L3 = spzeros(2,2);
    # some alternate variables from Lysmer
    alph = mu*0.16666666666666666;
    bet = mu*0.16666666666666666;
    # the 16 entries of the 4x4 elemental stiffness matrices of Lysmer
    L1[1,1] = 2*alph*h;
    L3[1,1] = (6*bet/h);
    L1[1,2] = alph*h;
    L3[1,2] = -(6*bet/h);
    L1[2,1] = L1[1,2];
    L3[2,1] = L3[1,2];
    L1[2,2] = L1[1,1];
    L3[2,2] = L3[1,1];
    # assemble mass and stiffness matrices from elemental matrices
    if (ii == Nsolid)
        M[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] = M[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] + M1[1:1,1:1];
        Ka1[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] = Ka1[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] + L1[1:1,1:1];
        Ka3[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] = Ka3[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] + L3[1:1,1:1];
    else
        M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + M1;
        Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L1;
        Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L3;
    end
end
# Set a lower bound on the Love wave speed
lspd = minimum(vsv);
# find the eigenvalue closest to the upper-bound eigenvalue
x=zeros(size(M,1),Nmodes,fnum);
d=zeros(Nmodes,1,fnum);
@simd for o=1:fnum
    dp,xp=eigs((ω[o]*ω[o]*M)-Ka3,Ka1,nev=Nmodes,sigma=(ω[o]/lspd)^2);
    d[:,:,o]=real.(dp);
    x[:,:,o]=real.(xp)
end
# Normalize the eigenfunction
ev=zeros(Nsolid,Nmodes,fnum)
for o=1:fnum
iev=zeros(Nsolid,Nmodes)
fctr=zeros(Nmodes,1)
    for i=1:Nmodes
        N=1/(x[1:Nsolid,i,o]'*M*x[1:Nsolid,i,o])
        fctr[i] = N[1];
        @. iev[:,i] = x[1:Nsolid,i,o]*sqrt(fctr[i]);
    end
    ev[:,:,o]=iev
end
# NOW! Some results!
# The wavenumber used is d
# The  Computed phase velocity (c)
vpk=zeros(Nmodes,1,fnum);
for o=1:fnum
    vpk[:,:,o] = @. ω[o]/sqrt.(d[:,:,o]); #√
end
# The Computed group velocity (U)
vgk=zeros(Nmodes,1,fnum);
for o=1:fnum
    ivgk=zeros(size(vpk,1),1);
    for i=1:size(vpk,1)
        ivgk[i] = (transpose(x[1:Nsolid,i,o])*(2*sqrt.(d[i,:,o]).*Ka1)*x[1:Nsolid,i,o]/(2*ω[o]))*(1/((transpose(x[1:Nsolid,i,o]))*M*x[1:Nsolid,i,o]));
    end
    vgk[:,:,o]=ivgk
end

# Atenuation Aproximation
# Working with the eigenfunctions
@simd for o=1:fnum
# The eigenfunctions
V3=ev[:,:,o] # SH-wave eigenfunction
# Components of the Lagrangian
I0=sum(rhov .*(V3.*V3) .*hsolid,dims=1);
#δcδvs
A=(vsv.*rhov)/(I0'.*vgk[:,:,o]);
TERM2=((([diff(V3,dims=1); zeros(1,Nmodes)])./hsolid).*(1 ./sqrt.(d[:,:,o]))').^2;
TERM1=V3.*V3;
B=sum((TERM1+TERM2).*hsolid,dims=1);
δcδvs=A.*B;
#δcδvp: Love wave do not depend on Vp velocity
# Attenuation Term
Qterm=sum((δcδvs.*(vsv*0.001)).*(1 ./qsv),dims=1)
# 𝛄 gamma
#TREF=50; # Reference period 50 s
ωr=0.12566370614359174#2*pi/TREF # T= 50 s as a reference
𝛄=(ω[o] ./(2*(vpk[:,:,o].*vpk[:,:,o]))).*Qterm';
# Anaelastic Phase velocity
c=vpk[:,:,o].+((0.3183098861837907*log(ω[o]/ωr))*Qterm)'
# Anaelastic Group velocity
A=2 .-(vgk[:,:,o]./vpk[:,:,o]);
B=(c-vpk[:,:,o])./vpk[:,:,o];
C=2*𝛄.*vgk[:,:,o]/pi*ω[o];
U=vgk[:,:,o].*(1 .+(A.*B) .+C);
# Safe the Values
vpk[:,:,o]=c;
vgk[:,:,o]=U;
end

return vpk, vgk
end

## RAYLEIGH DISPERSION
" THIS IS THE MAIN RAYLEIGH WAVE DISPERSION FUNCTION"
##
function Rayleigh_Forwardsp(Nsolid::Int64,vsv::Array{Float64,1},vpv::Array{Float64,1},rhov::Array{Float64,1},f::Array{Float64,1},fnum::Int64,hsolid::Array{Float64,1},Nmodes::Int64,Nfluid::Int64,vpfv::Array{Float64,2},rhofv::Array{Float64,2},hfv::Array{Float64,2},qsv::Array{Float64,1})
# Spherical Correction:
# ER=6371;
# HH=[0; hsolid]
# SHELLS=[cumsum(HH[1:end-1].*0.001) cumsum(HH[2:end].*0.001)]
# RAD=ER .-SHELLS
# CORRV=(2*ER) ./(RAD[:,1]+RAD[:,2]);
# CORRRHO_R=CORRV.^-2.275;
# vpv=vpv.*CORRV
# vsv=vsv.*CORRV
# rhov=rhov.*CORRRHO_R
# Fluid part of the model
# Check for the number of nodes in the fluid, based on the number of elements
if (Nfluid > 0)
    Nnfo = Nfluid + 1;
else
    Nnfo = 0;
end
# make fluid portion of model
# make kappa of the fluid , the modulus
kappafv = @. rhofv*vpfv*vpfv;
# Angular frecuency
# Calculate make angular frequency
ω = 2*pi*f;
# Initialize some matrixes
# initialize some local matrices
L1 = spzeros(2,2);
L3 = spzeros(2,2);
M1 = spzeros(2,2);
# initialize the global matrix
Ka1 = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
Ka3 = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
M = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
# MATRIXES IN THE FLUID PART OF THE Model
if Nfluid>0
# for all elements
    for ii=1:Nfluid
        # grab grid interval of current element
        h = hfluid[ii];
        # grab material properties of current element
        rhof = rhofv[ii];
        kappaf = kappafv[ii];
        # make elemental mass matrix
        M1 = spzeros(2,2);
        M1[1,1] = h/(2*kappaf);
        M1[2,2] = h/(2*kappaf);
        # make elemental stiffness matrices
        L1 = spzeros(2,2);
        L3 = spzeros(2,2);
        # some alternate variables from Lysmer
        alph = 1/(6*rhof);
        bet = 1/(6*rhof);
        # the 4 entries of the 2x2 elemental stiffness matrices of Lysmer
        L1[1,1] = 2*alph*h;
        L3[1,1] = (6*bet/h);
        L1[1,2] = alph*h;
        L3[1,2] = -(6*bet/h);
        L1[2,1] = L1[1,2];
        L3[2,1] = L3[1,2];
        L1[2,2] = L1[1,1];
        L3[2,2] = L3[1,1];
        # assemble mass and stiffness matrices from elemental matrices
        M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + M1;
        Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L1;
        Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L3;

    end
end
# Correct before solid part of Model
M[1,1] = M[1,1]*2;
Ka1[1,1] = Ka1[1,1]*2;
Ka3[1,1] = Ka3[1,1]*2;
# For Solid Part of the Model
# make solid portion of model
# make mu and lambda
muv = rhov.*vsv.*vsv;
lamdav = rhov.*vpv.*vpv - 2*muv;
#initialize some matrices
Ka2 = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
L1 = spzeros(4,4);
L2 = spzeros(4,4);
L3 = spzeros(4,4);
M1 = spzeros(4,4);
# Loop for the solid part of the Model
 for ii=1:Nsolid
    # grab grid interval of current element
    h = hsolid[ii];
    #grab material properties of current element
    mu = muv[ii];
    lamda = lamdav[ii];
    # make elemental mass matrix
    M1 = spzeros(4,4);
    M1[1,1] = h*rhov[ii]*0.5;
    M1[2,2] = h*rhov[ii]*0.5;
    M1[3,3] = h*rhov[ii]*0.5;
    M1[4,4] = h*rhov[ii]*0.5;
    # make elemental stiffness matrices
    L1 = spzeros(4,4);
    L2 = spzeros(4,4);
    L3 = spzeros(4,4);
    # some alternate variables from Lysmer
    alph = ((2*mu)+lamda)*0.16666666666666666;
    bet = mu/6;
    theta = (mu+lamda)*0.25;
    psi = (mu-lamda)*0.25;
    # the 16 entries of the 4x4 elemental stiffness matrices of Lysmer
    L1[1,1] = 2*alph*h;
    L3[1,1] = (6*bet/h);
    L2[1,2] = 2*psi;
    L1[1,3] = alph*h;
    L3[1,3] = -(6*bet/h);
    L2[1,4] = 2*theta;
    L2[2,1] = L2[1,2];
    L1[2,2] = 2*bet*h;
    L3[2,2] = (6*alph/h);
    L2[2,3] = -2*theta;
    L1[2,4] = bet*h;
    L3[2,4] = -(6*alph/h);
    L1[3,1] = L1[1,3];
    L3[3,1] = L3[1,3];
    L2[3,2] = L2[2,3];
    L1[3,3] = L1[1,1];
    L3[3,3] = L3[1,1];
    L2[3,4] = -2*psi;
    L2[4,1] = L2[1,4];
    L1[4,2] = L1[2,4];
    L3[4,2] = L3[2,4];
    L2[4,3] = L2[3,4];
    L1[4,4] = L1[2,2];
    L3[4,4] = L3[2,2];
    # assemble mass and stiffness matrices from elemental matrices
    if ii == Nsolid
        M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + M1[1:2,1:2];
        Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + L1[1:2,1:2];
        Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + L2[1:2,1:2];
        Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + L3[1:2,1:2];
    else
        M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + M1;
        Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + L1;
        Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + L2;
        Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + L3;
    end
end
# Construct the coupling matrix
if (Nfluid > 0)
Cm = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
Cm[Nnfo,Nnfo+2] = 1;
Cm[Nnfo+2,Nnfo] = 1;
else
    Cm = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
end
# find the rayleigh/scholte wave speed which would exist if the solid model
# were a halfspace with the minimum model velocity
# this is a lower bound on the velocity that can be passed to EIGS, based
# on ARPACK
if (Nfluid > 0)
    msval  = minimum(vsv);
    msloc = argmin(vsv);
    mfval  = minimum(vpfv);
    mfloc = argmin(vpfv);
    vsmay = msval;
    vpmay = vpv[msloc];
    vpfmay = mfval;
    rhofmay = rhofv[mfloc];
    rhomay = rhov[msloc];
    rspd = Stoneley(vpmay,vsmay,vpfmay,rhofmay,rhomay);
else
    msval  = minimum(vsv);
    msloc = argmin(vsv);
    vsmay = msval;
    vpmay = vpv[msloc];
    # coefficients of rayleigh's polynomial
    t1 = 1/(vsmay^6);
    t2 = -8/(vsmay^4);
    t3 = ((24/(vsmay^2))-(16/(vpmay^2)));
    t4 = -16*(1-((vsmay/vpmay)^2));
    # rayleigh wave speed
    rspd = sqrt(minimum(real(PolynomialRoots.roots([t4, t3, t2, t1]))));
end
# Find the eigenvalue closest to the upper-bound eigenvalue
x=zeros(2*size(M,1),Nmodes,fnum);
d=zeros(Nmodes,1,fnum);
xp=zeros(2*size(M,1),Nmodes);
dp=zeros(Nmodes,1);
@simd for o=1:fnum
    dp,xp=eigs([spzeros((Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))) sparse(I,(Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))); ((ω[o]*ω[o]*M)-Ka3-(ω[o]*Cm)) Ka2],[sparse(I,(Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))) spzeros((Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))); spzeros((Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))) Ka1],nev=Nmodes,sigma=ω[o]/rspd);
    x[:,:,o] = real.(xp);
    d[:,:,o] = real.(dp);
end

# Normalize the eigenfunctions
ev=zeros(size((Nnfo+1):(Nnfo+(2*Nsolid)),1),Nmodes,fnum);
for o=1:fnum
    iev=zeros(size((Nnfo+1):(Nnfo+(2*Nsolid)),1),Nmodes);
    for e=1:Nmodes
        fctr = 1/ (((x[1:1:(Nnfo+(2*Nsolid)),e,o]')*M*x[1:1:(Nnfo+(2*Nsolid)),e,o])-(((x[1:1:(Nnfo+(2*Nsolid)),e,o]')*Cm*x[1:1:(Nnfo+(2*Nsolid)),e,o])/2*ω[o]));
        evp = x[1:1:(Nnfo+(2*Nsolid)),e,o]*sqrt(fctr)*sign(x[Nnfo+1,e,o]);
        # # return only the eigenvector in the solid
        iev[:,e] = evp[(Nnfo+1):(Nnfo+(2*Nsolid))];
    end
    ev[:,:,o]=iev
end
# NOW! Some results!
# The wavenumber used is d
# The  Computed phase velocity (c)
vpk=zeros(Nmodes,1,fnum);
for o=1:fnum
    vpk[:,:,o] = @. ω[o]/d[:,:,o];
end
# The Computed group velocity (U)
vgk=zeros(Nmodes,1,fnum);
for o=1:fnum
    ivgk=zeros(size(vpk,1),1);
    for i=1:size(vpk,1)
        ivgk[i] = (((x[1:1:(Nnfo+(2*Nsolid)),i,o]')*((2*d[i,1,o]*Ka1)-Ka2)*x[1:1:(Nnfo+(2*Nsolid)),i,o])/((2*ω[o]*((x[1:1:(Nnfo+(2*Nsolid)),i,o]')*M*x[1:1:(Nnfo+(2*Nsolid)),i,o]))-((x[1:1:(Nnfo+(2*Nsolid)),i,o]')*Cm*x[1:1:(Nnfo+(2*Nsolid)),i,o])));
    end
    vgk[:,:,o]=ivgk
end
# Atenuation Aproximation
# Working with the eigenfunctions
@simd for o=1:fnum
    # The eigenfunctions
V1=ev[2:2:end,:,o] # vertical Eigenfcuntions
V2=ev[1:2:end,:,o] # radial eigenfunctions
# Components of the Lagrangian
I0=sum(rhov.*(V1.*V1 .+ V2.*V2).*hsolid,dims=1);
#δc/δvp
A=(vpv.*rhov)/(I0'.*vgk[:,:,o]);
derivativeV1=(([diff(V1,dims=1); zeros(1,Nmodes)])./hsolid) .*(1 ./d[:,:,o])'
B=sum(((V2-derivativeV1).^2).*hsolid,dims=1)
δcδvp=A.*B;
#δc/δvs
A=(vsv.*rhov)/(I0'.*vgk[:,:,o]);
derivativeV1=(([diff(V1,dims=1); zeros(1,Nmodes)])./hsolid) .*(4 ./d[:,:,o])'
derivativeV2=(([diff(V2,dims=1); zeros(1,Nmodes)])./hsolid) .*(1 ./d[:,:,o])'
TERM1=(V1+derivativeV2).*(V1+derivativeV2)
TERM2=V2.*derivativeV1
B=sum((TERM1+TERM2).*hsolid,dims=1)
δcδvs=A.*B;
# Attenuation Term
Qterm=sum((δcδvs.*(vsv*0.001)).*(1 ./qsv),dims=1)
# 𝛄 gamma
ωr=0.12566370614359174#2*pi/TREF # T= 50 s as a reference
𝛄=(ω[o] ./(2*(vpk[:,:,o].*vpk[:,:,o]))).*Qterm';
# Anaelastic Phase velocity
c=vpk[:,:,o]+((0.3183098861837907*log(ω[o]/ωr))*Qterm)'
# Anaelastic Group velocity
A=2 .-(vgk[:,:,o]./vpk[:,:,o]);
B=(c-vpk[:,:,o])./vpk[:,:,o];
C=2*𝛄.*vgk[:,:,o]/pi*ω[o];
U=vgk[:,:,o].*(1 .+(A.*B) .+C);
# Safe the Values
vpk[:,:,o]=c;
vgk[:,:,o]=U;
end

# Return Values
return vpk, vgk
end

##
# function fix_SWF_input(TF,Tstep,T0,h,Model)
# T=collect(TF:-Tstep:T0);
# fks=1 ./T; # vector of  frequencies  at which the velocities are measured (Hz) Change T (s) range
# Step=1000; #step for the fluid
# Mesh=cumsum(h,dims=1); #from special h to the FEM mesh
# Mesh=Mesh[1:end-1]
#
# # Iterpolate all the properties
# iVp=LinearInterpolation(Model[:,1]*1000,  Model[:,2]*1000, extrapolation_bc=Flat());
# iVs=LinearInterpolation(Model[:,1]*1000,  Model[:,3]*1000, extrapolation_bc=Flat());
# iρ=LinearInterpolation(Model[:,1]*1000,  Model[:,4]*1000, extrapolation_bc=Flat());
# vpv=iVp(Mesh*1000);
# vsv=iVs(Mesh*1000);
# rhov=iρ(Mesh*1000);
#
# # Atenuaction Qs
# iQs=LinearInterpolation(Model[:,1]*1000,  Model[:,5], extrapolation_bc=Flat());
# qsv=iQs(Mesh*1000);
#
# # Radial Anysotropy
# iR=LinearInterpolation(Model[:,1]*1000,  Model[:,6], extrapolation_bc=Flat());
# qR=iR(Mesh*1000).*0.01;
# vsv_v=((3 .-qR)*0.3333333333333333).*vsv; #Vertical velocity in realtion to Radial Anisitropy
# vsv_h=(1 .+ 0.6666666666666666.*qR).*vsv; #Horizontal velocity equal to total velocity
#
#
# # FEM MODEL
# # The model begins with the number of elements
# # SOLID PART (the code can consider the water layer)
# Nsolid = size(vpv,1); # number of elements in solid
# hsolid =[diff(Mesh); Model[end,1]-Mesh[end,1]]*1000;
# # Fluid Part
# Nfluid = 0; # number of elements in fluid
# hfluid = Step*ones(1,Nfluid); # grid spacing of mesh (meters)
# vpfluid = 1500; rhofluid = 1030;
# vpfv = vpfluid*ones(1,Nfluid);
# rhofv = rhofluid*ones(1,Nfluid);
#
# # Spherical Correction:
# ER=6371;
# HH=[0; hsolid]
# SHELLS=[cumsum(HH[1:end-1].*0.001) cumsum(HH[2:end].*0.001)]
# RAD=ER .-SHELLS
# CORRV=(2*ER) ./(RAD[:,1]+RAD[:,2]);
#
# #Apply to R wave data
# CORRRHO_R=CORRV.^-2.275;
# vpv=vpv.*CORRV
# vsv_v=vsv_v.*CORRV
# rhovR=rhov.*CORRRHO_R
#
# #Apply to L wave data
# CORRRHO_L=CORRV.^-5;
# vsv_h=vsv_h.*CORRV
# rhovL=rhov.*CORRRHO_L
#
# # frecuency
# fnum=size(fks,1);
#
# # other useful Components
#
# return fks,Mesh,vpv,vsv,vsv_v,vsv_h,rhovR,rhovL,qsv,qR,Nsolid,hsolid,Nfluid,hfluid,vpfluid,vpfv,rhofv,T,fnum
# end
