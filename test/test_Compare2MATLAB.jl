#Start creating vars for function: 
printstyled("creating func vars \n",color=:cyan)
#using Profile

global G
G=pathof(FischerNewton);
G=splitdir(G); #remove file name
G=G[1];
G=splitdir(G); #out of src
G=G[1];

#If we read mats in as text files:
using DelimitedFiles

#Would be nicer to do programatically
if Sys.iswindows()
	G=string(G,"\\test\\")
else
	G=string(G,"/test/")
end
GA=string(G,"Matricies-A.txt")
Gb=string(G,"Matricies-b.txt")
Gx=string(G,"Matricies-x.txt")
Gy=string(G,"Matricies-y.txt")
GDn=string(G,"Matricies-Dn.txt")
GDss=string(G,"Matricies-Dss.txt")
GDds=string(G,"Matricies-Dds.txt")
GTn=string(G,"Matricies-Tn.txt")
GTss=string(G,"Matricies-Tss.txt")
GTds=string(G,"Matricies-Tds.txt")
A=readdlm(GA, ',', Float64)
b=readdlm(Gb, ',', Float64)
xMATLAB=readdlm(Gx, ',', Float64)
yMATLAB=readdlm(Gy, ',', Float64)
DnMATLAB=readdlm(GDn, ',', Float64)
DssMATLAB=readdlm(GDss, ',', Float64)
DdsMATLAB=readdlm(GDds, ',', Float64)
TnMATLAB=readdlm(GTn, ',', Float64)
TssMATLAB=readdlm(GTss, ',', Float64)
TdsMATLAB=readdlm(GTds, ',', Float64)


#If we use MAT:
#=
if Sys.iswindows()
    G=string(G,"\\test\\Matricies.mat")
else
	G=string(G,"/test/Matricies.mat")
end
println(G)
#Inputs - Loading influence matricies
using MAT
file = matopen(G)
A=read(file, "A")
b=read(file, "b")
#MATLAB result
xMATLAB=read(file, "x")
yMATLAB=read(file, "y")
DnMATLAB=read(file, "Dn")
DssMATLAB=read(file, "Dss")
DdsMATLAB=read(file, "Dds")
TnMATLAB=read(file, "Tn")
TssMATLAB=read(file, "Tss")
TdsMATLAB=read(file, "Tds")
=#





printstyled("Vars loaded -> to fischerNewton func \n",color=:cyan)
#@profile (x)=FischerNewton.fischer_newton(A,b);
#Profile.print(format=:tree)#(format=:flat) (sortedby=:count)

@time (x)=FischerNewton.fischer_newton(A,b);
printstyled("Out of func \n",color=:cyan)

y = A*x+b;

printstyled("Checking complimentarity conditions are met \n",color=:cyan)

#X*Y=0 
using LinearAlgebra
xyMatlab=LinearAlgebra.dot(xMATLAB,yMATLAB)
println("sum(x.*y) MATLAB result")
@info xyMatlab
xyJulia=LinearAlgebra.dot(x,y)
println("sum(x.*y) Julia result")
@info xyJulia
if abs(xyJulia)>abs(xyMatlab)
	printstyled("Julia is not as good as MATLAB at meeting comp conditions \n",color=:light_red)
else
	printstyled("Julia better than MATLAB at meeting comp conditions \n",color=:green)
end

#y>0
if any(y.<0)
	maxneg=maximum(y[y.<0])
	if abs(maxneg)>1e-10 #some tolerance allowed
	@info maxneg
		error("negative y values")
	end
end
#x>0
if any(x.<0)
	maxneg=maximum(x[x.<0])
	@info maxneg
	if abs(maxneg)>1e-10 #some tolerance allowed
		error("negative x values")
	end
	
end

printstyled("Complimentarity conditions are met \n",color=:light_green)
printstyled("Checking residuals relative to MATLAB's result \n",color=:cyan)

#@info typeof(x[:]) 
#@info size(x[:])
#@info size(xMATLAB[:])
#@info x[1:10] xMATLAB[1:10]
#@info x[end-10:end] xMATLAB[end-10:end]
#@info y[end-10:end] yMATLAB[end-10:end]

(xRes)=FischerNewton.CompResiduals(x,xMATLAB);
(yRes)=FischerNewton.CompResiduals(y,yMATLAB);
#println("Values of residuals")
#@info xRes 
#@info yRes


# Creating some lengths we can use for extraction of sub matricies/vectors
# ne being the number of elements in one edge a submatrix/vector
ne=length(x)/5; 
L1=1:ne; 			L1=Int64.(L1); #Convert 2 int
L2=ne+1:2*ne;		L2=Int64.(L2); #Convert 2 int
L3=2*ne+1:3*ne;		L3=Int64.(L3); #Convert 2 int
L4=3*ne+1:4*ne;		L4=Int64.(L4); #Convert 2 int
L5=4*ne+1:5*ne;		L5=Int64.(L5); #Convert 2 int

#Extracting sub parts of vectors: 
#[x]
x1=x[L1];
x2=x[L2];
x3=x[L3];
x4=x[L4];
x5=x[L5];
#[y]
y1=y[L1];
y2=y[L2];
y3=y[L3];
y4=y[L4];
y5=y[L5];
                                                  
#Extracting slip                                      
Dn  = y1;                            
Dss = y2.-x4;     
Dds = y3.-x5;
#Extracting traction  
Tn=-x1;
Tss=y4.-x2;
Tds=y5.-x3;

println("Values of residuals (Disp)")
(DnRes)=FischerNewton.CompResiduals(Dn,DnMATLAB);
(DssRes)=FischerNewton.CompResiduals(Dss,DssMATLAB);
(DdsRes)=FischerNewton.CompResiduals(Dds,DdsMATLAB);
@info DnRes
@info DssRes
@info DdsRes

println("Values of residuals (Traction)")
(TnRes)=FischerNewton.CompResiduals(Tn,TnMATLAB);
(TssRes)=FischerNewton.CompResiduals(Tss,TssMATLAB);
(TdsRes)=FischerNewton.CompResiduals(Tds,TdsMATLAB);
@info TnRes
@info TssRes
@info TdsRes

if DnRes>1E-5
	error("Residual Dn too high")
end
if DssRes>1E-5
	error("Residual Dss too high")
end
if DdsRes>1E-5
	error("Residual Dds too high")
end
if TnRes>1E-5
	error("Residual Tn too high")
end
if TssRes>1E-5
	error("Residual Tss too high")
end
if TdsRes>1E-5
	error("Residual Tds too high")
end
printstyled("Residuals appear OK (according to the arbitary limits set) \n",color=:light_green)

