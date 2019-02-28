#Start creating vars for function: 
println("creating func vars")



#Inputs - Loading influence matricies
using MAT
file = matopen(raw"C:\Users\timmm\Desktop\FisherNewtonJulia\Matricies.mat")
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

println("Vars loaded -> to fischerNewton func")

#using Profile
@time (x)=FischerNewton.fischer_newton(A,b);
#poop
#Profile.print(format=:flat)# (sortedby=:count)

y = A*x+b;

println("Out of func")

@info typeof(x[:]) 
@info size(x[:])
@info size(xMATLAB[:])
@info x[1:10] xMATLAB[1:10]
@info x[end-10:end] xMATLAB[end-10:end]
@info y[end-10:end] yMATLAB[end-10:end]

(xRes)=FischerNewton.CompResiduals(x,xMATLAB);
(yRes)=FischerNewton.CompResiduals(y,yMATLAB);
println("Values of residuals")
@info xRes 
@info yRes


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

if DnRes>1E-6
	error("Residual too high")
end
if DssRes>1E-6
	error("Residual too high")
end
if DdsRes>1E-6
	error("Residual too high")
end
if TnRes>1E-6
	error("Residual too high")
end
if TssRes>1E-6
	error("Residual too high")
end
if TdsRes>1E-6
	error("Residual too high")
end

println("MaxDn")
@info maximum(Dn)