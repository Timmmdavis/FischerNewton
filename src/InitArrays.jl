function InitArrays(n)
#n - length of one part of the vector b (5*n vector for 3D)

N    = n; # Number of variables
flag = 1;

#--- Make sure we got good working default values -------------------------
x = zeros(N);
max_iter = floor(N/2);
tol_rel = 0.0001;
tol_abs = 10*eps(); # Order of 10th of numerical precision seems okay
lambda = 1.00; # No support for penalization in other solvers

#--- Make sure all values are valid ---------------------------------------
max_iter = max(max_iter,1.0);
tol_rel  = max(tol_rel,0.0);
tol_abs  = max(tol_abs,0.0);
x.= max.(x,0.0);

#--- Here comes a bunch of magic constants --------------------------------

#Using the 'zero' solver as described in the original function
#solver="zero";   

#doesn't make much diff
# h       = 1e-7;    # Fixed constant used to evaluate the directional detivative - not used in zero solver
alpha   = 0.5;     # Step reduction parameter for projected Armijo backtracking line search
beta    = 0.001;   # Sufficent decrease parameter for projected Armijo backtracking line search
gamma   = 1e-28;   # Perturbation values used to fix near singular points in derivative
rho     = eps();     # Descent direction test parameter used to test if the Newton direction does a good enough job.

#--- Setup values need while iterating ------------------------------------

err = [Inf];         # Current error measure
iter= 1;           # Iteration count
nabdx 	=[0.]
old_err =err;
test=[0.0];
grad_f=[0.0];
f_k=[0.0]


#vectors
y 	=zeros(N);
phi =copy(y);
phiT=zeros(1,N);
phi_k 	=copy(y); 
phi_kT 	=zeros(1,N);
phi_l 	=copy(y); 
phiM 	=copy(y);
dx 	=copy(y);  
absdx 	=copy(y); 
y_k =copy(y);
xdxtau 	=copy(y);
x_k 	= copy(y);
nabla_phi=similar(phiT);


I=zeros(Int64, N)
I=convert(Array{Bool,1},I) #convert to bool
J = zeros(N,N); #spzeros

Steps = 1:1:N::Int64
II=repeat(Steps,1,N);

Jsubs=copy(J)

#totaltime1=0.;
#totaltime2=0.;
#tic=time()		
#toc=time()
#println("Elapsed time")
#println(toc-tic)

useSparse=0

if useSparse==1
	#SparseMat stuff:
	JJ=transpose(II);
	JJ[Abad].=0; #Bad bits are zero (dropped later)
	#Preallocate some vectors to work with
	ISml=zeros(N^2)
	JSml=zeros(N^2)
	VSml=zeros(N^2)
end

Vects=Vectors(y, phi, phi_k, phi_l, phiM, dx, absdx, y_k, xdxtau, x_k, x)
Arrys=Arrays(err, nabdx, test, grad_f, f_k)
Flts=Float(alpha,beta,gamma,rho,max_iter,tol_rel,tol_abs,lambda)
Ints=Int(iter,N,flag,useSparse)

if useSparse==1
	Mats=Matricies(J, Jsubs, JJ, ISml, JSml, VSml, phiT, phi_kT, nabla_phi)
else
	arry= zeros(2,2)::Array
	Mats=Matricies(J, Jsubs,arry,arry,arry,arry, phiT, phi_kT, nabla_phi)
end
Bls=Bools(I)
IntArrys=IntArrays(II)

return Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys

end