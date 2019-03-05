function fischer_newton(A,b)
#Note this is an edited version of the code NUM4LCP.m
#Its recommended you use the original found at: 
#(Niebe 2016) = Niebe, S. and Erleben, K., 2015. Numerical methods for linear complementarity problems in physics-based animation. Synthesis Lectures on Computer Graphics and Animation, 7(1), pp.1-159.
# Copyright 2017, Tim Davis, The University of Aberdeen/Potsdam
# Copyright 2011, Kenny Erleben, DIKU
# Copyright 2012, Michael Andersen, DIKU
#
#	Original copyright below (Tim Davis, 2017): 
##
#Copyright (c) 2011 Kenny Erleben code.google.com/p/num4lcp/
#
#This software is provided 'as-is', without any express or implied
#warranty. In no event will the authors be held liable for any damages
#arising from the use of this software.
#
#Permission is granted to anyone to use this software for any purpose,
#including commercial applications, and to alter it and redistribute it
#freely, subject to the following restrictions:
#
#    1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
#
#    2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
#
#    3. This notice may not be removed or altered from any source
#    distribution.
##

if length(size(b))>1
	println("b should not have two dimensions, reshaping for you...")
	b=reshape(b,length(b));
end

N    = length(b); # Number of variables
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

#doesn't make much diff
h       = 1e-7;    # Fixed constant used to evaluate the directional detivative
alpha   = 0.5;     # Step reduction parameter for projected Armijo backtracking line search
beta    = 0.001;   # Sufficent decrease parameter for projected Armijo backtracking line search
gamma   = 1e-28;   # Perturbation values used to fix near singular points in derivative
rho     = eps();     # Descent direction test parameter used to test if the Newton direction does a good enough job.

#--- Setup values need while iterating ------------------------------------

convergence = []; # Used when profiling to measure the convergence rate

err     = Inf;         # Current error measure
#x         				# Current iterate
iter    = 1;           # Iteration count

solver="zero";   

println("check all arrays below are used")
# Init vars early on
J = zeros(N,N); #spzeros 
old_err=err;
#vectors
y=zeros(N);
phi=copy(y);
phiT=zeros(1,N);
phi_k=copy(y); 
phi_kT=zeros(1,N);
phi_l=copy(y); 
phiM=copy(y);
dx=copy(y);   
S=copy(y);
y_k=copy(y);
nabla_phi=similar(phiT);
totaltime1=0.;
totaltime2=0.;
test=[0.0];
grad_f=[0.0];
I=zeros(Int64, N)
I=convert(Array{Bool,1},I) #convert to bool

Steps = 1:1:N::Int64
II=repeat(Steps,1,N);
JJ=transpose(II);

ASpc=sparse(A);
(AI,AJ,AV)=findnz(ASpc);
ASpcbad=findall(iszero, ASpc);
II[ASpcbad].=0;
JJ[ASpcbad].=0; #Bad bits are zero (drop later)
#ISml=zeros(N^2);
#JSml=zeros(N^2);
#VSml=zeros(N^2);

#tic=time()	
#toc=time()
#println("Elapsed time")
#println(toc-tic)
	
Indx=0;
while (iter <= max_iter )
    #println("Looping")
	
	
	#y = A*x (pre allocated y)
	mul!(y,A,x) 
	y.=y.+b
	
	#--- Test all stopping criteria used ------------------------------------
	phi = phi_lambda!(y, x, lambda,phi,phi_l);         # Calculate fischer function
	
	old_err = err;
	for i=1:length(phi); phiT[i]=phi[i]; end #transpose
	err     = 0.5*dot(phiT,phi);       # Natural merit function
	err=err[1]; #extract
	
	if (abs(err-old_err) / abs(old_err)) < tol_rel  # Relative stopping criteria
		flag = 3;
		break;
	end
	if err < tol_abs   # Absolute stopping criteria
		flag = 4;
		break;
	end

	#--- Solve the Newton system --------------------------------------------
	#--- First we find the Jacobian matrix. We wish to avoid singular points
	#--- in Jacobian of the Fischer function
	S.= (abs.(phi).<gamma) .& (abs.(x).<gamma);  # Bitmask for singular indices
	I.= (S.==false);        
		

		
	#Function that creates sparse MAT J	
	#J1=WorkOnJ(J,A,x,y,I)
	println("Precompute J total time")
	singleloopt=@elapsed J=WorkOnJ_FastBigMats(A,x,y,I,II,JJ)	
	totaltime1=totaltime1+singleloopt;
	println(totaltime1)
	
	#If you want to compare the outputs of the methods above:
	 # (I1,J1,V1)=findnz(J1);
	 # (I2,J2,V2)=findnz(J2);
	  # if !isequal(V1,V2)
		  # @info size(V1) size(V2)
		  # error("Not eq")
	  # end	
	
	if min(size(A,1),50)/2<30
		restart=min(size(A,1),50)/2;  
	else
		restart = 30; 
	end  
	restart	= convert(Int64,restart)

	dx.=dx.*0; #Reset Newton direction
	phiM.=.-phi;
	
	dxSubset=dx[I];
	JSubset=J[I,I];
	phiMSubset=phiM[I]
	
	println("Total elapsed solver time")
	singleloopt=@elapsed IterativeSolvers.gmres!(dxSubset,JSubset,phiMSubset,tol=1e-6,restart=restart, initially_zero=true,maxiter=10);
	totaltime2=totaltime2+singleloopt;
	println(totaltime2)
	
	dx[I]=dxSubset;
	
	#IterativeSolvers.idrs!(dx,Jsp,phiM,s=8); #8 is good
	#IncompleteLU.LU = ilu(J, τ = 0.1);
	#IterativeSolvers.bicgstabl!(dx,J,phiM,2,Pl = LU); #not good. 
	#singleloopt=@elapsed FUNC

	
	
	# Test if the search direction is smaller than numerical precision. 
	# That is if it is too close to zero.
	if maximum(abs.(dx)) < eps()
		flag = 5;
		# Rather than just giving up we may just use the gradient direction
		# instead. However, I am lazy here!
		#  dx = nabla_phi'
		break;
	end
	
	#@info size(phiT) size(J)
	# Test if we have dropped into a local minimia if so we are stuck
	#nabla_phi = phiT*J;
	mul!(nabla_phi,phiT,J) 		
	if norm(nabla_phi) < tol_abs
		flag = 6;
		break;
	end
	

	# Test if our search direction is a 'sufficient' descent direction
	nabdx=nabla_phi*dx
	if  nabdx[1]  > -rho*dot(dx',dx)
		flag = 7;
		# Rather than just giving up we may just use the gradient direction
		# instead. However, I am lazy here!
		#  dx = nabla_phi'
		break;
	end
	

	
	#--- Armijo backtracking combined with a projected line-search ---------
	tau     = 1.0;                  # Current step length
	f_0     = err;
	grad_f= beta*dot(nabla_phi,dx);
	x_k     = x;
		
	while true #Inf loop, escapes when a break is performed
		x_k.=max.(0.0,x.+dx.*tau) #x_k   = max.(0.0,x + dx*tau); #non negativity
		#y_k   = A*x_k + b; (pre allocated y_k)
		mul!(y_k,A,x_k) 
		y_k.=y_k.+b
		phi_k = phi_lambda!( y_k, x_k, lambda,phi_k,phi_l );	
		
		for i=1:length(phi_k); phi_kT[i]=phi[i]; end #transpose
		f_k = 0.5*dot(phi_kT,phi_k);       # Natural merit function

		
		# Perform Armijo codition to see if we got a sufficient decrease
		test=f_0.+ tau*grad_f;
		if ( f_k <= test)
			break;
		end

		# Test if time-step became too small
		if tau*tau < gamma
			break;
		end	

		tau = alpha*tau;
	end #end of while lp 1. 

	
	# Update iterate with result from Armijo backtracking
	x.= x_k;

	# Increment the number of iterations
	iter = iter + 1;	

end #end of while lp2. 


	
if iter >= max_iter
	flag = 8;
	iter = iter - 1;
end		

return(x); #, err, iter, flag, convergence)

end #Func end


function fischer(y,x)
# Auxiliary function used by the Fischer-Newton method
# Copyright 2011, Kenny Erleben

phi  = (y.^2 + x.^2).^0.5 - y - x;
return(phi)
end

function fischer!(y,x,phi)
# Auxiliary function used by the Fischer-Newton method
# Copyright 2011, Kenny Erleben

for i=1:length(y)
	phi[i]= (y[i]^2 + x[i]^2)^0.5 - y[i] - x[i];
end

return(phi)
end


function phi_lambda(a,b,lambda)
## CCK NCP-function, a convex composition of Fischer-Burmeister and CCK NCP
#
#   Input
#       a -> A column vector size = (n,1)
#       b -> A column vector size = (n,1)
#       l -> A fixed lambda value used to weight the input.
#
#   Output
#       phi_l -> A column vector with the result of the Fischer-Burmeister
#                NCP function, with size = (n,1)

zers=zeros(size(a));
phi_l = lambda*fischer(a,b)+(1.0-lambda)*(max.(a,zers).*max.(b,zers));

return(phi_l)
end

function phi_lambda!(a,b,lambda,phi,phi_l)
## CCK NCP-function, a convex composition of Fischer-Burmeister and CCK NCP
#
#   Input
#       a -> A column vector size = (n,1)
#       b -> A column vector size = (n,1)
#       l -> A fixed lambda value used to weight the input.
#
#   Output
#       phi_l -> A column vector with the result of the Fischer-Burmeister
#                NCP function, with size = (n,1)

#updates phi in place
fischer!(a,b,phi)
for i=1:length(phi)
	phi_l[i] = lambda*phi[i]+(1.0-lambda)*(max(a[i],0.0)*max(b[i],0.0));
end

return(phi_l)

end
