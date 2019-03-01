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

# #Speedup
# if (~issparse(A))
    # #If not sparse they might be single, this func is not going to work if so.
    # if (isa(A,"single"))
        # A=double(A);
    # end
    # if (isa(b,"single"))
        # b=double(b);
    # end
  # A = sparse(A);	# Make sure M is sparse
# end
# b = full(b(:)); 	# Make sure q is a column vector


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
solver = "random"; #Gets switched to zero lower down...
profile = false;
#switch lower(solver)
#	case "penalized"
#		lambda = 0.95;
#	otherwise
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
zers=copy(J)
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
dx=zeros(N);
totaltime=0.;
S=copy(dx);
phi=copy(dx);
test=[0.0];
grad_f=[0.0];
I=zeros(Int64, N)
I=convert(Array{Bool,1},I) #convert to bool
p=copy(y);
q=copy(y);
x2=copy(y);
y2=copy(y);

#tic=time()	
#toc=time()
#println("Elapsed time")
#println(toc-tic)	

while (iter <= max_iter )
    #println("Looping")
	
	#y = A*x (pre allocated y)
	mul!(y,A,x) 
	y.=y.+b
	

	
	#--- Test all stopping criteria used ------------------------------------
	phi     = phi_lambda!(y, x, lambda,phi,phi_l);         # Calculate fischer function
	old_err = err;
	for i=1:length(phi); phiT[i]=phi[i]; end #transpose
	err     = 0.5*dot(phiT,phi);       # Natural merit function
	err=err[1]; #extract
	
	if profile
		convergence = [convergence err];
	end
	
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
	I.= (S.==false);                     		# Bitmask for non-singular indices
	
	Indx=findall(I)

	#@code_warntype
	WorkOnJ(J,A,x,y,Indx)

	
	if min(size(A,1),50)/2<30
		restart=min(size(A,1),50)/2;  
	else
		restart = 30; 
	end  
	restart	= convert(Int64,restart)

	#println("Creating sparse MAT")
	Jsp=SparseArrays.sparse(J); 
	#dropzeros!(Jsp)
	#println("Sparse MAT created")
	
	dx.=dx.*0; #Reset Newton direction
	
	#println("Into solver")
	phiM.=.-phi;
	#@time phiM=sparsevec(phiM)
	
	##Solvers:
	IterativeSolvers.gmres!(dx,Jsp,phiM,initially_zero=true,restart=restart,tol=1e-6);
	#singleloopt=@elapsed IterativeSolvers.idrs!(dx,Jsp,phiM,s=8); #8 is good
	#singleloopt=@elapsed IterativeSolvers.bicgstabl!(dx,Jsp,phiM,10); #NOT CONSISTENT!!!!! errors every few runs with test
	#singleloopt=@elapsed FUNC
	#totaltime=totaltime+singleloopt;
	#println(totaltime)
	
	#tic=time()	

	
	############################################

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
	
	
	#Reset J to zeros
	#J[Indx,Indx].=0.0; #fill!(J,0.0)
	for i=eachindex(Indx)
		 for j=eachindex(Indx)
			 J[Indx[i],Indx[j]]=0.0;
		 end
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
		if ( f_k[1] <= test[1])
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

	#println("Elapsed time")
	#toc=time()
	#println(toc-tic)	
	
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
