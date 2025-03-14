function fischer_newton(A,b,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
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

#Extract out of mat
alpha=Flts.alpha
beta=Flts.beta
gamma=Flts.gamma
rho=Flts.rho
max_iter=Flts.max_iter
tol_rel=Flts.tol_rel
tol_abs=Flts.tol_abs
lambda=Flts.lambda

iter=Ints.iter
N=Ints.N
flag=Ints.flag
useSparse=Ints.useSparse

y=Vects.y 
phi=Vects.phi
phi_k=Vects.phi_k
phi_l=Vects.phi_l
phiM=Vects.phiM
dx=Vects.dx
absdx=Vects.absdx
y_k=Vects.y_k
xdxtau=Vects.xdxtau
x_k=Vects.x_k
x=Vects.x

err=Arrys.err
nabdx=Arrys.nabdx
test=Arrys.test
grad_f=Arrys.grad_f
f_k=Arrys.f_k
old_err=Arrys.old_err

J=Mats.J
Jsubs=Mats.Jsubs
JJ=Mats.JJ
ISml=Mats.ISml
JSml=Mats.JSml
VSml=Mats.VSml
phiT=Mats.phiT
phi_kT=Mats.phi_kT
nabla_phi=Mats.nabla_phi

I=Bls.I
II=IntArrys.II

Indx=0;
while (iter <= max_iter )
	Indx+=1;
    #println("LoopNo")
	#println(Indx)
	
	#y = A*x+b 
	mul!(y,A,x) 
	for i=1:N
		y[i]+=b[i]
	end

	#--- Test all stopping criteria used ------------------------------------
	phi = phi_lambda!(y, x, lambda,phi,phi_l);         # Calculate fischer function
	
	old_err[1] = err[1];
	for i=1:length(phi); phiT[i]=phi[i]; end #transpose
	# Natural merit function
	mul!(err,phiT,phi)
	err[1]=err[1]*0.5     

	if (abs(err[1]-old_err[1]) / abs(old_err[1])) < tol_rel  # Relative stopping criteria
		flag = 3;
		break;
	end
	if err[1] < tol_abs   # Absolute stopping criteria
		flag = 4;
		break;
	end

	#--- Solve the Newton system --------------------------------------------
	#--- First we find the Jacobian matrix. We wish to avoid singular points
	#--- in Jacobian of the Fischer function
	for i=1:N
		# Bitmask for singular indices
		test1=abs(phi[i])<gamma
		test2=abs(x[i])<gamma
		if test1 & test2
			I[i]=false    
		else
			I[i]=true
		end
	end
	
	#Function that creates Matrix J (sparse if using 2nd func)	
	if useSparse==0
		J=WorkOnJ(J,A,x,y,I,II)
	else
		J=WorkOnJ_Sparse(A,x,y,I,II,JJ,ISml,JSml,VSml,N)	
	end
	#If you want to compare the outputs of the two functions above:
	#=
	(I1,J1,V1)=findnz(J);
	(I2,J2,V2)=findnz(J2);
	if !isequal(V1,V2)
	   @info size(V1) size(V2)
	   error("Not eq")
	end	
	=#
	
	if min(size(A,1),50)/2<30
		restart=min(size(A,1),50)/2;  
	else
		restart = 30; 
	end  
	restart	= convert(Int64,restart)
	
	for i=1:N
		dx[i]=0.0 #Reset Newton direction
		phiM[i]=-phi[i];	
	end

	#Adding to preexisting mat first 1:n rows (memoryless) then doing a view. 
	#Krylov is faster if we view ordered subsection of the matrix
	icount=0
	jcount=0
	n=0
	for i=1:N
		if I[i]==true
			n+=1
		end
	end
	for i=1:length(I)
		if I[i]==true
			icount+=1
		end
		for j=1:length(I)
			if I[j]==true
				jcount+=1
			end
			if I[i]==true && I[j]==true
				Jsubs[icount,jcount]=J[i,j];
			end
		end
		jcount=0; #reset
	end

	dxSubset=view(dx,I);
	phiMSubset=view(phiM,I)
	JSubset=view(Jsubs,1:n,1:n);

	###1 IterativeSolvers
	#dxSubset=IterativeSolvers.gmres!(dxSubset,JSubset,phiMSubset,tol=1e-6,restart=restart, initially_zero=true,maxiter=10*restart);
	
	####2 KrylovKit
	#alg = GMRES( krylovdim = restart, maxiter = 5, tol = 1e-6)
	#dxSubset, info = @inferred linsolve(JSubset,phiMSubset,dxSubset, alg)
	
	###3 Krylov.jl.git#gmres
	##DEFINE OUT OF LOOP global dxSubset
	#for i = 1 : restart	
	#	(dxSubset,stats) = Krylov.gmres(JSubset, phiMSubset, rtol=1e-6, x0=dxSubset, itmax=10)
	#end
	
	###4 Krylov DqGmres
	##Some parameters atol::Float64=1.0e-8 rtol::Float64=1.0e-6 itmax::Int=0
	dqgmres_tol = 1.0e-6
	(dxSubset,stats) = Krylov.dqgmres(JSubset, phiMSubset,memory=restart,itmax =10*restart,rtol =dqgmres_tol)

	
	dx[I]=dxSubset;


	# Test if the search direction is smaller than numerical precision. 
	# That is if it is too close to zero.
	for i=1:N
		absdx[i]=abs(dx[i])
	end
	if maximum(absdx) < eps()
		flag = 5;
		# Rather than just giving up we may just use the gradient direction
		# instead. However, I am lazy here!
		#  dx = nabla_phi'
		break;
	end
	
	# Test if we have dropped into a local minimia if so we are stuck
	#nabla_phi = phiT*J;
	mul!(nabla_phi,phiT,J) 		
	if norm(nabla_phi) < tol_abs
		flag = 6;
		break;
	end
	
	# Test if our search direction is a 'sufficient' descent direction
	mul!(nabdx,nabla_phi,dx)
	if  nabdx[1]  > -rho*(dx'*dx)
		flag = 7;
		# Rather than just giving up we may just use the gradient direction
		# instead. However, I am lazy here!
		#  dx = nabla_phi'
		break;
	end
		
	#--- Armijo backtracking combined with a projected line-search ---------
	tau= 1.0;                  # Current step length
	f_0     = err[1];
	grad_f= beta*dot(nabla_phi,dx);


	while true #Inf loop, escapes when a break is performed
		
		for i=1:N
			xdxtau[i]=x[i]+dx[i]*tau
			#non negativity
			if xdxtau[i]<=0.
				x_k[i]=0.
			else
				x_k[i]=xdxtau[i]
			end
		end

		#y_k   = A*x_k + b; (pre allocated y_k)
		mul!(y_k,A,x_k) 
		for i=1:N
				y_k[i]+=b[i]
		end
		phi_k=phi_lambda!( y_k, x_k, lambda,phi_k,phi_l );		
		
		for i=1:length(phi_k); phi_kT[i]=phi_k[i]; end #transpose
		f_k= 0.5*dot(phi_kT,phi_k);       # Natural merit function
		
		# Perform Armijo codition to see if we got a sufficient decrease
		test=f_0+ tau*grad_f;
		if ( f_k <= test)
			break;
		end
		# Test if time-step became too small
		if tau*tau < gamma	
			break;
		end	
		
		tau= alpha*tau;
		
	end #end of while lp 1. 
	
	# Update iterate with result from Armijo backtracking
	for i=1:N
		x[i]=x_k[i]
	end
	# Increment the number of iterations
	iter = iter + 1;	
	
end #end of while lp. 

	
if iter >= max_iter
	flag = 8;
	iter = iter - 1;
end		

return(x); 

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


function f(dx,I) 
   dxSubset=view(dx,I);
   return dxSubset
end