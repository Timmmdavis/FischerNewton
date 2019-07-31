function ResetInitArrays(Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys);

#Resetting if we have already called fischer newton
Ints.iter= 1
Ints.flag= 1;

#Flts.lambda=1.00

for i=1:Ints.N

	Vects.y[i] =0.;
	Vects.phi[i]=0.
	Vects.phi_k[i]=0.;
	Vects.phi_l[i]=0.;
	Vects.phiM[i]=0.;
	Vects.dx[i]=0.;
	Vects.absdx[i]=0.;
	Vects.y_k[i]=0.;
	Vects.xdxtau[i]=0.;
	Vects.x_k[i]=0.;
	Vects.x[i]=0.;
	Bls.I[i]=false

	Mats.phiT[1,i]=0.
	Mats.phi_kT[1,i]=0.
	Mats.nabla_phi[1,i]=0.

end

Arrys.err[1]=Inf
Arrys.old_err[1]=Arrys.err
Arrys.nabdx[1]=0.;
Arrys.test[1]=0.0;
Arrys.grad_f[1]=0.0
Arrys.f_k[1]=0.0

fill!(Mats.J, 0.0)
fill!(Mats.JSubs, 0.0)

if Ints.useSparse==1
	for i=1:Ints.N^2
		Mats.ISml[i]=0.
		Mats.JSml[i]=0.
		Mats.VSml[i]=0.
	end
end

return Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys

end