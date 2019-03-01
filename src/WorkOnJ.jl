function WorkOnJ(J,A::Array{Float64,2},x::Array{Float64},y::Array{Float64},Indx)



	# for i=1:length(Indx)
		# local idxi=Indx[i];
		# cons=(((y[idxi]^2) +(x[idxi]^2))^0.5);
		# p[idxi]=(x[idxi]./cons)-1.0;
		# q[idxi]=(y[idxi]./cons)-1.0;
		# for j=1:length(Indx)
			# local idxj=Indx[j];
			# J[idxi,idxj]=A[idxi,idxj].*q[idxi];
		# end
	# end
	# # #and add p to diagonal
	# # for i=1:length(Indx)
		# # local idxi=Indx[i];
		# # J[idxi,idxi] +=p[i];
	# # end
	
	# #J[Indx,Indx].=A[Indx,Indx].*q[Indx];
	# J[Indx,Indx]=J[Indx,Indx]+Diagonal(p[Indx])
	
	for i=eachindex(Indx)
		local idxi=Indx[i];
		cons=(((y[idxi]^2) +(x[idxi]^2))^0.5);
		p=(x[idxi]./cons)-1.0;
		q=(y[idxi]./cons)-1.0;
		for j=eachindex(Indx)
			local idxj=Indx[j];
			J[idxi,idxj]=A[idxi,idxj].*q;
		end
		J[idxi,idxi] +=p; #can use J[CartesianIndex(idxi,idxi)]
	end

end