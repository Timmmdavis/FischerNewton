function WorkOnJ(J,A::Array{Float64,2},x::Array{Float64},y::Array{Float64},I)
	
	Indx=findall(I)
	fill!(J, 0.0)
	#println("Resetting for saftey")
	#J=zeros(size(J));
	
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
	J=SparseArrays.sparse(J); 
	return(J)
	
end