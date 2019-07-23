function WorkOnJ(J,A::Array{Float64,2},x::Array{Float64},y::Array{Float64},I,II)

	#J defined as: J = zeros(N,N); #spzeros 
	fill!(J, 0.0)
	
	for i=1:length(I); #Threads.@threads
		if I[i]==false 
			continue
		end

		cons=(((y[i]^2) +(x[i]^2))^0.5);
		p=(x[i]./cons)-1.0;
		q=(y[i]./cons)-1.0;

		for j=1:length(I); #Threads.@threads
			if I[j]==false 
				continue
			end	
			if II[i,j]==0;
				#A way of working out how to find these without the loop over every entry would be far superior. 
				continue
			end
			J[i,j]=A[i,j]*q;
		end
		J[i,i] +=p; #can use J[CartesianIndex(idxi,idxi)]
	end
	#J=SparseArrays.sparse(J); 
	return(J)
	

end