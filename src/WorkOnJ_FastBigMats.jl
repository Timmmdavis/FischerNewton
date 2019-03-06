function WorkOnJ_FastBigMats(A,x,y,Flag,II,JJ,ISml,JSml,VSml) #,ISml,JSml,VSml
#A - inf mat
#xy - vectors
#Flag - bool vectors
#II JJ - matricies that mat sparse locations
#ISml etc - vectors to fill, push! on these is as fast...

	Indx=findall(Flag)
	
	I=view(II,Indx,Indx);
	J=view(JJ,Indx,Indx);
	
	#If we dont pre allocate these
	#VSml = Float64[]
	#ISml = Float64[]
	#JSml = Float64[]
	#sizehint!(VSml, length(Indx)^2)
	#sizehint!(ISml, length(Indx)^2)
	#sizehint!(JSml, length(Indx)^2)
	
	
	#Bits of A that are bad
	#Reset J to zeros	
	counter=0;
	for i=1:length(Indx)
		local idxi=Indx[i];
		cons=((y[idxi]^2) +(x[idxi]^2))^0.5;
		p=x[idxi]/cons-1.0;
		q=y[idxi]/cons-1.0;
		for j=1:length(Indx)
			local idxj=Indx[j];		
			if idxi==idxj
				if p==0 && A[idxi,idxj]==0
					continue
				end
				##If we dont pre allocate these
				#push!(VSml,(A[idxi,idxj]*q)+p)
				#push!(ISml,idxi)
				#push!(JSml,idxj)
				counter+=1;				
				VSml[counter]=A[idxi,idxj]*q+p;
				ISml[counter]=idxi;
				JSml[counter]=idxj;
				continue
			end		
			if I[i,j]==0;
				#A way of working out how to find these without the loop over every entry would be far superior. 
				continue
			end
			##If we dont pre allocate these			
			#push!(VSml,A[idxi,idxj]*q)
			#push!(ISml,idxi)
			#push!(JSml,idxj)
			counter+=1;			
			VSml[counter]=A[idxi,idxj]*q;
			ISml[counter]=idxi;
			JSml[counter]=idxj;			
		end

	end	
	#Getting the correct parts...
	Iview=view(ISml,1:counter);
	Jview=view(JSml,1:counter);
	Vview=view(VSml,1:counter);
	
	Jsp=sparse(Iview, JSml[1:counter], Vview,length(x),length(x));
	return(Jsp)
end	