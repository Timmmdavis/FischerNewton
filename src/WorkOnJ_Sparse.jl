function WorkOnJ_Sparse(A,x,y,Flag,II,JJ,ISml,JSml,VSml,N) #,ISml,JSml,VSml
#A - inf mat
#xy - vectors
#Flag - bool vectors
#II JJ - matricies that mat sparse locations
#ISml etc - vectors to fill, push! on these is as fast...
	
	counter=0;
	for i=1:length(Flag); 
		if Flag[i]==false 
			continue
		end
		cons=((y[i]^2) +(x[i]^2))^0.5;
		p=x[i]/cons-1.0;
		q=y[i]/cons-1.0;
		for j=1:length(Flag); 
			if Flag[j]==false 
				continue
			end	
			if i==j
				if p==0 && A[i,j]==0
					continue
				end
				counter+=1;				
				VSml[counter]=A[i,j]*q+p;
				ISml[counter]=i;
				JSml[counter]=j;
				continue
			end		
			if A[i,j]==0;
				#A way of working out how to find these without the loop over every entry would be far superior. 
				continue
			end
			counter+=1;			
			VSml[counter]=A[i,j]*q;
			ISml[counter]=i;
			JSml[counter]=j;			
		end

	end	
	#Getting the correct parts...
	Iview=view(ISml,1:counter);
	Jview=view(JSml,1:counter);
	Vview=view(VSml,1:counter);
	
	Jsp=sparse(Iview, Jview, Vview,N,N);
	return(Jsp)
end	