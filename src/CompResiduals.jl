function CompResiduals(a,b)
#Sum of residuals
Res=0.;
for i=1:length(a)
	Res=Res+abs(a[i]-b[i]);
end #Func end

return Res
end

