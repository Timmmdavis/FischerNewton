function string_as_varname(s::String,v::Any)
	s=Symbol(s)
	eval(:($s=$v))
end