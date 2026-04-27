% Sets default value for the field_str in the struct s
% s: 		 	struct
% field_str: 	string
% field_value:	default value of field
function s = defaultField(s, field_str, field_value)

if isfield(s, field_str)
	if isempty( getfield(s, field_str) )
		s = setfield(s, field_str, field_value);
	end
else
	s = setfield(s, field_str, field_value);
end



