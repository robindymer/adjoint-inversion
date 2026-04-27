% Helper function that saves space in the input files
function domain = domainStruct(domain, def, bc, faultBoundaryGroups)
	domain.def = def;
	domain.bc = bc;
	domain.faultBoundaryGroups = faultBoundaryGroups;
	domain.surfaceBoundaryGroups = [];
end