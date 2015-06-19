
# `object` and `new_class` should have the same classes
changeClassName = function(object, new_class) {
	new_object = new(new_class)
	for(sn in slotNames(object)) {
		slot(new_object, sn) = slot(object, sn)
	}
	return(new_object)
}


compare_unit = function(u1, u2) {
	x1 = convertUnit(u1, "mm", valueOnly = TRUE)
	x2 = convertUnit(u2, "mm", valueOnly = TRUE)
	ifelse(x1 > x2, 1, ifelse(x1 < x2, -1, 0))
}
 
