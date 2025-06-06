# Do a raytrace over positions on the focal plane from y=0->1
n = NSUR() # Number of surfaces

PRINT "X position (mm) | Angle (deg)"
FOR fpinc, 0, 215, 10.75
	RAYTRACE 0,fpinc/215,0,0
	PRINT fpinc, ",",ACOS(RAYN(n))*180.0/3.14159265359
NEXT

