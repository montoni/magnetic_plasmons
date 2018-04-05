def particle_placement(rij,numPart):
	Loc = []
	ex = rij*math.sqrt(3)/2
	ey = rij/2.
	'''Place first six particles around [0,0]'''
	theta = 2*math.pi/6
	beta = math.pi/6
	for num in range(6):
		Loc.append([rij*math.cos(theta*num + beta),rij*math.sin(theta*num + beta)])

	'''Do the next few sets of four particles extending to the right'''
	'''Need to define a new center for each set of four'''

	
