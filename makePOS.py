import numpy as np

def make_cell(ln,xw=1,zw=1,xedges=0,yedges=0):
	#add doctests and docstring
	"""A function which stacks monolayers into a simulation cell of correct size. note xw and zw are the number of unit cells width, if edges=True/or if edge types are given??? then terminate them and add padding to the cell
	"""


	#THIS IS WRONG! YOU NEED TO EDIT THIS SO IT WORKS IT OUT IN ABSOLUTE TERMS AND THEN CHANGES EVERYTHING TO BE RELATIVE POSITION, WHICH IS WHAT POSCAR IS. OR JUST USE THE ABSOLUTE ONES IN YOUR POSCAR FILES!


	#tests
	assert type(ln)==int, 'please insert an integer number of layers'
	assert ln>0, 'please insert a positive number of layers'


	atom_type="P"
	atom_number=4


	#define unit cell ranges here!! Remember to use padding! all defined from converged converged bulk pe-d2 vdw!
	axes=3 #i.e.  a,b,c
	a=3.334611
	b=(0.602094-0.397906)*10.476599+10
	c=4.394790

	#add an input to the poscar for monolayer and pull this out by default. #defined from converged bulk pe-d2 vdw
	x1=0
	x2=0.500000*a
	y1=10 
	y2=b
	z1=0.918361*c+10*bool(zw-1)
	z2=0.581639*c+10*bool(zw-1)
	z3=0.081639*c+10*bool(zw-1)
	z4=0.418361*c+10*bool(zw-1)

	#offset of layers #This is defined from converged bulk pe-d2 vdw
	xoffset=-0.500000*3.334611
	yoffset=(0.897906-0.397906)*10.476599 
	ygap=(0.897906-0.60209)*10.476599
	zoffset=0

	#build the monolayer array (with single occupancy)
	atoms_per_cell=4
	mono=np.array([[x1,y1,z1],[x2,y1,z2],[x1,y2,z3],[x2,y2,z4]])



	#FIX THIS. Its a bit of a long way around...
	def atom_stacker(mono,ln):
		#add tests and doctests. Clean this function up its a mess

		multiplier=np.repeat(np.array([[xoffset, yoffset, zoffset]]),atoms_per_cell,axis=0) #This is the shift in a monolayer upon stacking, e.g. a bilayer is one monolayer plus one monolayer offset by these coordinates. We do it 4 times to account for each of the four p atoms in unit cell

		#flatten to make each row of the array later a monolayer in the simulation cell
		mono=mono.flatten()
		multiplier=multiplier.flatten()
		layers=np.arange(ln).reshape(ln,1) #This makes a vector from zero to ln, will multiply the offset vector with this later
		multiplier=np.multiply(layers,multiplier)


		#Make the multilayer
		layers_offset=multiplier.flatten()
		monobase=np.tile(mono,ln) #Make an array with each row a monolayer of coordinates for the four atoms, ln layers thick
		atoms=monobase+layers_offset


		#reshape atoms to be only axes number of columns
		atoms=np.reshape(atoms,(atoms_per_cell*ln,axes)) #reshape(total numer of rows=number of layers*number of atoms per layer,no of coordinates i.e. x,y,z=3)
		layers[0]=1
		layers=np.repeat(layers,4,axis=0)
		mod_x_wrapped=np.floor((np.absolute(atoms[:,0,np.newaxis])/xoffset)%2)*np.absolute(xoffset)
		atoms[:,0]=mod_x_wrapped[:,0]
		return atoms


	#make the cell the right size
	def cell_expander(a,b,c,ln,xw=1,zw=1,xedges=0,zedges=0):
		#add tests and doctests
		y_height=10+10+ygap*(ln-1)+(b-10)*ln  #bottom padding + top padding+ total gaps between layers + total thickness of layers
		unit_cell=np.array([[a*xw+20*bool(xw-1),0,0],[0,y_height,0],[0,0,c*zw+20*bool(zw-1)]])

		return unit_cell


	#function to build the poscar string
	def build_poscar(atoms,unit_cell,ln,atom_type, atom_number,combo_rib_dim):
		#do doctests etc

		#CHECK these indentations work/are necessary. Make them do new lines!
		start_string=(str(atom_type)+"1"+'\n'+"1.0"+'\n')
		mid_string=('\n'+'\t'+str(atom_type)+'\n'+'\t'+
			str(atom_number*ln*combo_rib_dim)+'\n'+
			"Cartesian"+'\n')
		spaces=" "*9

		#turn the numpy arrays for unit cell and positions into strings
		unit_cell_string='\n\t\t'.join(spaces.join('%0.10f' %x for x in y) for y in unit_cell)
		cartesian_string='\n\t'.join(spaces.join('%0.10f' %x for x in y) for y in atoms)

		poscar_string=start_string+'\t\t'+unit_cell_string+mid_string+'\t'+cartesian_string

		return poscar_string


	#now need to print POSCAR out
	def write_poscar(poscar_string,filename):
		#add doctests etc

		assert type(filename)==str, "please enter a string as filename"

		filename=filename+".vasp"
		with open(filename,'w') as vasp_file:
			vasp_file.write(poscar_string)

		print("writing completed")
		return True



	def ribbon_maker(atoms,unit_cell,xw=1,zw=1):
		#make edges
		"""This function takes in an array of atom positions one until cell in the slab direction, therefore infinite direction, with ln stacks
		The function outputs an array of atom positions xw and/or zw unit cells wide, therefore making non-inifinite slabs.
		Note:the atom positions are shifted in the positive direction by 10A to account for padding, only if xw or zw !=0"""


		#add a duplicate of the atoms array with atoms in xw and zw shifted by an integer multiple of unit cell size for each duplicate
		no_atoms_stacked=atoms.shape[0]

		base=np.array([[]])
		for i in range(xw):
			j=0
			basex=np.array([[i,0,j]])
			try:
				base=np.append(base,basex,axis=0)
			except:
				base=np.array([[0,0,0]])
			for j in range(zw-1):
				basez=np.array([[i,0,j+1]])
				base=np.append(base,basez,axis=0)

		####this works okay

		base=np.repeat(base,no_atoms_stacked,axis=0) # this is now a valid multiplier with every combination of unit cell positions in the range you're interested in

		#now multiply each no_atoms_stacked set in base element wise by the atoms array. so e.g. the second set will be 1 0 0 added to the original so shift by one unit cell in x
		combos=zw*xw
####works okay to here
		atoms[:,0]=atoms[:,0]+10*bool(xw-1)#note we only do this for x because of the weird resetting method used in stacker to get the layer sliding, that resets the change that we
		#									would make at the start when we made the change for z
		atoms=np.tile(atoms,(combos,1))#now atoms is the same dimensions as base
		#now add base*unit cell dimensions to atoms
		xshift=(unit_cell[0,0]-20)/xw
		zshift=(unit_cell[2,2]-20)/zw
		atoms=atoms+base*np.array([[xshift,0,zshift]])

		#atoms[:,0]=atoms[:,0]+10*bool(xw-1)
		#atoms[:,2]=atoms[:,2]+10*bool(zw-1)
		#this is now all the valid positions


		#then figure out what to do with the edges

		return atoms



	def edge_modifier(atoms,xedge=0,zedge=0)
		#function to modify the edges of the ribbons
		#there are many types of edges so look into how they're defined and any common method before you start

		
		#symmetric edges?
		#orientation related types -> classes?

		#H-termination?


		return False


	#now call functions
	atoms=atom_stacker(mono,ln)
	unit_cell=cell_expander(a,b,c,ln,xw,zw)
	atoms=ribbon_maker(atoms,unit_cell,xw,zw)
	#x starts off shifted correctly, but then the next layer is shifted negatively instead of positively


	combo_rib_dim=xw*zw
	poscar_string=build_poscar(atoms,unit_cell,ln,atom_type, atom_number,combo_rib_dim)
	write_to_file=write_poscar(poscar_string,'test')

	return atoms


#this is notpart of any function
#NOTE! you need to change this so that user can input these parameters make_cell(ln,xw=1,zw=1,xedges=0,yedges=0)
ln=3
xw=1
zw=5

#you might be able to use symmetry to make this all alot easier/faster??
atoms=make_cell(ln,xw,zw)