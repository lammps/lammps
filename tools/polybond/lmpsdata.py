#
#
#	lmpsdata.py
#
#	For reading, writing and manipulating lammps data files
# 	For calculation of certain properties using lammps data files
# 	For creating VMD input text files using lammps data files
#	All x,y,z calculations assume the information includes image flags

class Lmpsdata:
	def __init__(self,file,atomtype):
		"""initiates lammps data structures"""
		self.atomtype=atomtype		
		self.keywords=[]
		self.atoms=[]
		self.angles=[]
		self.bonds=[]
		self.dihedrals=[]
		self.dipoles=[]
		self.impropers=[]
		self.masses=[]
		self.shapes=[]
		self.velocities=[]
		self.anglecoef=[]
		self.bondcoef=[]
		self.dihedralcoef=[]
		self.impropercoef=[]
		self.paircoef=[]
		self.read(file)

	def read(self,file):
		"""Reads in lammps data file
		Skips the first line of the file (Comment line)
		First reads the header portion of the file (sectflag=1)
			blank lines are skipped
			header keywords delineate an assignment
			if no header keyword is found, body portion begins
		Second reads the body portion of the file (sectflag=2)
			first line of a section has only a keyword
			next line is skipped 
			remaining lines contain values
			a blank line signifies the end of that section
			if no value is listed on a line than a keyword must be used
		File is read until the end
		The keywords read in are stored in self.keywords"""
		if file=='':
			print 'no file is given. Will have to build keywords and class structures manually'
			return
		sectflag=1
		f=open(file,'r')
		f.readline()
		for line in f:
			row=line.split()
			if sectflag==1:
				if len(row)==0:
					#skip line the line is blank
					pass
				elif len(row)==1:
					# Set Sectflag to 2, assume line is a keyword
					sectflag=2
					checkkey=1 #ensures keyword will be checked in body portion
					keyword=row 
				elif len(row)==2:
					if row[1]=='atoms':
						self.atomnum=row[0]
						self.keywords.append('atoms')
					elif row[1]=='bonds':
						self.bondnum=row[0]
						self.keywords.append('bonds')
					elif row[1]=='angles':
						self.anglenum=row[0]
						self.keywords.append('angles')
					elif row[1]=='dihedrals':
						self.dihedralnum=row[0]
						self.keywords.append('dihedrals')
					elif row[1]=='impropers':
						self.impropernum=row[0]
						self.keywords.append('impropers')
					else:
						# Set Sectflag to 2, assume line is a keyword
						sectflag=2
						checkkey=1 #ensures keyword will be checked in body portion
						keyword=row					
				elif len(row)==3:
					if row[1]=='atom' and row[2]=='types':
						self.atomtypenum=row[0]
						self.keywords.append('atom types')
					elif row[1]=='bond' and row[2]=='types':
						self.bondtypenum=row[0]
						self.keywords.append('bond types')
					elif row[1]=='angle' and row[2]=='types':
						self.angletypenum=row[0]
						self.keywords.append('angle types')
					elif row[1]=='dihedral' and row[2]=='types':
						self.dihedraltypenum=row[0]
						self.keywords.append('dihedral types')
					elif row[1]=='improper' and row[2]=='types':
						self.impropertypenum=row[0]
						self.keywords.append('improper types')
					else:
						# Set Sectflag to 2, assume line is a keyword
						sectflag=2
						checkkey=1 #ensures keyword will be checked in body portion
						keyword=row
				elif len(row)==4:
					if row[2]=='xlo' and row[3]=='xhi':
						self.xdimension=[row[0], row[1]]
						self.keywords.append('xlo xhi')
					elif row[2]=='ylo' and row[3]=='yhi':
						self.ydimension=[row[0], row[1]]
						self.keywords.append('ylo yhi')
					elif row[2]=='zlo' and row[3]=='zhi':
						self.zdimension=[row[0], row[1]]
						self.keywords.append('zlo zhi')
					else:
						# Set Sectflag to 2, assume line is a keyword
						sectflag=2
						checkkey=1 #ensures keyword will be checked in body portion
						keyword=row
				elif len(row)==5:
					if row[1]=='extra' and row[2]=='bond' and row[3]=='per' and row[4]=='atom':
						self.extrabonds=row[0]
						self.keywords.append('extra bond per atom')
					else:
						# Set Sectflag to 2, assume line is a keyword
						sectflag=2
						checkkey=1 #ensures keyword will be checked in body portion
						keyword=row
				elif len(row)==6:
					if row[3]=='xy' and row[4]=='xz' and row[5]=='yz':
						self.tilt=[row[0], row[1], row[2]]
						self.keywords.append('xy xz yz')
					else:
						# Set Sectflag to 2, assume line is a keyword
						sectflag=2
						checkkey=1 #ensures keyword will be checked in body portion
						keyword=row
				else:
					# set sectflag to 2, assume line is a keyword
					sectflag=2
					checkkey=1 #ensures keyword will be checked in body portion
					keyword=row				
			elif sectflag==2:
				if checkkey==1:
					if len(keyword)==1:
						if keyword[0]=='Atoms' or keyword[0]=='Velocities' or keyword[0]=='Masses' or\
						keyword[0]=='Shapes' or keyword[0]=='Dipoles' or keyword[0]=='Bonds' or\
						keyword[0]=='Angles' or keyword[0]=='Dihedrals' or keyword[0]=='Impropers':
							bodyflag=1
							blanknum=0							
							self.keywords.append(keyword[0])
							checkkey=0
						else:
							bodyflag=0
							checkkey=0
					elif len(keyword)==2:
						if row[1]=='Coeffs' and (row[0]=='Pair' or row[0]=='Bond' or row[0]=='Angle' or\
						row[0]=='Dihedral' or row[0]=='Improper'):#class 2 force field keywords not included
							bodyflag=1
							blanknum=0							
							self.keywords.append('{0} {1}'.format(keyword[0],keyword[1]))
							checkkey=0
						else:
							bodyflag=0
							checkkey=0
					else:
						#egnore line and set bodyflag=0
						bodyflag=0
						checkkey=0					
				if bodyflag==0: #bodyflag 0 means no body keyword has been found
					if len(row)==1:
						if row[0]=='Atoms' or row[0]=='Velocities' or row[0]=='Masses' or\
						row[0]=='Shapes' or row[0]=='Dipoles' or row[0]=='Bonds' or\
						row[0]=='Angles' or row[0]=='Dihedrals' or row[0]=='Impropers': 
							bodyflag=1
							blanknum=0
							keyword=row
							self.keywords.append(keyword[0])
						else:
							#egnore line
							pass
					elif len(row)==2:
						if row[1]=='Coeffs' and (row[0]=='Pair' or row[0]=='Bond' or row[0]=='Angle' or\
						row[0]=='Dihedral' or row[0]=='Improper'):#class 2 force field keywords not included	
							bodyflag=1
							blanknum=0
							keyword=row
							self.keywords.append('{0} {1}'.format(keyword[0],keyword[1]))
						else:
							#egnore line
							pass
					else:
						#egnore line
						pass
				elif bodyflag==1: #currently assumes 1 or more blank lines are between body data keywords
					if len(row)==0:
						blanknum+=1
						if blanknum>1:
							bodyflag=0
					elif len(keyword)==1:
						if keyword[0]=='Atoms':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.atoms.append(row)
						elif keyword[0]=='Velocities':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.velocities.append(row)
						elif keyword[0]=='Masses':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.masses.append(row)
						elif keyword[0]=='Shapes':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.shapes.append(row)
						elif keyword[0]=='Dipoles':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.dipoles.append(row)
						elif keyword[0]=='Bonds':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.bonds.append(row)
						elif keyword[0]=='Angles':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.angles.append(row)
						elif keyword[0]=='Dihedrals':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.dihedrals.append(row)
						elif keyword[0]=='Impropers':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.impropers.append(row)
						else:
							#egnore line and change bodyflag to 0
							bodyflag=0
					elif len(keyword)==2:
						if keyword[0]=='Pair' and keyword[1]=='Coeffs':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.paircoef.append(row)
						elif keyword[0]=='Bond' and keyword[1]=='Coeffs':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.bondcoef.append(row)
						elif keyword[0]=='Angle' and keyword[1]=='Coeffs':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.anglecoef.append(row)
						elif keyword[0]=='Dihedral' and keyword[1]=='Coeffs':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.dihedralcoef.append(row)
						elif keyword[0]=='Improper' and keyword[1]=='Coeffs':
							try:
								int(row[0])
							except ValueError:
								 keyword=row
								 checkkey=1
							else:
								self.impropercoef.append(row)
						else:
							#egnore line and change bodyflag to 0
							bodyflag=0
					else:
						#egnore line and change bodyflag to 0
						bodyflag=0
		f.close()
										
	def write(self,file,modflag):
		"""Write lammps data files using the lammps keywords and lammpsdata structures
		writes first line of the file as a Comment line
		if no modifications to any of the lammpsdata structures (modflag=0)
			Use Keywords to write lammpsdata structures directly 
		if modifications to any of the lammpsdata structures (modflag=1)
			Key Lammpsdata structures like atom numbers, coefficient numbers
			need to be modified to match the other modified lammpsdata structures
			This section will still use the keywords to write lammpsdata structures.
		For all modflags, the code will write data to the file until all of the
		keyword's data structures have been finished writing. 
		The keywords are stored in self.keywords"""
		f=open(file,'w')
		f.write('polymer data file\n')
		for row in self.keywords:
			if row=='atoms':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.atomnum,row)) 
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.atoms),row))
			elif row=='bonds':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.bondnum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.bonds),row))
			elif row=='angles':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.anglenum,row)) 
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.angles),row)) 
			elif row=='dihedrals':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.dihedralnum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.dihedrals),row))
			elif row=='impropers':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.impropernum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.impropers),row))
			elif row=='atom types':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.atomtypenum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.masses),row))
			elif row=='bond types':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.bondtypenum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.bondcoef),row))
			elif row=='angle types':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.angletypenum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.anglecoef),row))
			elif row=='dihedral types':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.dihedraltypenum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.dihedralcoef),row))
			elif row=='improper types':
				if modflag==0:
					f.write('{0} {1}\n'.format(self.impropertypenum,row))
				elif modflag==1:
					f.write('{0} {1}\n'.format(len(self.impropercoef),row))
			elif row=='xlo xhi':
				f.write('{0} {1} {2}\n'.format(self.xdimension[0],self.xdimension[1],row))
			elif row=='ylo yhi':
				f.write('{0} {1} {2}\n'.format(self.ydimension[0],self.ydimension[1],row))
			elif row=='zlo zhi':
				f.write('{0} {1} {2}\n'.format(self.zdimension[0],self.zdimension[1],row))
			elif row=='extra bond per atom':
				f.write('{0} {1}\n'.format(self.extrabonds,row))
			elif row=='xy xz yz':
				f.write('{0} {1} {2} {3}\n'.format(self.tilt[0],self.tilt[1],self.tilt[2],row))
			elif row=='Atoms':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.atoms:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Velocities':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.velocities:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Masses':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.masses:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Shapes':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.shapes:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Dipoles':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.dipoles:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Bonds':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.bonds:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Angles':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.angles:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Dihedrals':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.dihedrals:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
			elif row=='Impropers':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.impropers:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Pair Coeffs':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.paircoef:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Bond Coeffs':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.bondcoef:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Angle Coeffs':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.anglecoef:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Dihedral Coeffs':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.dihedralcoef:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			elif row=='Improper Coeffs':
				f.write('\n{0}'.format(row)) #new line between header and body portion or two body keywords
				f.write('\n')  #new line between body keyword and body data
				for line in self.impropercoef:
					f.write('\n') #creates a new line for adding body data
					for item in line:
						f.write(' {0}'.format(item)) #adds in each piece of body data with space imbetween
				f.write('\n') #allows space to be added between the end of body data and a new keyword
			else:
				pass
		f.close()
		
	def atomorder(self):
		"""Takes self.atoms and organizes the atom id from least to greatest.
		If the atom ids are already ordered this algorithm will do nothing."""
		current=range(len(self.atoms[0])) # initialize current [assumes self.atoms coloumn #'s does not change]
		for i in range(1,len(self.atoms)): #when i=0, self.atoms will not change; therefore its skipped 
			for k in range(len(self.atoms[i])):
				current[k]=self.atoms[i][k]
			for j in range(i-1,-1,-1):
				for k in range(len(self.atoms[j])):
					self.atoms[j+1][k]=self.atoms[j][k]
				if int(current[0]) > int(self.atoms[j][0]):
					for k in range(len(current)):
						self.atoms[j+1][k]=current[k]
					break
				elif j==0:
					for k in range(len(current)):
						self.atoms[j][k]=current[k]
					#dont need a break here because this is the last j value in the for loop	
	
	def addatoms(self, atoms, retlist=False):
		"""Appends atoms to self.atoms. Assumes atoms are written in the correct atomtype format.
		Change added atom numbers in self.atoms so self.atoms will be in increasing sequential order.
		If retlist is True return a list of the modified atom numbers otherwise exit."""
		initpos=len(self.atoms) #Store index of first appended atoms
		for item in atoms:
			self.atoms.append(range(len(item))) #initializing spots for the new atoms to go
		numberchange=[] #initiate number change where the changed atom numbers will be stored
		count=0
		for i in range(initpos,len(self.atoms)):
			for j in range(len(self.atoms[i])):
				if j==0:
					self.atoms[i][0]=str(i+1)
					numberchange.append(str(i+1))
				else:
					self.atoms[i][j]=atoms[count][j]
			count+=1
		if retlist: return numberchange
		else: return #need to redo this algorithm so their is no referencing going on.

	def adddata(self,data,keyword):
		"""Adds data to a keyword's existing data structure
		All body keywords are viable except Atoms which has its own method
		All header keywords will do nothing in this algrorithm because they have their own algorithm
		changes the body id number so id numbers will be in sequential order""" 
		if keyword=='Velocities':
			pos=len(self.velocities)
			for item in data:
				self.velocities.append(item) #adds data to existing data structure
				self.velocities[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Masses':
			pos=len(self.masses)
			for item in data:
				self.masses.append(item) #adds data to existing data structure
				self.masses[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Shapes':
			pos=len(self.shapes)
			for item in data:
				self.shapes.append(item) #adds data to existing data structure
				self.shapes[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Dipoles':
			pos=len(self.dipoles)
			for item in data:
				self.dipoles.append(item) #adds data to existing data structure
				self.dipoles[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Bonds':
			pos=len(self.bonds)
			for item in data:
				self.bonds.append(item) #adds data to existing data structure
				self.bonds[pos][0]=str(pos+1)
				pos+=1
		elif keyword=='Angles':
			pos=len(self.angles)
			for item in data:
				self.angles.append(item) #adds data to existing data structure
				self.angles[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Dihedrals':
			pos=len(self.dihedrals)
			for item in data:
				self.dihedrals.append(item) #adds data to existing data structure
				self.dihedrals[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Impropers':
			pos=len(self.impropers)
			for item in data:
				self.impropers.append(item) #adds data to existing data structure
				self.impropers[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Pair Coeffs':
			pos=len(self.paircoef)
			for item in data:
				self.paircoef.append(item) #adds data to existing data structure
				self.paircoef[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Bond Coeffs':
			pos=len(self.bondcoef)
			for item in data:
				self.bondcoef.append(item) #adds data to existing data structure
				self.bondcoef[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Angle Coeffs':
			pos=len(self.anglecoef)
			for item in data:
				self.anglecoef.append(item) #adds data to existing data structure
				self.anglecoef[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Dihedral Coeffs':
			pos=len(self.dihedralcoef)
			for item in data:
				self.dihedralcoef.append(item) #adds data to existing data structure
				self.dihedralcoef[pos][0]=str(pos+1) #changing body id number
				pos+=1
		elif keyword=='Improper Coeffs':
			pos=len(self.impropercoef)
			for item in data:
				self.impropercoef.append(item) #adds data to existing data structure
				self.impropercoef[pos][0]=str(pos+1) #changing body id number
				pos+=1
	
#	def modifiydata(data,keyword): Will not implement 
#	"""for modifying header data"""

	def changeatomnum(self,data,originalatoms,atomchanges,vflag=False):
		"""Takes data which contains atom ids 
		and alters those ids from the original atom id system in originalatoms 
		to the new atom id system in atomchanges."""
		#set up boolean array to match the size of data and with all True values
#		print atomchanges
		array=booleanarray(len(data),len(data[0]),True)
#		print originalatoms
		if vflag: # for changing atomnumbers for velocities
			for i in range(len(originalatoms)): #len of originalatoms should match len of atomchanges
				if originalatoms[i][0]==atomchanges[i]: continue
				for j in range(len(data)):
					for k in range(1):
						if data[j][k]==originalatoms[i][0]: #checks if the atom id in data
						#matches the atom id in origianalatoms
 							if array.getelement(j,k): #if the boolean array is true
								#set the boolean array to False and 
								#change the atom id in data to atomchanges
								array.setelement(j,k,False)
								data[j][k]=atomchanges[i]
								break
								#the change of boolean array to False ensures
								#the atom id in data will only be changed once.
		else: # for changing atom numbers for everything else
			for i in range(len(originalatoms)): #len of originalatoms should match len of atomchanges
				if originalatoms[i][0]==atomchanges[i]: continue
				for j in range(len(data)):
					for k in range(2,len(data[j])):
						if data[j][k]==originalatoms[i][0]: #checks if the atom id in data
						#matches the atom id in origianalatoms
 							if array.getelement(j,k): #if the boolean array is true
								#set the boolean array to False and 
								#change the atom id in data to atomchanges
								array.setelement(j,k,False)
								data[j][k]=atomchanges[i]
								break
								#the change of boolean array to False ensures
								#the atom id in data will only be changed once.
#		print data
		return data

	def deletebodydata(self,keyword):
		"""Changes a keyword's class data structure to an empty structure"""
		if keyword=='Velocities':
			self.velocities=[]
		elif keyword=='Masses':
			self.masses=[]
		elif keyword=='Shapes':
			self.shapes=[]
		elif keyword=='Dipoles':
			self.dipoles=[]
		elif keyword=='Bonds':
			self.bonds=[]
		elif keyword=='Angles':
			self.angles=[]
		elif keyword=='Dihedrals':
			self.dihedrals=[]
		elif keyword=='Impropers':
			self.impropers=[]
		elif keyword=='Pair Coeffs':
			self.paircoef=[]
		elif keyword=='Bond Coeffs':
			self.bondcoef=[]
		elif keyword=='Angle Coeffs':
			self.anglecoef=[]
		elif keyword=='Dihedral Coeffs':
			self.dihedralcoef=[]
		elif keyword=='Improper Coeffs':
			self.impropercoef=[]
		elif keyword=='Atoms':
			self.atoms=[]

	def extractmolecules(self,molecule):
		"""Takes the variable molecule and
		extracts the individual molecules' data back into lmpsdata.
		This extraction takes place through a 4 step process.
		Step 1: Use a molecule's keywords to alter the lmpsdata data structures to empty.
				To accomplish this procedure use the lmpsdata method deletebodydata.
		Step 2: Add the molecules' atoms to lmpsdata's atoms using the method addatoms.
				Return a list of atom id changes for each molecule.
		Step 3: Utilize each molecules list of atom id changes to change their data's atom id numbers.
				Uses changeatomnum and returns the altered molecule's data.
		Step 4: Add the altered molecules' data to lmpsdata's data using the method adddata"""
		#Use molecule index 0's keyword to change the equivalent lmpsdata structures to empty
		print 'extracting the molecules back to data'
		for keyword in molecule[0].keywords: # step 1
			self.deletebodydata(keyword)
		atomchanges=[] #initializing step 2
		for i in range(len(molecule)): # step 2
			atomchanges.append(self.addatoms(molecule[i].atoms,True))
		for i in range(len(molecule)): # step 3
			for keyword in molecule[i].keywords:
				if keyword=='Angles':
					molecule[i].angles=self.changeatomnum(molecule[i].angles,molecule[i].atoms,atomchanges[i])						
				if keyword=='Bonds':
					molecule[i].bonds=self.changeatomnum(molecule[i].bonds,molecule[i].atoms,atomchanges[i])						
				elif keyword=='Dihedrals':
					molecule[i].dihedrals=self.changeatomnum(molecule[i].dihedrals,molecule[i].atoms,atomchanges[i])
				elif keyword=='Velocities':
					molecule[i].velocities=self.changeatomnum(molecule[i].velocities,molecule[i].atoms,atomchanges[i],True)
				elif keyword=='Impropers':
					molecule[i].impropers=self.changeatomnum(molecule[i].impropers,molecule[i].atoms,atomchanges[i])				
		for i in range(len(molecule)): # step 4
			for keyword in molecule[i].keywords:
				if keyword=='Angles':
					self.adddata(molecule[i].angles,keyword)					
				elif keyword=='Bonds':
					self.adddata(molecule[i].bonds,keyword)						
				elif keyword=='Dihedrals':
					self.adddata(molecule[i].dihedrals,keyword)
				elif keyword=='Velocities':
					self.adddata(molecule[i].velocities,keyword)
				elif keyword=='Impropers':
					self.adddata(molecule[i].impropers,keyword)
				
	def density(self,ringsize,init,final,file=' '):
		"""Create spherical shells from the initial radius to the final radius
		The spherical shell's thickness is defined by ringsize
		Calculate the mass of particles within each spherical shell
		Calculate the volume of each spherical shell
		Then calculate the density in each spherical shell
		Current code is written with the assumption the spheres are centered around the origin
		If a file is listed place the results in the variable file otherwise return the results
		The distance calculations are written assuming image flags are in the data file"""
		rho=[]
		r=[init] #initial radius to examine
		radius=init+ringsize 
		while radius<final: #finding the rest of the radii to examine 
			r.append(radius)
			radius+=ringsize			
		for radius in r:
			volume=4.0/3.0*3.1416*((radius+ringsize)**3 - radius**3)
			mass=0.0
			for row in self.atoms:
				if self.atomtype=='angle' or self.atomtype=='atomic' or self.atomtype=='bond' or\
				self.atomtype=='charge' or self.atomtype=='colloid' or self.atomtype=='electron' or\
				self.atomtype=='full' or self.atomtype=='granular' or self.atomtype=='molecular' or\
				self.atomtype=='peri':
					l=len(row)-1
					dist=float(row[l-3])**2+float(row[l-4])**2+float(row[l-5])**2
					if dist<(radius+ringsize)**2 and dist>=radius**2:						
						if self.atomtype=='atomic' or self.atomtype=='charge' or\
						self.atomtype=='colloid' or self.atomtype=='electron' or\
						self.atomtype=='granular' or self.atomtype=='peri':
							type=int(row[1])-1
						else:	
							type=int(row[2])-1
						mass+=float(self.masses[type][1])
						#assumes self.masses is written in atomtype order
				elif self.atomtype=='dipole':
					l=len(row)-1
					dist=float(row[l-6])**2+float(row[l-7])**2+float(row[l-8])**2
					if dist<(radius+ringsize)**2 and dist>=radius**2:
						type=int(row[1])-1
						mass+=float(self.masses[type][1])
						#assumes self.masses is written in atomtype order
				elif self.atomtype=='ellipsoid':
					l=len(row)-1
					dist=float(row[l-7])**2+float(row[l-8])**2+float(row[l-9])**2 
					if dist<(radius+ringsize)**2 and dist>=radius**2:
						type=int(row[1])
						mass+=float(self.masses[type][1])
						#assumes self.masses is written in atomtype order
				elif self.atomtype=='hybrid':
					dist=float(row[2])**2+float(row[3])**2+float(row[4])**2
					if dist<(radius+ringsize)**2 and dist>=radius**2:
						type=int(row[1])
						mass+=float(self.masses[type][1])
						#assumes self.masses is written in atomtype order
			rho.append(mass/volume)			
		if file==' ':
			return r,rho
		else:
			f=open(file,'w')		
			for i in range(len(r)):
				f.write('{0}  {1}\n'.format(r[i],rho[i]))
			f.close()
			return
			
	def createxyz(self,file, routine='mass', values=None):
		"""Two possible routines one to use the masses from data and the other to use the atom type and values supplied by the user.
		The mass version is assessed by setting the routine to 'mass' which is the default method.
		The other version is assessed by setting the routine to 'atomtype'.
		The other version takes values which is a list containing the value the user wants to assign those atomtypes to.
		The atomtypes of the values in the list will start at 1 even if no atoms in molecule use 1.
		This makes it easier to find the atomtype and assign the value from the list
		All atom data is assumed to have image flags in the data."""
		f=open(file,'w')
		f.write('{0}\n'.format(len(self.atoms)))
		f.write('atoms\n')
		if routine=='mass':
			for line in self.atoms:
				if self.atomtype=='angle' or  self.atomtype=='bond' or self.atomtype=='full' or\
				self.atomtype=='molecular':				
					type=int(line[2])-1 #3rd position		
				else:
					type=int(line[1])-1 #2nd position
				mass=self.masses[type][1]
				if mass=='12.0107': elementnum=6
				elif mass=='15.9994': elementnum=8
				elif mass=='26.9815':elementnum=13
				else:
					print 'no matching mass value. The method will exit'
					return
				l=len(line)-1
				if self.atomtype=='angle' or self.atomtype=='atomic' or self.atomtype=='bond' or\
				self.atomtype=='charge' or self.atomtype=='colloid' or self.atomtype=='electron' or\
				self.atomtype=='full' or self.atomtype=='granular' or self.atomtype=='molecular' or\
				self.atomtype=='peri':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-5],line[l-4],line[l-3]))
				elif self.atomtype=='dipole':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-8],line[l-7],line[l-6]))
				elif self.atomtype=='ellipsoid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-9],line[l-8],line[l-7]))
				elif self.atomtype=='hybrid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[2],line[3],line[4]))
			f.close()
		elif routine=='atomtype':
			for line in self.atoms:
				if self.atomtype=='angle' or  self.atomtype=='bond' or self.atomtype=='full' or\
				self.atomtype=='molecular':				
					type=int(line[2])-1 #3rd position		
				else:
					type=int(line[1])-1 #2nd position
				elementnum=values[type]
				l=len(line)-1
				if self.atomtype=='angle' or self.atomtype=='atomic' or self.atomtype=='bond' or\
				self.atomtype=='charge' or self.atomtype=='colloid' or self.atomtype=='electron' or\
				self.atomtype=='full' or self.atomtype=='granular' or self.atomtype=='molecular' or\
				self.atomtype=='peri':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-5],line[l-4],line[l-3]))
				elif self.atomtype=='dipole':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-8],line[l-7],line[l-6]))
				elif self.atomtype=='ellipsoid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-9],line[l-8],line[l-7]))
				elif self.atomtype=='hybrid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[2],line[3],line[4]))		
			f.close()


class booleanarray:
	"""A class that stores boolean values in a list of lists."""	
	def __init__(self,rownum, colnum, initval): 		
		""" initialize a list of lists (array) with
		rownum correspondinig to the number of lists in the list and 
		colnum corresponding to the number of elements in the list's list.
		initval is the value the list of lists will be initialized with.
		initval should be a boolean value."""
		# initializing self.array
		self.array=[]
		for i in range(rownum):
			self.array.append(range(colnum))
		# setting elements of self.array to initval
		for i in range(rownum):
			for j in range(colnum):
				self.setelement(i,j,initval)
	
	def setelement(self,rownum, colnum, value):
		"""Assigns value to the list of lists (array) element at rownum and colnum."""
		self.array[rownum][colnum]=value
	
	def getelement(self,rownum,colnum):
		"""Returns element in the list of lists (array) at rownum and colnum."""
		return self.array[rownum][colnum]

def atomdistance(a,b,atomtype):
	"""Returns the distance between atom a and b.
	Atomtype is the style the atom data is written in.
	All atom data is assumed to have image flags in the data.
	Atom a and atom b are assumed to be the same atomtype."""
	from math import sqrt
	if atomtype=='angle' or atomtype=='atomic' or atomtype=='bond' or\
	atomtype=='charge' or atomtype=='colloid' or atomtype=='electron' or\
	atomtype=='full' or atomtype=='granular' or atomtype=='molecular' or\
	atomtype=='peri':
		l=len(a)-1
		dist=sqrt((float(a[l-3])-float(b[l-3]))**2\
		+(float(a[l-4])-float(b[l-4]))**2\
		+(float(a[l-5])-float(b[l-5]))**2)
	elif atomtype=='dipole':
		l=len(a)-1
		dist=sqrt((float(a[l-6])-float(b[l-6]))**2\
		+(float(a[l-7])-float(b[l-7]))**2\
		+(float(a[l-8])-float(b[l-8]))**2)
	elif atomtype=='ellipsoid':
		l=len(a)-1
		dist=sqrt((float(a[l-7])-float(b[l-7]))**2\
		+(float(a[l-8])-float(b[l-8]))**2\
		+(float(a[l-9])-float(b[l-9]))**2) 
	elif atomtype=='hybrid':
		dist=sqrt((float(a[2])-float(b[2]))**2\
		+(float(a[3])-float(b[3]))**2\
		+(float(a[4])-float(b[4]))**2)
	return dist

def distance(a,coord,atomtype):
	"""Returns the distance between atom a and coord.
	Atomtype is the style the atom data is written in.
	All atom data is assumed to have image flags in the data."""
	from math import sqrt #need to alter this slightly
	if atomtype=='angle' or atomtype=='atomic' or atomtype=='bond' or\
	atomtype=='charge' or atomtype=='colloid' or atomtype=='electron' or\
	atomtype=='full' or atomtype=='granular' or atomtype=='molecular' or\
	atomtype=='peri':
		l=len(a)-1
		dist=sqrt((float(a[l-3])-coord[2])**2\
		+(float(a[l-4])-coord[1])**2\
		+(float(a[l-5])-coord[0])**2)
	elif atomtype=='dipole':
		l=len(a)-1
		dist=sqrt((float(a[l-6])-coord[2])**2\
		+(float(a[l-7])-coord[1])**2\
		+(float(a[l-8])-coord[0])**2)
	elif atomtype=='ellipsoid':
		l=len(a)-1
		dist=sqrt((float(a[l-7])-coord[2])**2\
		+(float(a[l-8])-coord[1])**2\
		+(float(a[l-9])-coord[0])**2) 
	elif atomtype=='hybrid':
		dist=sqrt((float(a[2])-coord[0])**2\
		+(float(a[3])-coord[1])**2\
		+(float(a[4])-coord[0])**2)
	return dist


class particlesurface:
	def __init__(self,particle,cutoff,atomid,atomtype,shape='sphere'):
		"""Builds a particle surface with a specific shape from a particle
		The atoms chosen from the surface will have the specific atomid
		atomid will be given in terms of an integer and not a string."""
		self.particle=particle.atoms
		self.atomtype=atomtype
		self.cutoff=cutoff
		self.surface=[]
		self.particlelocator=[] #stores surface particle locations with respect to the particle
		self.deletedatoms=[]
		if shape=='sphere': self.createspheresurf(atomid)
	
	def createspheresurf(self,atomid):
		"""Assumes sphere is centered around 0,0,0.
		Finds maximum distance particle with atomid.
		Than all atoms within cutoff distance of max distance
		will be included into self.surface
		Also will build a list of where the surface particles are 
		in relationship to the particle.
		Will create a booleanarray to keep track of any changes made to the surface."""		
		particledist=[]		
		for row in self.particle: #Stores all of the distances of the particle in particledist
			if self.atomtype=='angle' or  self.atomtype=='bond' or self.atomtype=='full' or\
			self.atomtype=='molecular':				
				type=int(row[2]) #3rd position		
			else:
				type=int(row[1]) #2nd position
			if type==atomid: #only calculates the distance if the atom has the correct atomid number
				particledist.append(distance(row,[0,0,0],self.atomtype))
			else:
				particledist.append(0.0)
		self.maxdist=max(particledist)# finds the max dist.		
		for i in range(len(particledist)):
			if particledist[i]>=(self.maxdist-self.cutoff):
				self.particlelocator.append(i)
				self.surface.append(self.particle[i])
		self.bool=booleanarray(len(self.surface),1,True)

	def getsurfatom(self,rownum):
		return self.surface[rownum]

	def getbool(self,rownum):
		return self.bool.getelement(rownum,0)

	def setbool(self,rownum,value):
		self.bool.setelement(rownum,0,value)

	def removesurfatom(self,rownum):
		"""note this does not actually remove the surface atom.
		Instead this adds a value to the list deletedatoms.
		This list is later used in extractparticle to actually delete those atoms"""
		self.deletedatoms.append(self.particlelocator[rownum])

	def extractparticle(self):
		""""deletesatoms from particle from the list of deletedatoms 
		and than returns the particle."""
		self.particle=self.deleteatoms(self.particle,self.deletedatoms)		
		return self.particle

	def deleteatoms(self,structure,rows):
		"""delete atoms from particle and shifts the structure down"""
		new=[]
		#multiple copying of b to the rows being replaced.
		if rows==[]:
			for line in structure:
				new.append(line) #if no rows need replacing copy structure
			return new
		for i in range(rows[0]):
			new.append(structure[i]) #copy structure to new until the first replaced row
		count=0
		for i in range(rows[0],len(structure)-len(rows)):# to replace rows and shift undeleted rows over
			for j in range(i+1+count,len(structure)):
				for val in rows:
					if val==j:
						count+=1
						break
				if val==j:continue
				else:
					new.append(structure[j])
					break
		return new
		
	def addatom(self,atomtype,charge=None,moleculenum=None):
		"""Adds an atom to the particle surface between maxdist and the cutoff distance.
		This atom is stored in self.particle and self.surface.
		In the current coding the particle is centered at (0,0,0).
		This method has a required input of atomtype and optional inputs of charge and moleculenum.
		This method currently does not support correctly self.atomtype of dipole, electron, ellipsoid, granular, peri or hybrid
		Will use the image flags 0,0,0.
		Note: The atomtype here is different from the atomtype attribute as part of this class
		Note: If the added atom will be required for bonding later than you will need to extract the particle data
		and rebuild the surface, because this method doesn't include updates to the required variables due to some coding issues"""
		import math, random
		# Checks to make sure self.atomtype is not set to an unsupported value
		if self.atomtype=='dipole' or self.atomtype=='electron' or\
		self.atomtype=='ellipsoid' or self.atomtype=='peri' or self.atomtype=='hybrid':
			print self.atomtype, 'is not currently supported. But, you can add support by adding the needed functionality and inputs'
			return
		
		print 'adding an atom to the surface'
		# Randomly assigns the position of a new atom in a spherical shell between self.cutoff and self.maxdist
		foundpoint=False
		while (foundpoint==False):
			x=random.uniform(-self.maxdist,self.maxdist)
			y=random.uniform(-self.maxdist,self.maxdist)
			try:
				z=random.uniform(math.sqrt(self.cutoff**2-x**2-y**2),math.sqrt(self.maxdist**2-x**2-y**2))
			except ValueError:
				continue
			distance=math.sqrt(x**2+y**2+z**2)
			if distance<self.maxdist and distance>=self.maxdist-self.cutoff:
				foundpoint=True
		
		# initialize new row in self.particle and self.surface
		self.particle.append([])
		self.surface.append([])
		
		pposition=len(self.particle)-1 #particle position
		sposition=len(self.surface)-1 #surface position
		
		# Adds the atom id number to the new row
		self.particle[pposition].append(str(pposition+1))
		self.surface[sposition].append(str(pposition+1))
		
		# Adds the atomtype and the moleculenum if needed to the new row 
		if self.atomtype=='angle' or  self.atomtype=='bond' or self.atomtype=='full' or\
		self.atomtype=='molecular':				
			self.particle[pposition].append(str(moleculenum))
			self.surface[sposition].append(str(moleculenum))
			self.particle[pposition].append(str(atomtype))
			self.surface[sposition].append(str(atomtype))							
		else:
			self.particle[pposition].append(str(atomtype))
			self.surface[sposition].append(str(atomtype))
			
		# Adds the charge if needed to the new row
		if self.atomtype=='charge' or self.atomtype=='full':
			self.particle[pposition].append(str(charge))
			self.surface[sposition].append(str(charge))
			
		# Adds the atom's position to the new row
		self.particle[pposition].append(str(x))
		self.surface[sposition].append(str(x))
		self.particle[pposition].append(str(y))
		self.surface[sposition].append(str(y))
		self.particle[pposition].append(str(z))
		self.surface[sposition].append(str(z))
		
		# Adds the atom's image flags to the new row
		self.particle[pposition].append('0')
		self.surface[sposition].append('0')
		self.particle[pposition].append('0')
		self.surface[sposition].append('0')
		self.particle[pposition].append('0')
		self.surface[sposition].append('0')		
			
		
	def createxyz(self,file,data,routine='mass', values=None):
		"""This shows the particle surface. To show the particle surface after bonding has occurred,
		you will need to extract the particle than reinsert the particle into the class and use createxyz.
		Two possible routines one to use the masses from data and the other to use the atom type and values supplied by the user.
		The mass version is assessed by setting the routine to 'mass' which is the default method.
		The other version is assessed by setting the routine to 'atomtype'.
		The other version takes values which is a list containing the value the user wants to assign those atomtypes to.
		The atomtypes of the values in the list will start at 1 even if no atoms in molecule use 1.
		This makes it easier to find the atomtype and assign the value from the list
		All atom data is assumed to have image flags in the data."""
		f=open(file,'w')
		f.write('{0}\n'.format(len(self.surface)))
		f.write('atoms\n')
		if routine=='mass':
			for line in self.surface:
				if self.atomtype=='angle' or  self.atomtype=='bond' or self.atomtype=='full' or\
				self.atomtype=='molecular':				
					type=int(line[2])-1 #3rd position		
				else:
					type=int(line[1])-1 #2nd position
				mass=data.masses[type][1]
				if mass=='12.0107': elementnum=6
				elif mass=='15.9994': elementnum=8
				elif mass=='26.981539':elementnum=13
				else:
					print 'no matching mass value. The method will exit'
					return 
				l=len(line)-1
				if self.atomtype=='angle' or self.atomtype=='atomic' or self.atomtype=='bond' or\
				self.atomtype=='charge' or self.atomtype=='colloid' or self.atomtype=='electron' or\
				self.atomtype=='full' or self.atomtype=='granular' or self.atomtype=='molecular' or\
				self.atomtype=='peri':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-5],line[l-4],line[l-3]))
				elif self.atomtype=='dipole':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-8],line[l-7],line[l-6]))
				elif self.atomtype=='ellipsoid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-9],line[l-8],line[l-7]))
				elif self.atomtype=='hybrid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[2],line[3],line[4]))
			f.close()
		elif routine=='atomtype':
			for line in self.surface:
				if self.atomtype=='angle' or  self.atomtype=='bond' or self.atomtype=='full' or\
				self.atomtype=='molecular':				
					type=int(line[2])-1 #3rd position		
				else:
					type=int(line[1])-1 #2nd position
				elementnum=values[type]
				l=len(line)-1
				if self.atomtype=='angle' or self.atomtype=='atomic' or self.atomtype=='bond' or\
				self.atomtype=='charge' or self.atomtype=='colloid' or self.atomtype=='electron' or\
				self.atomtype=='full' or self.atomtype=='granular' or self.atomtype=='molecular' or\
				self.atomtype=='peri':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-5],line[l-4],line[l-3]))
				elif self.atomtype=='dipole':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-8],line[l-7],line[l-6]))
				elif self.atomtype=='ellipsoid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-9],line[l-8],line[l-7]))
				elif self.atomtype=='hybrid':
					f.write('{0} {1} {2} {3}\n'.format(elementnum,line[2],line[3],line[4]))		
			f.close()

def molecules(data,init,final,processors, method='all'):
	""" There are two ways to initialize the class Lmpsmolecule.
	Method controls which way is chosen.
	Acceptable values for Method are 'all', 'atom'.
	The default method is 'all'"""
	molecule=[]
	from multiprocessing import Pool
	p=Pool(processors)
	for i in range(init,final+1):
		molecule.append(p.apply_async(Lmpsmolecule,(i,data,method,)))
	for i in range(len(molecule)):
		molecule[i]=molecule[i].get()
	p.close()
	p.join()
	return molecule

class Lmpsmolecule: #Technically should be a meta class but written as a separate class for easier coding.
	def __init__(self,moleculenum,data,method):
		"""initiates lammps molecule structures
		and than extract the appropriate molecular structures from the base class data"""
# ***** Lammps Molecule Structures		
		self.keywords=['Atoms']#include atoms here since this will be in every instance of this class
		self.atoms=[] #extract
		self.angles=[] #extract if keyword is there
		self.bonds=[] #extract if keyword is there
		self.dihedrals=[] #extract if keyword is there
		self.impropers=[] #extract if keyword is there
		self.velocities=[] #extract if keyword is there
# *****	Molecule number
		self.moleculenum=moleculenum
# ***** Extracting the moleculue's atom information from data 
		self.extract(data,method)
# ***** If methods is 'all' then self.extract will extract all of the molecule's information from data
		
	def extract(self,data,method):
		"""uses base class data to extract the molecules atoms 
		and any other molecule data from molecule moleculenum. This extraction is done in a two step process.
		Step 1: extract data.atoms with moleculenum and renumber self.atoms beginning from 1
				store extracted atom numbers from data.atoms in a list called temp
		Step 2: pass temp and data to method changeatomnum 
				which extracts the rest of the Lammps Molecule structures from data
				and changes the atomnumbers in those structures to match the ones already in self.atoms"""
		#checking to make sure atomtype is a valid molecule
		#if not a valid molecule print error and exit method
		if data.atomtype=='full' or data.atomtype=='molecular':
# *****	Sets the molecule's atomtype		
			self.atomtype=data.atomtype
		else:
			print "not a valid molecular structure"
			return
		#extract the molecule self.moleculenum from data.atom to self.atom		 
		atomnum=1
		temp=[]
		for i in range(len(data.atoms)):
			if int(data.atoms[i][1])==self.moleculenum:
				temp.append(data.atoms[i][0]) #store extracted atomnumbers into temp
				self.atoms.append(data.atoms[i]) #store extracted atom into self.atom
				self.atoms[atomnum-1][0]=str(atomnum) #change extracted atom to the correct atomnumber
				atomnum+=1 #update atomnumber
		#extract the rest of the Lammps Molecule structures from data using changeatomnum
		if method=='atom':
			return
		else:
			self.changeatomnum(temp,data)
		
	def changeatomnum(self,temp,data=' '):
		"""changes the atomnumbers in Lmpsmolecule structures angles, bonds, dihedrals, impropers,
		and velocities in a three step process. 
		If data is defined extract the rest of the keywords used for Lmpsmolecule.
		Step 1A:If data is defined copy from data one of the above structures.
				If data is not defined skip this step.
		Step 1B:If data is defined remove the rows of copy which do not contain any values from temp
		Step 2: If data is defined extract from copy to an above Lmpsmolecule structure and remove the extracted rows
				If data is not defined skip this step.
		Step 3: Use temp and self.atoms to convert the atomnumbers in Lmpsmolecule structures
				to the correct values. 
				temp and self.atoms line up, so temp[i] and self.atoms[i][0]
				correspond to the current Lmpsmolecule structures atom numbers 
				and correct values
				will use booleanarray class to ensure that lmpsmolecule's structure values are altered only once."""
		#Step 1 and Step 2
		if data!=' ':
			#extract molecular keywords from data.keywords
			for keywords in data.keywords:
				if keywords=='Angles' or keywords=='Bonds' or keywords=='Dihedrals' or\
				keywords=='Impropers' or keywords=='Velocities':
					self.keywords.append(keywords)
			for keywords in self.keywords:
				if keywords!='Atoms': print 'extracting the data from', keywords
				if keywords=='Angles':
					copy=self.copy(data.angles) #Step 1A
					dnr=[] #dnr means do not remove a 0 means remove and a 1 means keep
					for j in range(len(copy)): #Step 1B
						dnr.append(0) #adds 0 to all list elements of dnr					
					for item in temp: #finds copied structure that has temp values 					
						for j in range(len(copy)):
							for i in range(2,len(copy[j])):
								if copy[j][i]==item: 
									dnr[j]=1 #changes jth list element of dnr to 1								
									break
					remove=[]
					for j in range(len(dnr)): # finds dnr values that are still 0
						if dnr[j]==0:
							remove.append(j) #and appends their index value to the list remove  
					copy=self.deleterows(copy,remove) #removes all unneeded rows from copy				
					structnum=1
					print 'the length of data is', len(copy)
					for item in temp: #Step 2
						found=[]
						for j in range(len(copy)): 
							for i in range(2,len(copy[j])):
								if copy[j][i]==item:
									found.append(j)
									self.angles.append(copy[j])
									l=len(self.angles)-1
									self.angles[l][0]=str(structnum)
									structnum+=1
									break
				#		print 'the item is', item
				#		print 'found is', found
						copy=self.deleterows(copy,found)
				elif keywords=='Bonds':
					copy=self.copy(data.bonds) #Step 1A
					dnr=[] #dnr means do not remove a 0 means remove and a 1 means keep
					for j in range(len(copy)): #Step 1B
						dnr.append(0) #adds 0 to all list elements of dnr					
					for item in temp: #finds copied structure that has temp values 					
						for j in range(len(copy)):
							for i in range(2,len(copy[j])):
								if copy[j][i]==item: 
									dnr[j]=1 #changes jth list element of dnr to 1								
									break
					remove=[]
					for j in range(len(dnr)): # finds dnr values that are still 0
						if dnr[j]==0:
							remove.append(j) #and appends their index value to the list remove  
					copy=self.deleterows(copy,remove) #removes all unneeded rows from copy				
					structnum=1
					print 'the length of data is', len(copy)
					for item in temp: #Step 2
						found=[]
						for j in range(len(copy)): 
							for i in range(2,len(copy[j])):
								if copy[j][i]==item:
									found.append(j)
									self.bonds.append(copy[j])
									l=len(self.bonds)-1
									self.bonds[l][0]=str(structnum)
									structnum+=1
									break
					#	print 'the item is', item
					#	print 'found is', found
						copy=self.deleterows(copy,found)						
				elif keywords=='Dihedrals':
					copy=self.copy(data.dihedrals) #Step 1A
					dnr=[] #dnr means do not remove a 0 means remove and a 1 means keep
					for j in range(len(copy)): #Step 1B
						dnr.append(0) #adds 0 to all list elements of dnr					
					for item in temp: #finds copied structure that has temp values 					
						for j in range(len(copy)):
							for i in range(2,len(copy[j])):
								if copy[j][i]==item: 
									dnr[j]=1 #changes jth list element of dnr to 1								
									break
					remove=[]
					for j in range(len(dnr)): # finds dnr values that are still 0
						if dnr[j]==0:
							remove.append(j) #and appends their index value to the list remove  
					copy=self.deleterows(copy,remove) #removes all unneeded rows from copy				
					structnum=1
					print 'the length of data is', len(copy)
					structnum=1
					for item in temp: #Step 2
						found=[]
						for j in range(len(copy)): 
							for i in range(2,len(copy[j])):
								if copy[j][i]==item:
									found.append(j)
									self.dihedrals.append(copy[j])
									l=len(self.dihedrals)-1
									self.dihedrals[l][0]=str(structnum)
									structnum+=1
									break
					#	print 'the item is', item
					#	print 'found is', found
						copy=self.deleterows(copy,found)
				elif keywords=='Impropers':
					copy=self.copy(data.impropers) #Step 1B
					dnr=[] #dnr means do not remove a 0 means remove and a 1 means keep
					for j in range(len(copy)): #Step 1B
						dnr.append(0) #adds 0 to all list elements of dnr					
					for item in temp: #finds copied structure that has temp values 					
						for j in range(len(copy)):
							for i in range(2,len(copy[j])):
								if copy[j][i]==item: 
									dnr[j]=1 #changes jth list element of dnr to 1								
									break
					remove=[]
					for j in range(len(dnr)): # finds dnr values that are still 0
						if dnr[j]==0:
							remove.append(j) #and appends their index value to the list remove  
					copy=self.deleterows(copy,remove) #removes all unneeded rows from copy				
					structnum=1
					print 'the length of data is', len(copy)
					for item in temp: #Step 2
						found=[]
						for j in range(len(copy)): 
							for i in range(2,len(copy[j])):
								if copy[j][i]==item:
									found.append(j)
									self.impropers.append(copy[j])
									l=len(self.impropers)-1
									self.impropers[l][0]=str(structnum)
									structnum+=1
									break
					#	print 'the item is', item
					#	print 'found is', found									
						copy=self.deleterows(copy,found)
				elif keywords=='Velocities':
					copy=self.copy(data.velocities) #Step 1
					dnr=[] #dnr means do not remove a 0 means remove and a 1 means keep
					for j in range(len(copy)): #Step 1B
						dnr.append(0) #adds 0 to all list elements of dnr					
					for item in temp: #finds copied structure that has temp values 					
						for j in range(len(copy)):
							for i in range(1):
								if copy[j][i]==item: 
									dnr[j]=1 #changes jth list element of dnr to 1								
									break
					remove=[]
					for j in range(len(dnr)): # finds dnr values that are still 0
						if dnr[j]==0:
							remove.append(j) #and appends their index value to the list remove  
					copy=self.deleterows(copy,remove) #removes all unneeded rows from copy				
					structnum=1
					print 'the length of data is', len(copy)
					for item in temp: #Step 2
						found=[]
						for j in range(len(copy)): 
							for i in range(1):
								if copy[j][i]==item:
									found.append(j)
									self.velocities.append(copy[j])
									l=len(self.velocities)-1
									self.velocities[l][0]=str(structnum)
									structnum+=1
									break
					#	print 'the item is', item
					#	print 'found is', found								
						copy=self.deleterows(copy,found)
									
		#Step 3
		for keywords in self.keywords:
			if keywords!='Atoms': print 'altering data structure values for', keywords
			if keywords=='Angles':				
				copy=self.copy(self.angles)
				array=booleanarray(len(copy),len(copy[0]),True) #creating a booleanarray with true values
				for i in range(len(temp)):
					if temp[i]==self.atoms[i][0]: continue
					for j in range(len(copy)): 
						for k in range(2,len(copy[j])):
							if copy[j][k]==temp[i]:
								if array.getelement(j,k): #if the boolean array is true
									self.angles[j][k]=self.atoms[i][0]
									array.setelement(j,k,False)
									break 								
			elif keywords=='Bonds':
				copy=self.copy(self.bonds)
				array=booleanarray(len(copy),len(copy[0]),True) #creating a booleanarray with true values
				for i in range(len(temp)):
					if temp[i]==self.atoms[i][0]: continue
					for j in range(len(copy)): 
						for k in range(2,len(copy[j])):
							if copy[j][k]==temp[i]:
								if array.getelement(j,k): #if the boolean array is true
									self.bonds[j][k]=self.atoms[i][0]
									array.setelement(j,k,False)
									break						
			elif keywords=='Dihedrals':
				copy=self.copy(self.dihedrals)
				array=booleanarray(len(copy),len(copy[0]),True) #creating a booleanarray with true values 
				for i in range(len(temp)):
					if temp[i]==self.atoms[i][0]: continue
					for j in range(len(copy)): 
						for k in range(2,len(copy[j])):
							if copy[j][k]==temp[i]:
								if array.getelement(j,k): #if the boolean array is true
									self.dihedrals[j][k]=self.atoms[i][0]
									array.setelement(j,k,False)
									break
			elif keywords=='Impropers':
				copy=self.copy(self.impropers)
				array=booleanarray(len(copy),len(copy[0]),True) #creating a booleanarray with true values
				for i in range(len(temp)):
					if temp[i]==self.atoms[i][0]: continue
					for j in range(len(copy)): 
						for k in range(2,len(copy[j])):
							if copy[j][k]==temp[i]:
								if array.getelement(j,k): #if the boolean array is true
									self.impropers[j][k]=self.atoms[i][0]
									array.setelement(j,k,False)
									break
			elif keywords=='Velocities':
				copy=self.copy(self.velocities)
				array=booleanarray(len(copy),len(copy[0]),True) #creating a booleanarray with true values
				for i in range(len(temp)):
					if temp[i]==self.atoms[i][0]: continue
					for j in range(len(copy)): 
						for k in range(1):
							if copy[j][k]==temp[i]:
								if array.getelement(j,k): #if the boolean array is true
									self.velocities[j][k]=self.atoms[i][0]
									array.setelement(j,k,False)
									break
	
	def copy(self,structure):
		"""copies structure to same and returns same."""
		same=[]
		for row in structure:
			same.append(row)
		return same
		
	def deleterows(self,structure,rows): # run through this {structure is list of strings and rows is list
	#of numbers}
		"""delete rows in a structure and shifts the structure up
		rows must be in increasing order for this algorithm to work correctly"""
		new=[]
		#multiple copying of b to the rows being replaced.
		if rows==[]:
			for line in structure:
				new.append(line) #if no rows need replacing copy structure
			return new
		for i in range(rows[0]):
			new.append(structure[i]) #copy structure to new until the first replaced row
		count=0
		for i in range(rows[0],len(structure)-len(rows)):# to replace rows and shift undeleted rows over
			for j in range(i+1+count,len(structure)):
				for val in rows:
					if val==j:
						count+=1
						break
				if val==j:continue
				else:
					new.append(structure[j])
					break
		return new
										
	def modifyatom(self,atomnumber,column,newvalue):
		"""modifies self.atom[atomnumber-1][column] to newvalue.
		*note: newvalue needs to be a string
		if column is 0 than changeatomnum method might need runnning for all atoms
		that have had their atomnumber(column=0) modified.
		changeatomnum method is not ran in this method when column=0"""
		self.atoms[atomnumber-1][column]=newvalue		


	def deleteatoms(self,atomnumbers,atomid): 
		"""Algorithm to find all atoms bonded to atomnumbers in the direction of atomid.
		Atomid can be a list or a single integer. The single integer corresponds to atom's atomtype value.
		The list corresponds to the atom's atomid values. The only difference between both cases is in the top portion of code. 
		When all bonded atoms have been found; ie: (the modified return values from findbonds yields an empty list);
	 	delete those bonded atoms and all molecule structures which contain those atoms.
		Convert the list of bonded atoms into rows and delete those rows from molecule.atoms"""
		print 'finding atoms to delete'
		bondedatoms=[]	
		# 1st iteration
		nextatoms=self.findbonds(atomnumbers)		
		try: # tests whether atomid is a list or a single integer
			atomid[0]
		except TypeError:
		#need to remove atoms from nextatoms which dont have the proper atom id
			testflag=True
			i=0
			while testflag:
				if i>len(nextatoms)-1: break #For the case where i becomes greater than the len of nextatoms break loop
				if int(self.atoms[nextatoms[i]-1][2])!=atomid:#uses the atom id for the row in atoms and than checks atom id
					del nextatoms[i] #delets the atom at i 
					i-=1 #and than decreases i by 1 so next atom in the list will line up when i is increased
				i+=1 #increase i by 1
		else:
		#need to remove atoms from nextatoms which dont have the proper atom id
			testflag=True
			i=0
			while testflag:
				if i>len(nextatoms)-1: break #For the case where i becomes greater than the len of nextatoms break loop
				keep=False
				for id in atomid:				
					if int(self.atoms[nextatoms[i]-1][0])==id: #checking if atomid is in next atom
						keep=True #keep this atom
						break
				if not keep:
					del nextatoms[i] #delets the atom at i 
					i-=1 #and than decreases i by 1 so next atom in the list will line up when i is increased
				i+=1 #increase i by 1
		
		#append next atoms into bondedatoms
		#copy next atoms into prevatoms
		prevatoms=[]
		for atom in nextatoms:
			bondedatoms.append(atom)
			prevatoms.append(atom)
			
							
		#2nd iteration
		if prevatoms==[]:
			print 'no bonds were found in first iteration that had atomid criteria'
			return
		
		nextatoms=self.findbonds(prevatoms)		
		#need to remove atoms from nextatoms which are in atomnumbers
		for atom in atomnumbers:
			for i in range(len(nextatoms)):
				if nextatoms[i]==atom: #checking if atom is in next atom
					del nextatoms[i] #delete the atom at i
					break
			if nextatoms==[]: break #all bonds from find bonds have already been added to bondedatoms
		
		#append next atoms into bondedatoms
		#copy next atoms into prevatoms
		prevatoms=[]
		for atom in nextatoms:
			bondedatoms.append(atom)
			prevatoms.append(atom)

		#iterative process for finding the rest of the atoms bonded to atomnumbers in the direction of atomid.
		while prevatoms!=[]:
			nextatoms=self.findbonds(prevatoms)
			#need to remove atoms from nextatoms which are in the prevatoms
			for atom in bondedatoms:
				for i in range(len(nextatoms)):
					if nextatoms[i]==atom: #checking if atom is in next atom
						del nextatoms[i] #delete the atom at i
						break
				if nextatoms==[]: break #all bonds from find bonds have already been added to bondedatoms
			
			#append next atoms into bondedatoms
			#copy next atoms into prevatoms
			prevatoms=[]
			for atom in nextatoms:
				bondedatoms.append(atom)
				prevatoms.append(atom)
	
		print 'the atoms to delete are', bondedatoms
		print 'deleting atoms from structures'
		#delete bonded atoms from the molecule's structures except the atom structure
		for i in range(1,len(self.keywords)): #goes through all keywords except the atom keyword
			print self.keywords[i]
			if self.keywords[i]=='Angles':
				rows=self.findatomnumbers(bondedatoms,self.angles,False)
				rows=self.listorder(rows) #to order the rows in increasing order
				self.angles=self.deleterows(self.angles,rows) #requires that rows be in increasing order
			elif self.keywords[i]=='Bonds':
				rows=self.findatomnumbers(bondedatoms,self.bonds,False)
				rows=self.listorder(rows) #to order the rows in increasing order
				self.bonds=self.deleterows(self.bonds,rows)#requires that rows be in increasing order
			elif self.keywords[i]=='Dihedrals':
				rows=self.findatomnumbers(bondedatoms,self.dihedrals,False)
				rows=self.listorder(rows) #to order the rows in increasing order
				self.dihedrals=self.deleterows(self.dihedrals,rows) #requires that rows be in increasing order
			elif self.keywords[i]=='Impropers':	
				rows=self.findatomnumbers(bondedatoms,self.impropers,False)
				rows=self.listorder(rows) #to order the rows in increasing order
				self.impropers=self.deleterows(self.impropers,rows) #requires that rows be in increasing order
			elif self.keywords[i]=='Velocities':
				rows=self.findatomnumbers(bondedatoms,self.velocities,True)
				rows=self.listorder(rows) #to order the rows in increasing order
				self.velocities=self.deleterows(self.velocities,rows) #requires that rows be in increasing order
 		
		print 'Atoms'
		#convert bondedatoms from atom numbers to row numbers
		for i in range(len(bondedatoms)):
			bondedatoms[i]-=1
		bondedatoms=self.listorder(bondedatoms) #to order the row numbers in increasing order
		#delete bonded atoms (row numbers) from the atom structure
		self.atoms=self.deleterows(self.atoms,bondedatoms) #requires that row numbers be in increasing order
		
	def findatomnumbers(self,atomnumbers,structure,vflag): #need to read through this algorithm
		"""Algorithm to find atomnumbers in a molecule structure except.
		Atoms structure is not handled in here.
		Returns a list of the rows in which the atomnumbers are contained in the molecule structure"""
		rows=[]
		if vflag: #for handling velocity structure
			for atom in atomnumbers:
				for i in range(len(structure)):			
					if int(structure[i][0])==atom:
						#duplicate rows for vflag=True are not possible.
						rows.append(i)
						break
		else: #for handling all other structures except atoms
			for atom in atomnumbers:
				for i in range(len(structure)):
					for j in range(2,len(structure[i])):
						if int(structure[i][j])==atom:
							#need to make sure duplicate value of rows are not being added
							if rows==[]:
								rows.append(i) #appends the row number
							else:
								#checking for duplicate values of rows
								duplicateflag=0
								for k in range(len(rows)):
									if rows[k]==i:
										duplicateflag=1
										break
								if duplicateflag==0: #if no duplicates adds bond number
									rows.append(i) #appends the row number
							break
		print 'the finished row is', rows
		return rows

	def listorder(self,struct):
		"""Takes struct and organizes the list from least to greatest.
		If the list is already ordered this algorithm will do nothing."""
		if len(struct)==1: return struct #with the length at 1; there is only one element and therefore nothing to order
		for i in range(1,len(struct)): #when i=0, struct will not change; therefore its skipped
			copy=struct[i]			
			for j in range(i-1,-1,-1):
				struct[j+1]=struct[j]
				if copy >struct[j]:
					struct[j+1]=copy
					break
				elif j==0:
					struct[j]=copy
	#dont need a break here because this is the last j value in the for loop
		print 'the organized row is', struct
		return struct
					
				
					
	def findbonds(self,atomnumbers):
		"""Algorithm to find all atoms bonded to atomnumbers.
		Returns a list of atomnumbers"""
		#finds the bonds in which the atomnumbers are located
		bondids=[]
		for i in range(len(atomnumbers)):
			for j in range(len(self.bonds)):
				for k in range(2,len(self.bonds[j])):
					if int(self.bonds[j][k])==atomnumbers[i]:
						if bondids==[]:
							bondids.append(int(self.bonds[j][0])) #appends the bond number
						else:
							#checking for duplicates of bondids
							duplicateflag=0
							for l in range(len(bondids)):
								if bondids[l]==int(self.bonds[j][0]):
									duplicateflag=1
									break
							if duplicateflag==0: #if no duplicates adds bond number
								bondids.append(int(self.bonds[j][0]))
						break
		
		#Using bondids find the atoms bonded to atomnumbers
		bondedatoms=[]
		from math import fabs
		for id in bondids:
			for atoms in atomnumbers:
				found=False
				for i in range(2,len(self.bonds[id-1])):
					if int(self.bonds[id-1][i])==atoms:
						j=int(fabs(i-5)) #switches the index from the atomnumber location to the other location
						bondedatoms.append(int(self.bonds[id-1][j])) #appends the atomnuber at the other location
						found=True
						break
				if found==True: break
		return bondedatoms	 

	def findparticlebondingpoints(self,particle,atomid,cutoffdistance,bondnumber):
		"""Particle is a particlesurfaceobject, atomid is an int.
		Finds atoms with atomid in the molecule which are less than the cutoffdistance from atoms in the particle.
		The found atoms and particles locations are stored in a 3 dimensional list called possiblebonds
		The first index corresponds with the molecule's atom
		The second index correspons with the molecule/particle combination
		The third index corresponds with whether the value is the molecule or the particle
		Always bonds the two ends of the molecule that meet cutoff requirement.
		All other possible bonds are randomly chosen until the required bondnumbers are met.
		After every bond is chosen, the particle object's boolean list is updated,	and possiblebonds is updated.
		The update to possiblebonds involves removing the row from which the bonded molecule's atom is located
		and also removing the particle atom and it's corresponding bonded atom from other rows of possiblebonds. 
		The final bonds are all stored in self.bonding as a 2 dimensional list"""
		possiblebonds=[]
		for i in range(len(self.atoms)):
			if int(self.atoms[i][2])==atomid:
				row=[]
				for j in range(len(particle.surface)):
					#assumes molecule and particle have same atomtype if not true than atomdistance will give bad value
					if atomdistance(self.atoms[i],particle.surface[j],self.atomtype)<=cutoffdistance:
						if particle.getbool(j): row.append([i,j]) #if not bonded than can form a possible bond
				if row!=[]:possiblebonds.append(row) #need to correct this....
		
		#initiate section which assigns bonds into bonding information
		bonds=0
		self.bondinginformation=[] #initiate new class member

		#Checks to see if no bonds are possible [2 possible cases]
		if possiblebonds==[]: 
			print 'no possible bonds can be formed'
			return
		if bondnumber==0:
			print 'bondnumber is 0; so, no possible bonds can form'
			return	

		#section which assigns a bond to the first molecule atom which can be bonded.	
		self.bondinginformation.append(self.particlebondinginfo(particle,possiblebonds,0))
		bonds+=1
		if bonds==bondnumber: return
		del possiblebonds[0] #deletes possible bonds to the first molecule which can be bonded
		l=len(self.bondinginformation)-1
		possiblebonds=self.updatepossiblebonds(possiblebonds,self.bondinginformation[l][1]) #updates possiblebonds 
		#by removing any bonds which contain the newly bonded particle.
		
		if possiblebonds==[]:return
		#section which finds the last molecule atom which can be bonded
		l=len(possiblebonds)-1
		while possiblebonds[l]==[]:# to find the last molecule atom which can be bonded
			if possiblebonds==[]:return
			del possiblebonds[l] #since there are no more possible bonds in this location delete
 			l-=1 #go to next possible spot where the last molecule atom could be bonded

		#section which assigns a bond to the last molecule atom which can be bonded.
		self.bondinginformation.append(self.particlebondinginfo(particle,possiblebonds,l))			 
		bonds+=1
		if bonds==bondnumber: return
		del possiblebonds[l] #deletes possible bonds to the last molecule which can be bonded
		l=len(self.bondinginformation)-1
		possiblebonds=self.updatepossiblebonds(possiblebonds,self.bondinginformation[l][1])
		
		if bondnumber-bonds>=len(possiblebonds): #the rest of the bonds are assigned in order
			while possiblebonds!=[]:
				if possiblebonds[0]==[]: #if row of possible bonds is empty then delete the row
					del possiblebonds[0]				
				else: #else find a bond in the row of possible bonds 
					self.bondinginformation.append(self.particlebondinginfo(particle,possiblebonds,0))
				#dont need to update bonds since possiblebonds will become empty at the same time or before bonds
				#is equal to bondnumber
					del possiblebonds[0]
					l=len(self.bondinginformation)-1
					possiblebonds=self.updatepossiblebonds(possiblebonds,self.bondinginformation[l][1])
			return
		else: #the rest of the bonds are assigned randomly
			from random import randint #use to randomly choose an index to bond.
			while bonds<bondnumber:
				if possiblebonds==[]:break #exits while loop when possiblebonds has become an empty set.
				#this ensures that in the case there are not egnough viable bonds from possiblebonds to
				#reach the bondnumber. Than, the while loop will not become infinite.
				l=len(possiblebonds)-1
				i=randint(0,l)
				if possiblebonds[i]==[]: #if row of possible bonds is empty then delete the row
					del possiblebonds[i]
		 		else: #else find a bond in the row of possible bonds 
					self.bondinginformation.append(self.particlebondinginfo(particle,possiblebonds,i))
					bonds+=1
					del possiblebonds[i]
					l=len(self.bondinginformation)-1
					possiblebonds=self.updatepossiblebonds(possiblebonds,self.bondinginformation[l][1])
			return

	def particlebondinginfo(self,particle,possiblebonds,index1):
		"""Takes list of possiblebonds and assigns a bond from possiblebonds[index1].
		Then assigns false to particle.setbool at the bonding location 
		to ensure no more bonds can form with that particle.
		Returns the assigned bond."""
		from random import randint #use to randomly choose 0 and 1 where 1 is keep.
		for i in range(len(possiblebonds[index1])):
			if i==len(possiblebonds[index1])-1: #automattically bonds under this condition
				break								
			else: #bond has random chance of forming
				if randint(0,1)==1: #bond forms
					break
				else: #no bond forms
					continue		
		particle.setbool(possiblebonds[index1][i][1],False)#makes sure no more bonds can form
		return possiblebonds[index1][i]

	def updatepossiblebonds(self,possiblebonds,particlenum):
		"""finds the particle number in the remaining possible bonds than deletes that bonding information. 
		This algorithm assumes that particlenum can exist only once in each row of possiblebonds."""
		for i in range(len(possiblebonds)):
			for j in range(len(possiblebonds[i])):
				if possiblebonds[i][j][1]==particlenum:
					del possiblebonds[i][j]
					break #go to next row of possiblebonds
		return possiblebonds
	
	def bondtoparticle(self,particle,atomid,newid,newcharge):
		"""Bonds molecule to particle. Particle is a particlesurface object. 
		Moves the molecule's atoms bonded to the particle to the particle's position.
		Alters the atomtype of the molecule's atoms bonded to the particle to newid.
		Alters the charge of the molecule's atoms bonded to the particle to newcharge 	
		Because bonds have formed between the molecule's atoms and the particle's atoms,
		atoms with atomid on the molecule need to be removed otherwise the molecule's atoms will have to many bonds.
		self.deleteatoms will take care of deleting the atoms atomid and 
		atoms down the polymer chain in the direction away from the molecule/particle bond.
		Now remove the bonded particle atoms from the particle surface becaused the bonded molecule atoms
		have replaced those particle atoms and their bonds with the rest of the particle.
		Newid must be a string. Atomid can now be a list or an integer.
		The list is a list of atom's id values and the single integer is atom's atomtype"""
		#moves the molecule's atoms bonded to the particle to the particle's position
		#need to do this step first before the atom id's and row indexes get out of sync 
		#which will occur after atoms get deleted.
		print 'modifying the bonded molecules atom information'
		for i in range(len(self.bondinginformation)):
			#patom=particle.getsurfatom(self.bondinginformation[i][1])
			#if self.atomtype=='molecular':#3-5->x,y,z[molecular]
			#	for j in range(3,6):
			#		self.modifyatom(self.bondinginformation[i][0]+1,j,patom[j])#uses the atomid here
			#elif self.atomtype=='full':#4-6->x,y,z[full]
			#	for j in range(4,7):		
			#		self.modifyatom(self.bondinginformation[i][0]+1,j,patom[j])#uses the atomid here
			#alters the atomtype to newid of the molecule's atoms bonded to the particle.
			self.modifyatom(self.bondinginformation[i][0]+1,2,newid) #uses the atomid here
			if self.atomtype=='full':
				# Alters the charge of the molecule's atoms bonded to the particle to newcharge
				self.modifyatom(self.bondinginformation[i][0]+1,3,newcharge)			

		#This is a separate loop so atom id's and row indexes for previous steps wont get out of sync		
		#create atomnumbers to begin deletion process.
		atomnumbers=[] 
		for i in range(len(self.bondinginformation)):
			atomnumbers.append(self.bondinginformation[i][0]+1)#using actual atomnumbers rather than the row index
						
		#Call algorithm to find all atoms bonded to atomnumbers in the direction of atomid.
		#Than delete those atoms and all molecule structures which contain those atoms.
		self.deleteatoms(atomnumbers,atomid)
		
		print 'beginning deletion process of the surface atoms for which the molecule atoms have replaced'
		#Goes through the bondinginformation and superficially removes the surfaceatom
		for i in range(len(self.bondinginformation)):
			particle.removesurfatom(self.bondinginformation[i][1]) #uses the row number
			#Allows the particle extract method used outside this class to remove this atom.

	def createxyz(self,file,data,routine='mass', values=None):
		"""Two possible routines one to use the masses from data and the other to use the atom type and values supplied by the user.
		The mass version is assessed by setting the routine to 'mass' which is the default method.
		The other version is assessed by setting the routine to 'atomtype'.
		The other version takes values which is a list containing the value the user wants to assign those atomtypes to.
		The atomtypes of the values in the list will start at 1 even if no atoms in molecule use 1.
		This makes it easier to find the atomtype and assign the value from the list
		All atom data is assumed to have image flags in the data."""
		f=open(file,'w')
		f.write('{0}\n'.format(len(self.atoms)))
		f.write('atoms\n')
		if routine=='mass':
			for line in self.atoms:
				type=int(line[2])-1
				mass=data.masses[type][1]
				if mass=='12.0107': elementnum=6
				elif mass=='15.9994': elementnum=8
				elif mass=='26.9815':elementnum=13
				else:
					print 'no matching mass value. The method will exit'
					return
				l=len(line)-1
				f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-5],line[l-4],line[l-3]))
			f.close()
		elif routine=='atomtype':
			for line in self.atoms:
				type=int(line[2])-1
				elementnum=values[type]
				l=len(line)-1
				f.write('{0} {1} {2} {3}\n'.format(elementnum,line[l-5],line[l-4],line[l-3]))				
			f.close()
