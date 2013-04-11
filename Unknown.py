import scipy.optimize as scop
class Unknown(object):
	def __init__(self,name,sterm=0,mterm=1,reciprocal=False):
		self.name = name
		self.reciprocal = reciprocal
		self.sterm=sterm
		self.mterm=mterm
		
	def __repr__(self):
		return "%s((%s*Unknown(%s))+%s)" % ("1./" if self.reciprocal else "",repr(self.mterm),repr(self.name),repr(self.sterm))
	
	def __str__(self):
		return "%s((%s*%s)+%s)" % ("1./" if self.reciprocal else "",str(self.mterm),str(self.name),str(self.sterm))
	
	def __mul__(self, other):
		oterm = (1./other) if self.reciprocal else other
		return Unknown(name=self.name,sterm=self.sterm*oterm,mterm=self.mterm*oterm,reciprocal=self.reciprocal)
	
	def __div__(self, other):
		return self.__mul__(1./other)
	
	def __add__(self, other):
		if self.reciprocal:
			# This is where it gets hairy
			denom = Unknown(name=self.name,sterm=self.sterm,mterm=self.mterm,reciprocal=True)
			return Unknown(name=self.name,sterm=(self.sterm*other+1)*denom,mterm=(self.mterm*other)*denom,reciprocal=False)
		else:
			return Unknown(name=self.name,sterm=self.sterm+other,mterm=self.mterm,reciprocal=False)
		self.sterm = self.sterm + other
		return self
	
	def __sub__(self, other):
		return self.__add__(-other)
		
	def __rmul__(self, other):
		return self.__mul__(other)
	
	# Division is not commutative!
	def __rdiv__(self, other):
		return Unknown(name=self.name,sterm=self.sterm,mterm=self.mterm,reciprocal = not self.reciprocal).__rmul__(other)
	
	def __radd__(self, other):
		return self.__add__(other)
	
	def __rsub__(self, other):
		return self.__neg__().__radd__(other)
	
	def __neg__(self):
		return Unknown(name=self.name,sterm=-self.sterm,mterm=-self.mterm,reciprocal=self.reciprocal)
	
	#Evaluate one level deep
	def __call__(self, **variables):
		tterm = Unknown(self.name)
		mdict = {}
		sdict = {}
		for n,v in variables.iteritems():
			if self.name == n:
				tterm = v
			if type(self.mterm) == Unknown and n in self.mterm:
				mdict[n] = v
			if type(self.sterm) == Unknown and n in self.sterm:
				sdict[n] = v
		term = (self.sterm(**sdict) if sdict else self.sterm) + (tterm * (self.mterm(**mdict) if mdict else self.mterm))
		return 1/term if self.reciprocal else term
	
	def __contains__(self, key):
		ret = False
		if self.name == key:
			ret = True
		if type(self.mterm) == Unknown:
			ret |= self.mterm.__contains__(key)
		if type(self.sterm) == Unknown:
			ret |= self.sterm.__contains__(key)
		return ret
	
	def findaroot(self, **guesses):
		args = sorted(guesses.keys())
		init_vals = array([guesses[a] for a in args])
		def redict(argorder,values):
			return {a:v for a,v in zip(args,values)}
		
		def solvable(array_values):
			return self.__call__(**redict(args,array_values))
			
		def sstatus(x,f):
			print "Trying %s to get %f" % (str(x),self.__call__(**redict(args,x)))
		
			
		result = scop.root(solvable,init_vals,method='excitingmixing',tol=1e-10)
		if result.success:
			ret = redict(args,result.x)
		else:
			ret = {}
		return ret
	
	def collect_additive_terms(self,name):
		if self.reciprocal:
			raise ValueError("Can't collect additive terms of inverted Unknown: %s" % repr(self))
		else:
			usterm = type(self.sterm) == Unknown
			sub_terms = self.sterm.collect_additive_terms(name) if usterm else (0,0)
			return ((self.mterm if name==self.name else 0)+sub_terms[0],(0 if usterm else self.sterm)+sub_terms[1]+(Unknown(name=self.name,mterm=self.mterm) if name!=self.name else 0))
		
	# not sure if this is safe with other variables...
	def collapse_on(self, name):
		m,s = self.collect_additive_terms(name)
		return Unknown(name=name,mterm=m,sterm=s)
	
	

