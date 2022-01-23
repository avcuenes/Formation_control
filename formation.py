import numpy as np



class formation:
	
	def __init__(self, x_i: list,y_i: list, x1_h: int, y1_h: int, numUAV: int, angle: float, mes: float, form: int):
		"""
		This class create formation respect to given final and final position of leader agent , angle of formation , distance between vehicle
		, number of agent and formation type

		:x_i: İnitial x position of agent
		:y_i: İnitial y position of agent
		:x1_h: Final or target x position of leader agent
		:y1_h: Final or target y position of leader agent
		:numUAV: Number of agent
		:angle: Angle of formation
		:mes: Distance between agents
		:form: Formation type ===== #circle # triangle # square

		:return: Final position of all agent
		"""
		self.x_i = x_i #
		self.y_i = y_i
		self.x1_h = x1_h
		self.y1_h = y1_h
		self.numUAV = numUAV 
		self.angle = angle 
		self.mes = mes
		self.form = form

	def formation_cal(self):
		"""
		This function is main function and calculate formation position of the other agent respect to agent
		This calculation is based on graph theory

		:return: Final positon of agent respect to formation type and leader final position
		"""
		A, L = self.lagran()
		if self.form ==1:
			R_x, R_y = self.triangle()
		elif self.form ==2:
			R_x, R_y = self.square()
		elif self.form ==3:
			R_x, R_y = self.circle()
		else:
			R_x, R_y = self.triangle()

		S_x = np.matmul(A,R_x)
		S_y = np.matmul(A,R_y)
		x = self.x_i
		y = self.y_i
		e = 0.25
		for j in range(0,30):
			theta = self.angle*np.pi/180
			for i in range(0,self.numUAV):
				L_x = (np.matmul(L, x))
				L_y = (np.matmul(L, y))
				diag_x = np.diagonal(S_x)
				diag_y = np.diagonal(S_y)
				pos = np.array(((diag_x[i]), (diag_y[i])))
				c, s = np.cos(theta), np.sin(theta)
				R = np.array(((c,-s), (s,c)))
				pos = np.matmul(R, pos)
				x[i] = x[i]-1*e*(L_x[i]-pos[0])
				y[i] = y[i]-1*e*(L_y[i]-pos[1])
			x[0]=self.x1_h
			y[0]=self.y1_h
		return x,y

	def lagran(self):
		"""
		This function create adjancy matrix and lagragian matrix on graph theory
		This function calculate two matrix this matrix keep in relative state of the other agent respect to leader

		:return: Adjancy matrix== A /// Lagragian Matrix === L
		"""
		
		A = np.zeros((self.numUAV, self.numUAV))
		L = np.zeros((self.numUAV, self.numUAV))
		
		for i in range(0,self.numUAV):
			a = 0
			if i <= 3:
				for j in range(0,self.numUAV):
					a = a +1
					if i==0:
						if j == 0:
							A[i][j] = 1
					else:
						if i == j :
							A[i][j] = 0
						else:
							if a <=3 and i<=3:
								A[i][j] = 1
								
							else:
								A[i][j] = 0
			else:
				for j in range(i-3,i+1):
					a = a +1
					if i==0:
						if j == 0:
							A[i][j] = 1
					else:
						if i == j :
							A[i][j] = 0
						else:
							if a <4:
								A[i][j] = 1
							else:
								A[i][j] = 0
		
		for i in range(0,self.numUAV):
			a = 0
			if i <= 3:
				for j in range(0,self.numUAV):
					a = a +1
					if i==0:
						if j == 0:
							L[i][j] = 0
					else:
						if i == j :
							L[i][j] = 3
						else:
							if a <4:
								L[i][j] = -1
							else:
								L[i][j] = 0
			else:
				for j in range(i-3,i+1):
					a = a +1
					if i==0:
						if j == 0:
							L[i][j] = 0
					else:
						if i == j :
							L[i][j] = 3
						else:
							if a <4:
								L[i][j] = -1
							else:
								L[i][j] = 0
		
		A[1][3] = 1
		L[1][3] = -1
		A[2][3] = 1
		L[2][3] = -1
		return A,L

	def triangle(self):
		"""
		This function calculate agents relative position to make formation type triangle

		:return: Relative x position == R_x /// Relative y positon == R_y
		"""
		x_f = [0,-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7]
		y_f = [0,1,1,2,2,3,3,4,4,5,5,6,6,7,7]

		R_x =np.zeros((self.numUAV,self.numUAV))
		R_y =np.zeros((self.numUAV,self.numUAV))

		for i in range(0,self.numUAV):
			for j in range(0,self.numUAV):
				R_x[i][j] = x_f[j]-x_f[i]
				R_y[i][j] = y_f[i]-y_f[j]
		R_x =self.mes*np.array(R_x)
		R_y =self.mes*np.array(R_y)
		
		return R_x,R_y

	def square(self):
		"""
		This function calculate agents relative position to make formation type square

		:return: Relative x position == R_x /// Relative y positon == R_y
		"""
		x_f = [0,-1,1,-2,2,-2,2,-2,2,-2,2,-2,2,-1,1]
		y_f = [0,0,0,0,0,-1,-1,-2,-2,-3,-3,-4,-4,-4,-4]

		R_x =np.zeros((self.numUAV,self.numUAV))
		R_y =np.zeros((self.numUAV,self.numUAV))

		for i in range(0,self.numUAV):
			for j in range(0,self.numUAV):
				R_x[i][j] = x_f[j]-x_f[i]
				R_y[i][j] = y_f[j]-y_f[i]
		R_x =self.mes*np.array(R_x)
		R_y =self.mes*np.array(R_y)

		return R_x,R_y
		
	def circle(self):
		"""
		This function calculate agents relative position to make formation type circle

		:return: Relative x position == R_x /// Relative y positon == R_y
		"""
		xi = 0
		yi= 1
		th = (360/self.numUAV)*np.pi/180
		x_f = []
		y_f = []
		ah = th
		for i in range(0,self.numUAV):
			pos = np.array((xi,yi))
			c,s = np.cos(ah),np.sin(ah)
			R = np.array(((c,-s),(s,c)))
			if (i%2)==0:
				ah = (i+1)*th
			else:
				ah = -(i+1)*th
			pos = np.matmul(R,pos)
			xi,yi = pos[0],pos[1]
			x_f.append(pos[0])
			y_f.append(pos[1])
		R_x =np.zeros((len(x_f),len(y_f)))
		R_y =np.zeros((len(x_f),len(y_f)))
		
		for i in range(0,len(x_f)):
			for j in range(0,len(y_f)):
				R_x[i][j] = x_f[j]-x_f[i]
				R_y[i][j] = y_f[j]-y_f[i]
		R_x =self.mes*np.array(R_x)
		R_y =self.mes*np.array(R_y)

		return R_x,R_y
		


