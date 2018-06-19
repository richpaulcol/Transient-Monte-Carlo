import numpy as np
import pylab as pp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation

def normal_dist(mean,variance,Min,Max):
	stddev = np.sqrt(variance)
	#x = np.linspace(mean-3*stddev,mean+3*stddev,100)
	x = np.arange(Min,Max,1./1000)
	y = 1./(np.sqrt(2 * np.pi * variance)) * np.exp(- (x-mean)**2 / (2*variance))
	return x,y

def auto_normal_dist(mean,variance):
	stddev = np.sqrt(variance)
	x = np.linspace(mean-6*stddev,mean+6*stddev,10000)
	#x = np.arange(Min,Max,1./1000)
	y = 1./(np.sqrt(2 * np.pi * variance)) * np.exp(- (x-mean)**2 / (2*variance))
	return x,y

class SubplotAnimation(animation.TimedAnimation):
	def __init__(self):
		fig = pp.figure(figsize=(10,6))
		self.States = np.load('States.npy')
		self.Covariances = np.load('Covariances.npy')
		self.times = np.load('Times.npy')
		self.IncomingTMass = np.load('Incoming.npy',)
		self.tau_s0s = np.load('tau_s.npy',)
		self.Q = np.load('Q.npy')
		self.T = np.load('Turbidity.npy')

		ax1 = pp.subplot2grid((4, 3), (0, 0),colspan = 2,rowspan = 2)
		ax1.plot(self.times,self.Q,'k--',label = 'Measured Flow')
		ax1.set_ylabel('Applied Shear (Pa)')
		self.line1 = Line2D([],[],color='red')
		self.line1a = Line2D([],[],color='red',alpha = 0.5)
		self.line1b = Line2D([],[],color='red',alpha = 0.5)
		ax1.add_line(self.line1)
		ax1.add_line(self.line1a)
		ax1.add_line(self.line1b)




		ax2 = pp.subplot2grid((4, 3), (2, 0),colspan = 2,rowspan = 2)
		ax2.plot(self.times,self.T,'k--',label = 'Measured Turbidity')
		ax2.set_ylabel('Turbidity (turbidity units)')
		ax2.set_xlabel('Time (time units)')
		self.line2 = Line2D([],[],color='red')
		self.line2a = Line2D([],[],color='red',alpha = 0.5)
		self.line2b = Line2D([],[],color='red',alpha = 0.5)
		
		ax2.add_line(self.line2)
		ax2.add_line(self.line2a)
		ax2.add_line(self.line2b)

		self.ax3 = pp.subplot2grid((4, 3), (0, 2))
		self.ax3.set_xlabel('Upstream Turbidity (turbidity units)')
		self.ax3.set_ylim(0,10)
		
		self.line3 = Line2D([],[],color='red')
		self.line3a = Line2D([],[],color = 'black',linestyle = '--')
		self.ax3.add_line(self.line3)
		self.ax3.add_line(self.line3a)

		self.ax4 = pp.subplot2grid((4, 3), (1, 2))
		self.ax4.set_xlabel('tau_s (N/m2)')
		self.ax4.set_ylim(0,10.)
		
		self.line4 = Line2D([],[],color='red')
		self.line4a = Line2D([],[],color = 'black',linestyle = '--')
		self.ax4.add_line(self.line4)
		self.ax4.add_line(self.line4a)

		self.ax5 = pp.subplot2grid((4, 3), (2, 2))
		self.ax5.plot([5,5],[0,1000000000],'k--')
		self.ax5.set_xlabel('alpha (mass/m2)')
		self.line5 = Line2D([],[],color='red')
		self.ax5.add_line(self.line5)

		self.ax6 = pp.subplot2grid((4, 3), (3, 2))
		self.ax6.plot([0.2,0.2],[0,100000000],'k--')
		self.ax6.set_xlabel('beta (1/s)')
		self.line6 = Line2D([],[],color='red')
		self.ax6.add_line(self.line6)
		
		pp.tight_layout()
		animation.TimedAnimation.__init__(self, fig, interval=50, blit=False)


	def _draw_frame(self, framedata):
		i = framedata

		self.line1.set_data(self.times[:i], self.States[-4,:i])
		self.line1a.set_data(self.times[:i],self.States[-4,:i]+self.Covariances[-4,-4,:i])
		self.line1b.set_data(self.times[:i],self.States[-4,:i]-self.Covariances[-4,-4,:i])
		self.line2.set_data(self.times[:i], self.States[-6,:i])
		self.line2a.set_data(self.times[:i],self.States[-6,:i]+self.Covariances[-6,-6,:i])
		self.line2b.set_data(self.times[:i],self.States[-6,:i]-self.Covariances[-6,-6,:i])
		
		x,y = auto_normal_dist(self.States[-5,i],self.Covariances[-5,-5,i]**2)
		self.line3.set_data(x,y)
		self.ax3.set_ylim(0,max(y))
		self.ax3.set_xlim(0,+0.5)
		self.line3a.set_data([self.IncomingTMass[i],self.IncomingTMass[i]],[0,max(y)])
		
		x,y = auto_normal_dist(self.States[-3,i],self.Covariances[-3,-3,i]**2)
		self.line4.set_data(x,y)
		self.ax4.set_ylim(0,max(y))
		self.ax4.set_xlim(0,5.0)
		self.line4a.set_data([self.tau_s0s[i],self.tau_s0s[i]],[0,max(y)])

		x,y = auto_normal_dist(self.States[-2,i],self.Covariances[-2,-2,i]**2)
		self.line5.set_data(x,y)
		self.ax5.set_ylim(0,max(y))
		self.ax5.set_xlim(0.,10.)
		


		x,y = auto_normal_dist(self.States[-1,i],self.Covariances[-1,-1,i]**2)
		self.line6.set_data(x,y)
		self.ax6.set_ylim(0,max(y))
		self.ax6.set_xlim(0.,1.0)
		
		self._drawn_artists = [self.line1,self.line1a,self.line1b, self.line2,self.line2a,self.line2b,self.line3,self.line3a,self.line4,self.line4a,self.line5,self.line6]

	def new_frame_seq(self):
        	return iter(range(self.times.size))

    	def _init_draw(self):
        	lines = [self.line1,self.line1a,self.line1b, self.line2,self.line2a,self.line2b,self.line3,self.line3a,self.line4,self.line4a,self.line5,self.line6]
       		for l in lines:
            		l.set_data([], [])
		
                              
ani = SubplotAnimation()
#ani.save('Alpha_Beta_Train_OnlinePODDS1.mp4')
pp.show()
