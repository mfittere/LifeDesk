import os as os
import sys as sys
from subprocess import Popen,PIPE
from inspect import getfile

try:
  import numpy as np
  import matplotlib
  if 'DISPLAY' not in os.environ:
      matplotlib.use('Agg')
  import matplotlib.pyplot as pl
  import scipy.signal
except ImportError:
  print "No module found: numpy matplotlib and scipy modules should be present to run sixdb"
  raise ImportError

class LifeDeskDB(object):
  def __init__(self,ltr_dir='.',ltr_file='lhc.ltr',plt_dir='.',input_param={},data=np.array([])):
    """creates the LifeDeskDB class object
    Parameters:
    -----------
    ltr_dir : directory with tracking data, default = '.'
    ltr_file: file with tracking data, default = 'lhc.ltr'
    plt_dir : directory to save plots, default = '.'
    input_param: dictionary containing input parameters
        'Monitoring IP': monitoring IP
        'steplen'      : step length = number of turns tracked per step
        'nsteps'       : number of steps, total number of turns tracked
                         is nturn = nstep * steplen
        'gamma'        : relativistic gamma
        'circlhc'      : LHC circumference (constant = 26658.8832)
    data: structured array containing the output data
      step       : step number
      nturn      : turn number
      time       : time [s]
      emit1      : hor. normalized emittance [mum]
      emit2      : vert. normalized emittance [mum]
      sigm       : bunch length [cm]
      intensity  : normalized beam intensity I_tot/I_0
      luminosity :
      lossrate   : normalized loss rate I_lost/I_0

    """
    if(ltr_dir=='.'): ltr_dir= os.getcwd()
    self.lifedeskenv={'ltr_dir':ltr_dir,'ltr_file':ltr_file,'plt_dir':plt_dir}
    self.input_param=input_param
    self.data = data
# dictionary which stores the lbls and units
    self._unit = {'step': ('step number',''),'nturn': ('number of turns',''),'time': ('time','[s]'),'emit1':('hor. emittance','[$\mu$m]'),'emit2':('vert. emittance','[$\mu$m]'),'sigm':('bunch length','[cm]'),'intensity': ('normalized beam intensity',''),'lossrate':('normalized loss rate','')}
  @classmethod
  def getdata(cls,ltr_dir='.',ltr_file='lhc.ltr',plt_dir='.'):
    '''create LifeDeskDB class object from dat
    in directory ltr_dir.
    Parameters:
    -----------
    ltr_dir : directory with tracking data, default = '.'
    ltr_file: file with tracking data, default = 'lhc.ltr
    '''
# check if ltr_dir/ltr_file exists
    if not os.path.isfile(os.path.join(ltr_dir,ltr_file)):
      print 'ERROR in getdata: file %s does not exist!'%(os.path.join(ltr_dir,ltr_file))
      return
    print '... getting data from %s'%ltr_file
# get path to lhcpost, output.sh files in LifeDesk directory
    script_path='%s/scripts'%(os.path.dirname(getfile(LifeDeskDB)))
# run lhcpost to get input parameters
# and output files: emit.txt,intensity.txt,lossrate.txt,luminosity.txt,sigm.txt
    if(ltr_dir=='.'): ltr_dir= os.getcwd()
    if(plt_dir=='.'): plt_dir= os.getcwd()
    if not os.path.isdir(plt_dir):
      os.makedirs(plt_dir)
# delete old output files
    for ff in 'emit.txt intensity.txt lossrate.txt luminosity.txt'.split():
      if os.path.isfile(os.path.join(ltr_dir,ff)):
        os.remove(os.path.join(ltr_dir,ff))
        print 'deleted file %s'%(ff)
# create new output files
    print '... calling script %s/lhcpost'%(script_path)
    p = Popen(['%s/lhcpost'%(script_path),ltr_dir,ltr_file,script_path], stdout=PIPE)
    stdout,stderr=p.communicate()
    if stderr != None:
      print "ERROR while executing command lhcpost %s %s:"%(ltr_dir,ltr_file)
      print stderr
      return
    check=True
    for ff in 'emit.txt intensity.txt lossrate.txt luminosity.txt'.split():
      if not os.path.isfile(os.path.join(ltr_dir,ff)):
        print "ERROR when calling output.sh: file %s has not been generated!"%ff
        check=False
    if check: print "... created emit.txt, intensity.txt, lossrate.txt, luminosity.txt"
# -- get in put parameters
    lout=stdout.split('\n')
    input_param={}
    for s in lout:
      if s.find('Monitoring IP')>-1:
        input_param['Monitoring IP'] = s.split()[-1]
      if s.find('steplen')>-1:
        input_param['steplen'] = float((s.split(',')[0]).split()[-1])
        input_param['nsteps']  = float((s.split(',')[1]).split()[-1])
      if s.find('gamma')>-1:
        input_param['gamma'] = float(s.split()[-1])
    input_param['circlhc'] = 26658.8832
# -- read in emit.txt,
    clight = 299792458 # speed of light
    gamma  = input_param['gamma']
    beta   = np.sqrt(1-1/gamma**2)
    emit1,emit2,step = (np.loadtxt('%s/emit.txt'%(ltr_dir))).T # rms emittance [cm]
    emit1 = emit1*(beta*gamma)*1.e4 # normalized emittance [mum]
    emit2 = emit2*(beta*gamma)*1.e4 # normalized emittance [mum]
    nturn            = step * input_param['steplen']
    time             = nturn*input_param['circlhc']/(beta*clight) # time [s]
    sigm             = np.loadtxt('%s/sigm.txt'%(ltr_dir))
    intensity        = np.loadtxt('%s/intensity.txt'%(ltr_dir))
    lumi             = np.loadtxt('%s/luminosity.txt'%(ltr_dir))
    loss             = np.loadtxt('%s/lossrate.txt'%(ltr_dir))
    data_flat = (np.column_stack((step,nturn,time,emit1,emit2,sigm,intensity,lumi,loss))).ravel()
    ftype=[('step',float),('nturn',float),('time',float),('emit1',float),('emit2',float),('sigm',float),('intensity',float),('luminosity',float),('lossrate',float)]
    data = data_flat.view(ftype)
    db=cls(ltr_dir,ltr_file,plt_dir,input_param,data)
    return db
  def set_env(self,ltr_dir=None,ltr_file=None,plt_dir=None):
    for n,v in zip(['ltr_dir','ltr_file','plt_dir'],[ltr_dir,ltr_file,plt_dir]):
      if v !=None: self.lifedeskenv[n]=v
  def print_env(self):
    print 'tracking directory: ltr_dir  = %s'%(self.lifedeskenv['ltr_dir'])
    print 'output file name  : ltr_file = %s'%(self.lifedeskenv['ltr_file'])
    print 'plot directory: plt_dir  = %s'%(self.lifedeskenv['plt_dir'])
  def plot_2d(self,xaxis='time',yaxis='emit1',color='b',lbl=None,alpha=1.0):
    """plot *xaxis* vs *yaxis*
    Parameters:
    -----------
    param   : hor. (mode=1) or vert. (mode=2) emittance
    tunit   : 'nturn': number of turns
              'time' : time [s]
    norm    : eps = beta*gamma*epsn
              false: rms emittance eps [mum]
              true : normalized emittance epsn [mum]
    color   : plot color
    """
    data_names = self.data.dtype.names
    # check if fields exist
    for n in xaxis,yaxis:
      if n not in data_names:
        print '%s not found in data'%n
        return 0
    x,y = self.data[xaxis],self.data[yaxis]
    pl.plot(x,y,linestyle='-',marker='o',color=color,label=lbl,alpha=alpha)
    pl.xlabel(r'%s %s'%self._unit[xaxis])
    pl.ylabel(r'%s %s'%self._unit[yaxis])
    pl.legend(loc='best',fontsize=12)
  def plot_all(self,color='b',lbl=None,title=None,export=None,alpha=1.0):
    """plots emittance, bunch length, intensity
    luminosity and loss rate vs time [s].
    
    Parameters:
    -----------
    color: plot color
    lbl: plot label
    export: export format, e.g. 'png'
    """
    for p in ['emit1','emit2','sigm','intensity','luminosity','lossrate']:
      pl.figure(p)
      try:
        self.plot_2d(xaxis='time',yaxis=p,color=color,lbl=lbl,alpha=1.0)
        # place the legend to the outside of the plot
        box = pl.gca().get_position()
        # if more than 4 entries, take two columns, otherwise one
        nlabels = len(pl.gca().get_legend_handles_labels()[1])
        if (nlabels <=4) : ncol = 1
        elif (nlabels > 5 and nlabels <=8): ncol = 2
        else: ncol = 3
        # legend on top
        pl.legend(bbox_to_anchor=(0., 1.07, 1.0, .102), loc=3,ncol=ncol, mode="expand", borderaxespad=0.,fontsize=12,title=title)
        pl.subplots_adjust(left=0.15, right=0.95, top=0.7, bottom=0.1)
#        # legend on right side
#        legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=12)
#        pl.subplots_adjust(left=0.15, right=0.8, top=0.1, bottom=0.1)
        if export != None:
          print '%s.%s'%(p,export)
          pl.savefig('%s/%s.%s'%(self.lifedeskenv['plt_dir'],p,export))
      except KeyError:
        print 'ERROR in plot_all: could not plot %s'%p
        print '   self.data[%s][0]=%s'%(p,self.data[p][0])
