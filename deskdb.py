import os as os
import sys as sys
from subprocess import Popen,PIPE
from inspect import getfile
from matplotlib import gridspec
import shutil as shutil

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

mycolors=['b','r','g','m','orange','pink','cyan','indigo','lime']
def colorrotate():
  c=mycolors.pop(0);mycolors.append(c)
  return c
def gaussian(x,mean,sigma):
#  return np.exp(-((x-mean)**2)/(2*sigma**2))/(np.sqrt(2*sigma**2*np.pi))
  return np.exp(-((x-mean)**2)/(2*sigma**2))

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
    hist: histograms (only filled if self.get_hist()
          is called)
      'header' : mean,sigma and norm of histograms
      'data'   : histogram data
    """
    if(ltr_dir=='.'): ltr_dir= os.getcwd()
    self.lifedeskenv={'ltr_dir':ltr_dir,'ltr_file':ltr_file,'plt_dir':plt_dir}
    self.input_param=input_param
    self.data = data
    self.hist = {}
# dictionary which stores the lbls and units
    self._unit = {'step': ('step number',''),'nturn': ('number of turns',''),'time': ('time','[s]'),'emit1':('hor. emittance','[$\mu$m]'),'emit2':('vert. emittance','[$\mu$m]'),'sigm':('bunch length','[cm]'),'intensity': ('normalized beam intensity',''),'lossrate':('normalized loss rate','')}
  @classmethod
  def getdata(cls,ltr_dir='.',ltr_file='lhc.ltr',plt_dir='.',verbose=True):
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
    if verbose: print '... getting data from %s'%ltr_file
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
        if verbose: print 'deleted file %s'%(ff)
# create new output files
    if verbose: print '... calling script %s/lhcpost'%(script_path)
    p = Popen(['%s/lhcpost'%(script_path),ltr_dir,ltr_file,script_path], stdout=PIPE)
    stdout,stderr=p.communicate()
    if verbose: print stdout
    if stderr != None:
      print "ERROR while executing command lhcpost %s %s:"%(ltr_dir,ltr_file)
      print stderr
      return
    check=True
    for ff in 'emit.txt intensity.txt lossrate.txt luminosity.txt'.split():
      if not os.path.isfile(os.path.join(ltr_dir,ff)):
        print "ERROR when calling output.sh: file %s has not been generated!"%ff
        check=False
    if check:
      print "... created emit.txt, intensity.txt, lossrate.txt, luminosity.txt"
# -- get input parameters
    lout=stdout.split('\n')
    input_param={}
    for s in lout:
      if s.find('Monitoring IP')>-1:
        input_param['Monitoring IP'] = s.split()[-1]
      if s.find('steplen')>-1:
#        input_param['steplen'] = float((s.split(',')[0]).split()[-1])
#        input_param['nsteps']  = float((s.split(',')[1]).split()[-1])
        input_param['steplen'] = [float(i) for i in (s.split()[1][:-1]).replace('(','').replace(')','').split(',')]
        input_param['nsteps']  = float((s.split(',')[-1]).split()[-1])
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
    nstep=len(step)/len(input_param['steplen'])
    nturn            = np.cumsum(np.array(input_param['steplen'] * nstep))
    time             = nturn*input_param['circlhc']/(beta*clight) # time [s]
    sigm             = np.loadtxt('%s/sigm.txt'%(ltr_dir))
    intensity        = np.loadtxt('%s/intensity.txt'%(ltr_dir))
    lumi             = np.loadtxt('%s/luminosity.txt'%(ltr_dir))
    loss             = np.genfromtxt('%s/lossrate.txt'%(ltr_dir),missing_values='************',filling_values=0)
    # for loop over sigm,lumi,loss,intensity didn't work -> fix this
    if len(lumi) == 0:
      lumi = np.array([float('nan') for i in range(len(nturn))])
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
  def get_hist(self,fn='lhc.hist'):
    """read in histogram data from file *fn*
    and store in self.hist with:
    'header' : header data containing
      mean : mean values
      sigm : standard deviation of gaussian fit
      norm : ?
      emit : calculated emittance
      each of these contains the values for
        x,x',y,y',z,dE/E [cm,rad,cm,rad,cm,1]
    """
    ltr_dir=self.lifedeskenv['ltr_dir']
    fna=os.path.join(ltr_dir,fn) #absolute path
    data = np.array([])
    header = {}
    if not os.path.isfile(fna) and self.hist=={}:
      print "ERROR: file %s does not exist!"%fn
    else:
      print "... getting histogram data from %s"%fn
      #read in histogram data
      print "reading data"
      ftype=[('val','f8'),('x','f8'),('px','f8'),('ax','f8'),('y','f8'),('py','f8'),('ay','f8'),('z','f8'),('pz','f8'),('az','f8'),('parity','S2'),('step','S100')]
      data = np.loadtxt(fna,comments='#',dtype=ftype)
      # read in header with statistical data
      print "reading header with statistical data"
      ff=open(fna,'r')
      header={}
      htype=[('mean','6f8'),('sigm','6f8'),('norm','6f8'),('emit','6f8')]
      for line in ff:
        if line.startswith('#'):
# go in exactly the order of the header, so that values are grouped correctly togeter
          if 'Step' in line:
            step=int(line.split()[1])
          if '|mean' in line: # otherwise is take 'Mean_values'
            mean=map(float,line.replace('|','').split()[1:])
          if '|sigm' in line: # otherwise is take 'Mean_values'
            sigm=map(float,line.replace('|','').split()[1:])
          if '|norm' in line: # otherwise is take 'Mean_values'
            norm=map(float,line.replace('|','').split()[1:])
          if '|emit' in line: # otherwise is take 'Mean_values'
            emit=map(float,line.replace('|','').replace(')','').replace('(','').split()[1:])
          if 'Val' in line: # end of header
# figure out read in of structured array!!!
            header[step]=np.array((mean,sigm,norm,emit),dtype=htype)
      self.hist['header']= header
      self.hist['data']  = data
  def plot_hist(self,fn='lhc.hist',plane='x',nstep=None,fit=True,log=True,res=True,verbose=False):
    """plot histogram of particle distribution.
    Plotrange is by default set to [-6,6] sigma,
    while histograms usually extend further.
    Parameters:
    fn:    filename with histogram data
    plane: plane, options are x,px,ax,y,py,ay,z,pz,az with
      ax,ay,az being the normalized amplitude
    nstep: list of steps for which the histogram is plotted,
      e.g. [1,5,10]. If nturn=None, then the first and last 
      turn are used
    fit: plot Gaussian fit
    res: if res=True the difference (residual)
         between the first histograms (n0) and
         the other histograms (n1,n2,...) is plotted
         in a separte subplot. The list of histogams
         is given by nstep=[n0,n1,...].
    """
    if self.hist == {}:
      print 'ERROR: no histogram data available! Run self.get_hist() to read in histogram data, if it exists.'
    else:
      # map plane to array in self.hist['header'][*]
      hist_keys={'x':0,'px':1,'y':2,'py':3,'z':4,'pz':5}
      lbl_keys={'x':'x','px':'p_x','y':'y','py':'p_y','z':'z','pz':'p_z'}
      if nstep == None:
        steps=self.hist['header'].keys()
        nstep=[steps[0],steps[-1]]
      if res and len(nstep)>1: 
        gs = gridspec.GridSpec(3, 1, height_ratios=[3,1,1])
        ax1 = pl.subplot(gs[1])
        ax2 = pl.subplot(gs[2])
      else: gs = gridspec.GridSpec(1, 1)
      ax0 = pl.subplot(gs[0])
      for step in nstep:
        data = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%step]
        steplen = self.input_param['steplen'][0] # number of turns per step
        fit = (self.hist['header'][step]) 
        width = data['val'][1]-data['val'][0]
        for p in list(plane):
          myplot=pl.plot(data['val'],data[p],ls='steps',label=r'$%s(\mathrm{turn \ %s)}$'%(lbl_keys[p],int(step*steplen)))
          c=myplot[-1].get_color()
          if verbose: print fit['norm'][hist_keys[p]]
          ax0.plot(data['val'],gaussian(data['val'],0,fit['norm'][hist_keys[p]]),linestyle='--',color=c)
          if step != nstep[0]:
            # intersept the histogram bins, plot only bins for which interseption exists
            data_step0 = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%nstep[0]]
            vals_intersept = np.array(list(set(data['val']) & set(data_step0['val'])))
            vals_step0 = np.array([ x in vals_intersept for x in data_step0['val'] ])
            vals = np.array([ x in vals_intersept for x in data['val'] ])
            if np.max(np.abs(data_step0['val'][vals_step0]-data['val'][vals]))>0:
              print 'ERROR: something went wrong! The bins of the data for step %s and %s do not agree'%(nstep[0],step)
            else:
              bins = data_step0['val'][vals_step0]
              residual = (data[p][vals]-data_step0[p][vals_step0])*100
              residual_log = data[p][vals]/data_step0[p][vals_step0]
              ax1.plot(bins,residual,linestyle='-',color=c,label=r'v=$%s$'%lbl_keys[p])
              ax2.plot(bins,residual_log,linestyle='-',color=c,label=r'v=$%s$'%lbl_keys[p])
          if res: axes = [ax0,ax1,ax2]
          else: axes = [ax0]
          for ax in axes:
            ax.set_xlim([-6,6])
            ax.set_xlabel(r'$\sigma$',fontsize=12)
            ax.grid(b=True)
          ax0.legend(loc='upper right',fontsize=12)
          ax0.set_ylabel(r'count')
          if log == True: ax0.set_yscale('log')
          if res:
            ax1.set_title(r'$\mathrm{residual = v(turn \ %s) - v(turn \ %s)}$'%(int(step*steplen),int(nstep[0]*steplen)))
            ax2.set_title(r'$\mathrm{residual = v(turn \ %s)/v(turn \ %s)}$'%(int(step*steplen),int(nstep[0]*steplen)))
            ax1.legend(loc='upper right',fontsize=12,ncol=2,columnspacing=0.4,handlelength=0.1)
            ax2.legend(loc='upper center',fontsize=12,ncol=2,columnspacing=0.4,handlelength=0.1)
            ax1.set_ylabel(r'residual [%]')
            ax2.set_ylabel(r'residual')
            ax2.set_yscale('log')
    pl.tight_layout()
  def mk_hist_video(self,fn='lhc.hist',nstep=None,fit=True,log=True,res=True,plt_dir=None,export=False,verbose=False,ylimhist=[1.e-4,1.1],ylimresdiff=[-15,15],ylimresrat=[1.e-1,15],delay=20):
    """make an animated gif of histograms for all
    planes with x,px in one plot, y,py in one etc..
    Plotrange is by default set to [-6,6] sigma,
    while histograms usually extend further.
    Parameters:
    fn:    filename with histogram data
    nstep: nstep = [start,stop], e.g. [1,5]. If nturn=None, 
      then all steps are plotted.
    fit: plot Gaussian fit
    res: if res=True the difference (residual)
         between the first histograms (n0) and
         the other histograms (n1,n2,...) is plotted
         in a separte subplot. The list of histogams
         is given by nstep=[n0,n1,...].
    export: If True do not delete png files
    ylimhist: ylim for histogram
    ylimresdiff: ylim for plot of residuals v(t)-v(t0)
    ylimresrat: ylim for plot of residuals v(t)/v(0)
    delay: delay for frames for convert 
    """
    if nstep == None: nstep = [1,self.input_param['nsteps']]
    if plt_dir == None: plt_dir = self.lifedeskenv['plt_dir']
    tmpdir=os.path.join(plt_dir,'tmp')
    if os.path.exists(tmpdir) == False: os.makedirs(tmpdir)
    for pp in [['x','px'],['y','py'],['z','pz']]:
      for step in range(nstep[0],nstep[1]+1): 
        pl.figure(pp[0],figsize=(8,8))
        pl.clf()
        self.plot_hist(nstep=[nstep[0],step],plane=pp,fit=fit,log=log,res=res)
        ax2,ax1,ax0=pl.gcf().get_axes()
        ax0.set_ylim(ylimhist)
        ax1.set_ylim(ylimresrat)
        ax2.set_ylim(ylimresdiff)
        fnpl=os.path.join(tmpdir,'%s_%s_%s.png'%(pp[0],pp[1],str(int(step)).zfill(len(str(pp[1]))+1)))
        pl.savefig(fnpl)
        if verbose: print '... save png %s'%(fnpl)
      cmd="convert -delay %s %s %s"%(delay,os.path.join(tmpdir,'%s*.png'%pp[0]),os.path.join(plt_dir,'%s_%s.gif'%(pp[0],pp[1])))
      os.system(cmd)
      if verbose: print '... creating .gif file with convert' 
    if export == False: shutil.rmtree(tmpdir)
  def plot_2d(self,xaxis='time',yaxis='emit1',color='b',lbl=None,title=None,alpha=1.0,linestyle='-',marker='o',indstep=None,verbose=False):
    """plot *xaxis* vs *yaxis*
    Parameters:
    -----------
    param   : hor. (mode=1) or vert. (mode=2) emittance
    tunit   : 'nturn': number of turns
              'time' : time [s]
    norm    : eps = beta*gamma*epsn
              false: rms emittance eps [mum]
              true : normalized emittance epsn [mum]
    indstep : make plots for each steplength, e.g. if
              different stplength are used, create a 
              different plot for each steplength
    color   : plot color
    """
    data_names = self.data.dtype.names
    # check if fields exist
    for n in xaxis,yaxis:
      if n not in data_names:
        print '%s not found in data'%n
        return 0
    x,y = self.data[xaxis],self.data[yaxis]
    pl.plot(x,y,linestyle=linestyle,marker=marker,color=color,label=lbl,alpha=alpha)
    pl.xlabel(r'%s %s'%self._unit[xaxis])
    pl.ylabel(r'%s %s'%self._unit[yaxis])
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
#    # legend on right side
#    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=12)
#    pl.subplots_adjust(left=0.15, right=0.8, top=0.1, bottom=0.1)
    pl.grid()
  def plot_all(self,color='b',lbl=None,title=None,export=None,alpha=1.0,linestyle='-',marker='o'):
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
        self.plot_2d(xaxis='time',yaxis=p,color=color,lbl=lbl,title=title,alpha=alpha,linestyle=linestyle)
        if export != None:
          print '%s.%s'%(p,export)
          pl.savefig('%s/%s.%s'%(self.lifedeskenv['plt_dir'],p,export),bbox_inches='tight')
      except KeyError:
        print 'ERROR in plot_all: could not plot %s'%p
        print '   self.data[%s][0]=%s'%(p,self.data[p][0])
