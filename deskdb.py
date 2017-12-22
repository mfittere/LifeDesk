import os
import sys
from subprocess import Popen,PIPE
from inspect import getfile
from matplotlib import gridspec
import shutil
from utilities import grep
from scipy.stats import norm

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

import numpy.lib.recfunctions as rfn
from matplotlib.colors import LogNorm

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
    loss: amplitudes and at which aperture lost (only
          filled if self.getloss() is called)
      'dist' : input distribution used
      'init' : structured array with
          losses in first 2 turns (adjustment
          of distribution)
      'all' : structured array with losses after 
          first 2 turns
    """
    ltr_dir=os.path.abspath(ltr_dir)
    if not os.path.isdir(ltr_dir):
      print('ERROR: Directory %s not found!'%ltr_dir)
      return
    self.lifedeskenv={'ltr_dir':ltr_dir,'ltr_file':ltr_file,'plt_dir':plt_dir}
    self.input_param=input_param
    self.data = data
    self.hist = {}
    self.loss = {}
# dictionary which stores the lbls and units
    # unnormalized
    self._unit = {'step': ('step number',''),'nturn': ('number of turns',''),'time': ('time','[s]'),'emit1':('hor. emittance','[$\mu$m]'),'emit2':('vert. emittance','[$\mu$m]'),'sigm':('bunch length','[cm]'),'intensity': ('normalized beam intensity',''),'lossrate':('normalized loss rate',''),'luminosity':('luminosity','')}
    # normalized to the initial value
    self._unitnorm = {'emit1':('hor. emittance $\epsilon_x/\epsilon_{x,0}$',''),'emit2':('vert. emittance $\epsilon_y/\epsilon_{y,0}$',''),'sigm':(r'bunch length $\sigma/\sigma_0$','')}
    for c in 'x y z'.split():
      self._unit[c]=('initial %s normalized'%(c),'[$\sigma$]')
      self._unit['p'+c]=('initial $p_%s$ normalized'%(c),'[$\sigma$]')
    for c in 'ax ay az'.split():
      plane=c.split('a')[-1]
      self._unit[c]=('initial normalized amplitude in (%s,p%s)'%(plane,plane),'[$\sigma$]')
    self._unit['ar']=('initial radius $r=\sqrt{a_x^2+a_y^2}$','[$\sigma$]')
  @classmethod
  def getdata(cls,ltr_dir='.',ltr_file='lhc.ltr',plt_dir='.',verbose=False):
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
    outputfiles='emit.txt intensity.txt lossrate.txt luminosity.txt'
# check if outputfiles exist and delete old output files if force = True
    check=True
    for ff in outputfiles.split():
      if os.path.isfile(os.path.join(ltr_dir,ff)):
        os.remove(os.path.join(ltr_dir,ff))
        if verbose: print 'deleted file %s'%(ff)
# create new output files
      if verbose: print '... calling script %s/lhcpost'%(script_path)
      fdevnull=open(os.devnull,'wb')
      p = Popen(['%s/lhcpost'%(script_path),ltr_dir,ltr_file,script_path], stdout=PIPE,stderr=fdevnull)
      fdevnull.close()
      stdout,stderr=p.communicate()
      if verbose: print stdout
      if stderr != None:
        print "ERROR while executing command lhcpost %s %s:"%(ltr_dir,ltr_file)
        print stderr
        return
      check=True
      for ff in outputfiles.split():
        if not os.path.isfile(os.path.join(ltr_dir,ff)):
          print "ERROR when calling output.sh: file %s has not been generated!"%ff
          check=False
      if check and verbose:
        print '... created %s'%outputfiles
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
      lumi = np.array([0 for i in range(len(nturn))])
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
  def get_hist(self,fn='lhc.hist',verbose=False):
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
      if verbose: print "ERROR: file %s does not exist!"%fn
    else:
      if verbose: print "... getting histogram data from %s"%fn
      #read in histogram data
      if verbose: print "reading data"
      ftype=[('val','f8'),('x','f8'),('px','f8'),('ax','f8'),('y','f8'),('py','f8'),('ay','f8'),('z','f8'),('pz','f8'),('az','f8'),('parity','S2'),('step','S100')]
      data = np.loadtxt(fna,comments='#',dtype=ftype)
      # read in header with statistical data
      if verbose: print "reading header with statistical data"
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
  def plot_hist_xyz(self,fn='lhc.hist',nstep=None,fit=True,
        log=True,ylimhist=[1.e-4,1.1],ylimresdiff=[-15,15],ylimresrat=[1.e-1,15],
        color='k',title=''):
    """plot histogram of particle distribution
    in physical coordinates for start and end turn.
    For the longitudinal the 3 step is used.
    Plotrange is by default set to [-6,6] sigma,
    while histograms usually extend further.
    Parameters:
    fn:    filename with histogram data
    fit: plot Gaussian fit
    log: switch on log scale for distribution plot
    nstep: either list with [nstart,nend] or 
           nturn=None, then the first and last 
           turn are used
    color: color used for plot
    title: figure title
    """
    if self.hist == {}:
      raise ValueError('ERROR: no histogram data available! Run '+
        'self.get_hist() to read in histogram data, if it exists.')
    # map plane to array in self.hist['header'][*]
    hist_keys={'x':0,'px':1,'y':2,'py':3,'z':4,'pz':5}
    lbl_keys={'x':'x','px':'p_x','y':'y','py':'p_y','z':'z','pz':'p_z','ax':'ax','ay':'ay','r':'r'}
    if nstep == None:
      steps=self.hist['header'].keys()
      nstep=[steps[0],steps[-1]]
    for plane in 'x','y','z':
      # generate the figure
      pl.figure('%s %s'%(title,plane),figsize=(5,6))
      pl.gcf().set_size_inches(5, 6, forward=True)
      gs = gridspec.GridSpec(3, 1, height_ratios=[2,1,1],hspace=0.07,
            left=0.15,right=0.98)
      ax0 = pl.subplot(gs[0])
      ax1 = pl.subplot(gs[1])
      ax2 = pl.subplot(gs[2])
      for step,alpha in zip(nstep,[0.2,1.0]):
        data = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%step]
        steplen = self.input_param['steplen'][0] # number of turns per step
        histfit = (self.hist['header'][step]) 
        width = data['val'][1]-data['val'][0]
        # get the data
        xdata = data['val']
        ydata = data[plane]
        myplot=ax0.plot(xdata,ydata,ls='steps',color=color[0],alpha=alpha,
                       label=r'$%s(\mathrm{turn \ %4.1e)}$'%(lbl_keys[plane],int(step*steplen)))
        if fit:
          ax0.plot(xdata,gaussian(xdata,0,histfit['norm'][hist_keys[plane]]),linestyle='-',color='r',label='Gaussian fit(turn %4.1e)'%(int(step*steplen)),alpha=alpha)
          ax1.plot(xdata,(ydata-gaussian(xdata,0,histfit['norm'][hist_keys[plane]]))*100,linestyle='-',color='r',alpha=alpha)
      # intersept the histogram bins, plot only bins for which interseption exists
      data_step0 = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%nstep[0]]
      xdata0=data_step0['val']
      ydata0 = data_step0[plane]
      data = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%nstep[1]]
      xdata = data['val']
      ydata = data[plane]
      vals_intersept = np.array(list(set(xdata) & set(xdata0)))
      vals_step0 = np.array([ x in vals_intersept for x in xdata0 ])
      vals = np.array([ x in vals_intersept for x in xdata ])
      if np.max(np.abs(xdata0[vals_step0]-xdata[vals]))>0:
        raise ValueError('Something went wrong! The bins of the data for step %s and %s do not agree'%(nstep[0],step))
      bins = xdata0[vals_step0]
      residual = (ydata[vals]-ydata0[vals_step0])*100
      ratio = ydata[vals]/ydata0[vals_step0]
      ax1.plot(bins,residual,linestyle='-',color=color)
      ax2.plot(bins,ratio,linestyle='-',color=color)
      for ax in [ax0,ax1,ax2]:
        ax.set_xlim([-6,6])
        ax0.legend(loc='lower center',fontsize=12)
        ax.grid(b=True)
      # remove tick labels from ax0,ax1 + add plot label to ax2
      ax0.set_xticklabels([])
      ax1.set_xticklabels([])
      ax2.set_xlabel(r'$\sigma$')
      ax0.set_ylabel(r'count')
      if log == True: ax0.set_yscale('log')
      ax1.legend(loc='lower left',fontsize=12,ncol=2,columnspacing=0.4,handlelength=0.1)
      ax2.legend(loc='upper center',fontsize=12,ncol=2,columnspacing=0.4,handlelength=0.1)
      ax1.set_ylabel('residual [\%]')
      ax2.set_ylabel(r'ratio')
      ax2.set_yscale('log')
      ax0.set_ylim(ylimhist)
      ax1.set_ylim(ylimresdiff)
      ax2.set_ylim(ylimresrat)
  def plot_hist(self,fn='lhc.hist',plane='x',nstep=None,fit=True,
        log=True,res=True,
        ylimhist=[1.e-4,1.1],ylimresdiff=[-15,15],ylimresrat=[1.e-1,15],
        color=None,verbose=False):
    """plot histogram of particle distribution.
    Plotrange is by default set to [-6,6] sigma,
    while histograms usually extend further.
    Parameters:
    fn:    filename with histogram data
    plane: plane, options are
             1) physical coordinates: x,px,y,py,z,pz
             2) normalized amplitudes: ax,ay,az
             3) radius: r=sqrt(ax**2+ay**2)
    nstep: list of steps for which the histogram is plotted,
      e.g. [1,5,10]. If nturn=None, then the first and last 
      turn are used
    fit: plot Gaussian fit
    res: if res=True the difference (residual)
         and the ratio (ratio) between the
         first histograms (n0) and the other 
         histograms (n1,n2,...) is plotted
         in a separte subplot. The list of histogams
         is given by nstep=[n0,n1,...].
    color: color of plot, either string or list of strings 
           with 2*len(plane)
    """
    if self.hist == {}:
      raise ValueError('ERROR: no histogram data available! Run '+
        'self.get_hist() to read in histogram data, if it exists.')
#    # generate the figuer
    pl.gcf().set_size_inches(5, 6, forward=True)
    # map plane to array in self.hist['header'][*]
    hist_keys={'x':0,'px':1,'y':2,'py':3,'z':4,'pz':5}
    lbl_keys={'x':'x','px':'p_x','y':'y','py':'p_y','z':'z','pz':'p_z','ax':'ax','ay':'ay','r':'r'}
    if nstep == None:
      steps=self.hist['header'].keys()
      nstep=[steps[0],steps[-1]]
    if res and len(nstep)>1: 
      gs = gridspec.GridSpec(3, 1, height_ratios=[2,1,1],hspace=0.07,
            left=0.15,right=0.98)
      ax1 = pl.subplot(gs[1])
      ax2 = pl.subplot(gs[2])
    else: gs = gridspec.GridSpec(1, 1)
    ax0 = pl.subplot(gs[0])
    if color is None or type(color) is str:
      color = [color]*len(plane)*2
    if len(color) != len(plane)*2:
      raise ValueError('color must be either None, string or length' +
        'of list of colors must be 2 x length of planes!')
    for step,al,lns in nstep:
      data = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%step]
      steplen = self.input_param['steplen'][0] # number of turns per step
      fit = (self.hist['header'][step]) 
      width = data['val'][1]-data['val'][0]
      if type(plane) is str:
        lplane = [plane]
      else:
        lplane = plane
      for p,cl in zip(lplane,color):
        # get the data
        xdata = data['val']
        if p == 'r':
          ydata = np.sqrt(data['ax']**2+data['ay']**2)
        else:
          ydata = data[p]
        myplot=ax0.plot(xdata,ydata,ls='steps',color=cl,
                       label=r'$%s(\mathrm{turn \ %s)}$'%(lbl_keys[p],int(step*steplen)))
        c=myplot[-1].get_color()
        if p in ['x','px','y','py','z','pz']:
          ax0.plot(xdata,gaussian(xdata,0,fit['norm'][hist_keys[p]]),linestyle='--',color=c)
          if verbose: print fit['norm'][hist_keys[p]]
        if step != nstep[0]:
          # plot difference in respect to Gaussian
          if p in ['x','px','y','py','z','pz']:
            ax1.plot(xdata,(ydata-gaussian(xdata,0,fit['norm'][hist_keys[p]]))*100,linestyle='--',color=c)
          # intersept the histogram bins, plot only bins for which interseption exists
          data_step0 = (self.hist['data'])[(self.hist['data'])['step']=='S_%s'%nstep[0]]
          xdata0=data_step0['val']
          if p == 'r':
            ydata0 = np.sqrt(data_step0['ax']**2+data_step0['ay']**2)
          else:
            ydata0 = data_step0[p]
          vals_intersept = np.array(list(set(xdata) & set(xdata0)))
          vals_step0 = np.array([ x in vals_intersept for x in xdata0 ])
          vals = np.array([ x in vals_intersept for x in xdata ])
          if np.max(np.abs(xdata0[vals_step0]-xdata[vals]))>0:
            print 'ERROR: something went wrong! The bins of the data for step %s and %s do not agree'%(nstep[0],step)
          else:
            bins = xdata0[vals_step0]
            residual = (ydata[vals]-ydata0[vals_step0])*100
            ratio = ydata[vals]/ydata0[vals_step0]
            ax1.plot(bins,residual,linestyle='-',color=c,label=r'$%s$'%lbl_keys[p])
            ax2.plot(bins,ratio,linestyle='-',color=c,label=r'$%s$'%lbl_keys[p])
        if res: axes = [ax0,ax1,ax2]
        else: axes = [ax0]
        ftl = 12 # fontsize legend
        fta = 12
        for ax in axes:
          if p in ['r','ax','ay','az']:
            ax.set_xlim([0,6])
            ax0.legend(loc='lower left',fontsize=ftl)
          else:
            ax.set_xlim([-6,6])
            ax0.legend(loc='lower center',fontsize=ftl)
          ax.grid(b=True)
        # remove tick labels from ax0,ax1 + add plot label to ax2
        ax0.set_xticklabels([])
        ax1.set_xticklabels([])
        ax2.set_xlabel(r'$\sigma$',fontsize=fta)
        ax0.set_ylabel(r'count',fontsize=fta)
        if log == True: ax0.set_yscale('log')
        if res:
#          ax1.set_title(r'$\mathrm{residual = v(turn \ %s) - v(turn \ %s)}$'%(int(step*steplen),int(nstep[0]*steplen)))
#          ax2.set_title(r'$\mathrm{ratio = v(turn \ %s)/v(turn \ %s)}$'%(int(step*steplen),int(nstep[0]*steplen)))
          ax1.legend(loc='lower left',fontsize=ftl,ncol=2,columnspacing=0.4,handlelength=0.1)
          ax2.legend(loc='upper center',fontsize=ftl,ncol=2,columnspacing=0.4,handlelength=0.1)
          ax1.set_ylabel('residual [\%]',fontsize=fta)
          ax2.set_ylabel(r'ratio',fontsize=fta)
          ax2.set_yscale('log',fontsize=fta)
          ax0.set_ylim(ylimhist)
          ax1.set_ylim(ylimresdiff)
          ax2.set_ylim(ylimresrat)
  def mk_hist_video(self,fn='lhc.hist',nstep=None,fit=True,
        log=True,res=True,plt_dir=None,export=False,verbose=False,
        ylimhist=[1.e-4,1.1],ylimresdiff=[-15,15],ylimresrat=[1.e-1,15],
        delay=20):
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
        self.plot_hist(nstep=[nstep[0],step],plane=pp,fit=fit,log=log,
                       res=res)
        ax2,ax1,ax0=pl.gcf().get_axes()
        ax0.set_ylim(ylimhist)
        ax1.set_ylim(ylimresrat)
        ax2.set_ylim(ylimresdiff)
        fnpl=os.path.join(tmpdir,'%s_%s_%s.png'%(pp[0],pp[1],
               str(int(step)).zfill(len(str(pp[1]))+1)))
        pl.savefig(fnpl)
        if verbose: print '... save png %s'%(fnpl)
      cmd="convert -delay %s %s %s"%(delay,
          os.path.join(tmpdir,'%s*.png'%pp[0]),
          os.path.join(plt_dir,'%s_%s.gif'%(pp[0],pp[1])))
      os.system(cmd)
      if verbose: print '... creating .gif file with convert' 
    if export == False: shutil.rmtree(tmpdir)
  def getloss(self,fndist=None,verbose=False):
    """get amplitudes and aperture for lost
    particles
    
    Parameters
    ----------
    fndist: path to input distribution
            if fndist = None (default) then the
            distribution name is taken from the
            ltr file and loaded from the distributions
            folder in LifeDesk
    verbose: verbose mode
    
    Returns:
    --------
    dictionary self.loss with
    'dist': input distribution
    'init': losses in first 2 turns
    'all': losses after first two turns

    parameters of structured array:
    x,px,y,py,z,pz: normalized coordinates in sigma
    weight: weight in distribution
    turn: turn number when particle got lost
    why: at which aperture it got lost, e.g. AX or AY
    """
    lifedeskdir=os.path.dirname(getfile(LifeDeskDB))
    inputfile=os.path.join(self.lifedeskenv['ltr_dir'],self.lifedeskenv['ltr_file'])
    # get input distribution file, check that only one is defined
    if fndist is None:
      fndist=grep('Distr_init',inputfile)
      if len(fndist)==0:
        print 'ERROR: no input distribution found!'
        return
      elif len(fndist)>1:
        print 'ERROR: more than one input distribution found in %s'%fn
        return
      else:
        fndist=(fndist[0].rstrip().split('/'))[-1]
        if verbose: print '... input distribution used %s'%fndist
      # get path to distribution and inilost.sh files in LifeDesk directory
      fndist='%s/distributions/%s'%(lifedeskdir,fndist)
    if not os.path.isfile(fndist):
      print('ERROR: distribution file %s not found!'%fndist)
      return
    script_path='%s/scripts'%(lifedeskdir)
    if verbose:
      print "... calling script %s/inilost.sh %s %s"%(script_path,
        inputfile,fndist)
    # go to study directory
    cwd = os.getcwd()
    os.chdir(self.lifedeskenv['ltr_dir'])
    p = Popen(['%s/inilost.sh'%(script_path),inputfile,fndist], stdout=PIPE,stderr=PIPE)
    stdout,stderr=p.communicate()
    os.chdir(cwd)
    if verbose:
      print stdout
    if stderr != None and 'grep' not in stderr:
      print "ERROR while executing command %s/inilost.sh %s %s"%(script_path,inputfile,fndist)
      print stderr
      return
    # path to input distribution
    self.loss['dist']=fndist
    # get the losses and amplitudes
    ftype=[('x','f8'),('px','f8'),('y','f8'),('py','f8'),('z','f8'),('pz','f8'),('weight','f8'),('nturn','f8'),('why','S100')]
    for fn,par in [('inilost.1.out','init'),('inilost.out','all')]:
      # get data
      data1 = np.loadtxt(os.path.join(self.lifedeskenv['ltr_dir'],fn),comments='#',dtype=ftype)
      # calculate time and amplitue
      clight = 299792458 # speed of light
      gamma  = self.input_param['gamma']
      beta   = np.sqrt(1-1/gamma**2)
      time   = data1['nturn']*self.input_param['circlhc']/(beta*clight) # time [s]
      ax     = np.sqrt(data1['x']**2 + data1['px']**2)
      ay     = np.sqrt(data1['y']**2 + data1['py']**2)
      az     = np.sqrt(data1['z']**2 + data1['pz']**2)
      ar     = np.sqrt(ax**2 + ay**2)
      data2=np.array(zip(time,ax,ay,az,ar),dtype=[('time','f8'),('ax','f8'),('ay','f8'),('az','f8'),('ar','f8')])
      self.loss[par] = rfn.merge_arrays((data1, data2), asrecarray=True, flatten=True)
  def getlifetime(self,sigma=1,verbose=False):
    """
    calculate the lifetime. fndist is assumed to be
    uniform in x,y and xp=yp=0. The lifetime is caluclated
    for a Gaussian distribution with sigma *sigma*
    by folding the uniform distribution with the Gaussian
    Parameter:
    ----------
    sigma: sigma of Gaussian distribution
    
    Returns:
    --------
    tau: lifetime [h]
    """
    lifedeskdir=os.path.dirname(getfile(LifeDeskDB))
    inputfile=os.path.join(self.lifedeskenv['ltr_dir'],self.lifedeskenv['ltr_file'])
    fndist=self.loss['dist']
    # format x,xp,y,yp,z,delta,weight,particle number
    dist0=np.loadtxt(fndist,skiprows=1)
    if np.any(dist0[:,1]) !=0:
      print('ERROR: xp must be 0!')
    if np.any(dist0[:,3]) !=0:
      print('ERROR: yp must be 0!')
    # the number of paricles is proportional to 
    # N0    = sum(x*gauss(x,sigma)*y*gauss(y,sigma))
    # Nlost(t) = sum(x_lost*gauss(x_lost,sigma)*y_lost*gauss(y_lost,sigma))
    # => N(t)=N0-Nlost(t)
    # the lifetime can be obtained with:
    # tau = t/(ln(N)-ln(N0))
    x0,y0=dist0[:,0],dist0[:,2]
    xlost,ylost = self.loss['all']['x'],self.loss['all']['y']
    N0=np.sum(norm.pdf(x0,0,sigma)*norm.pdf(y0,0,sigma))
    Nlost = np.sum(norm.pdf(xlost,0,sigma)*norm.pdf(ylost,0,sigma))
    N = N0-Nlost
    # get the time interval
    dt = self.data['time'][-1]-self.data['time'][0]
    return dt/(np.log(N0)-np.log(N))/(60*60)
  def plot_loss_2d(self,xaxis='ax',yaxis='time',bins=50,log=True):
    """make a 2d histogram of losses"""
    if self.loss == {}:
      print('Losses have not been generated or failed. Try self.getloss()!')
      return 'Failed'
    if len(self.loss['all'])==0:
      print('WARNING: No lost particles. Aborting!')
      return 'Failed'
    if log:
      pl.hist2d(self.loss['all'][xaxis],self.loss['all'][yaxis],bins=bins,norm=LogNorm())
    else:
      pl.hist2d(self.loss['all'][xaxis],self.loss['all'][yaxis],bins=bins)
    clb = pl.colorbar()
    clb.set_label('$\log(N_{\mathrm{lost}})$',fontsize=14)
    pl.xlabel(r'%s %s'%self._unit[xaxis])
    pl.ylabel(r'%s %s'%self._unit[yaxis])
  def plot_2d(self,xaxis='time',yaxis='emit1',norm=False,color=None,lbl=None,title=None,alpha=1.0,linestyle='-',marker='o',indstep=None,verbose=False):
    """plot *xaxis* vs *yaxis*
    Parameters:
    -----------
    param   : hor. (mode=1) or vert. (mode=2) emittance
    tunit   : 'nturn': number of turns
              'time' : time [s]
    norm    : normalize variable to initial value
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
    if norm: y=y/y[0]
    pl.plot(x,y,linestyle=linestyle,marker=marker,color=color,label=lbl,alpha=alpha)
    pl.xlabel(r'%s %s'%self._unit[xaxis])
    if norm: pl.ylabel(r'%s %s'%self._unitnorm[yaxis])
    else: pl.ylabel(r'%s %s'%self._unit[yaxis])
    # place the legend to the outside of the plot
    box = pl.gca().get_position()
    # if more than 4 entries, take two columns, otherwise one
    nlabels = len(pl.gca().get_legend_handles_labels()[1])
    if (nlabels <=4) :
      ncol = 1
    elif (nlabels >= 5 and nlabels <=8):
      ncol = 2
    else:
      ncol = 3
    # legend on top
    pl.legend(bbox_to_anchor=(0., 1.07, 1.0, .102), loc=3,ncol=ncol, mode="expand", borderaxespad=0.,fontsize=12,title=title)
    pl.subplots_adjust(left=0.15, right=0.95, top=0.68, bottom=0.1)
    pl.grid(b=True)

  def plot_all(self,color=None,lbl=None,title=None,export=None,alpha=1.0,linestyle='-',marker='o',fignum=''):
    """plots emittance, bunch length, intensity
    luminosity and loss rate vs time [s].
    
    Parameters:
    -----------
    color: plot color
    lbl: plot label
    title: title of legend, e.g. collision, beta*=6m
    export: export format, e.g. 'png'
    fignum: figure number (string) appended to the usual emit1,emit2 etc
    """
    for p in ['emit1','emit2','sigm','intensity','lossrate','luminosity']:
      if not np.isnan(self.data[p][0]):
        pl.figure(p+fignum,figsize=(8,6))
        try:
          self.plot_2d(xaxis='time',yaxis=p,color=color,marker=marker,lbl=lbl,title=title,alpha=alpha,linestyle=linestyle)
          if export != None:
            pl.savefig('%s/%s.%s'%(self.lifedeskenv['plt_dir'],p,export),bbox_inches='tight')
        except KeyError:
          print 'ERROR in plot_all 2: could not plot %s'%p
          print '   self.data[%s][0]=%s'%(p,self.data[p][0])
          pl.close(p)
      else:
        print 'ERROR in plot_all: could not plot %s'%p
        print '   self.data[%s][0]=%s'%(p,self.data[p][0])

