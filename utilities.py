import re as re

def grep(ss,fns):
  """linux grep command with
  uses re for regular expression searching

  grep ss fn
  Parameters:
  -----------
  fns: list of files
  ss: search string
  Returns:
  --------
  exp: list of all lines in which ss was 
       found
  """
  if type(fns) == str: fns = [fns]
  exp=[]
  for fn in fns:
    for line in open(fn,'r'):
      if re.search(ss,line):
        exp.append(line)
  return exp
