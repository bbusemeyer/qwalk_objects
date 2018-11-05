from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.core.periodic_table import Element
import qwalk_objects as obj
from xml.etree.ElementTree import ElementTree
import numpy as np
import os

# TODO the code revolving around _elements is bad and should be replaced.
# * _elements is a hidden dependency between geom and basis_section
# * new geom functions must include this dependency, but its not really documented how.

class CrystalWriter: 
  def __init__(self,options):
    #Geometry input.
    self.struct=None
    self.struct_input=None # Manual structure input.
    self.supercell=None #[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]  
    self.boundary='3d'

    #Electron model
    self.spin_polarized=True    
    self.xml_name="BFD_Library.xml"
    self.functional={'exchange':'PBE','correlation':'PBE','hybrid':0,'predefined':None}
    self.total_spin=0
    self.modisymm=None
    self.symmremo=False # Remove all symmetries for SCF.

    # Initial guess.
    self.initial_spins=[]
    self.guess_fort=None # or, path to a fort.79 or fort.9.
    self.guess_fort13=None # path to fort.13, activates reduced symmetry guess (GUESSYMP)
    self.spinedit=[] # Spins to flip from guess_fort

    # D or F occupation guess. If any is None it's ignored.
    self.majority_guess=None # Guess for the initial orbital occupation. See FDOCCUP
    self.minority_guess=None
    self.df_basis_element=None # must also know the basis order because crystal is a jerk.

    #Numerical convergence parameters
    self.basis_params=[0.2,2,3]
    self.basislines=None # Option to manually input basis lines.
    self.cutoff=0.0    
    self.kmesh=[8,8,8]
    self.gmesh=16
    self.tolinteg=[8,8,8,8,18]
    self.dftgrid=''
    
    #Memory
    self.biposize=100000000
    self.exchsize=10000000
    
    #SCF convergence parameters
    self.initial_charges={}    
    self.fmixing=80
    self.maxcycle=100
    self.edifftol=8
    self.levshift=[]
    self.diis=True
    self.diis_opts = [] # lines for DIIS.
    self.broyden=[]
    self.anderson=False
    self.smear=0.0001

    # Use the new crystal2qmc script. This should change soon!
    self.cryapi=True

    self.restart=False
    self.completed=False

    self.set_options(options)

    # Internal.
    self._elements=[] # List of elements set by geometry reading. 
    
  #-----------------------------------------------
    
  def set_struct_fromcif(self,cifstr,primitive=True):
    self.primitive=primitive
    self.cif=cifstr
    pystruct=CifParser.from_string(self.cif).get_structures(primitive=self.primitive)[0]
    self.struct=pystruct.as_dict()
  #-----------------------------------------------

  def set_struct_fromxyz(self,xyzstr):
    self.xyz=xyzstr
    self.struct=XYZ.from_string(xyzstr).molecule.as_dict()
    self.boundary="0d"

  #-----------------------------------------------

  def set_options(self, d):
    # TODO add check for charge vs spin.
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        raise AssertionError("Error: %s not a keyword for CrystalWriter"%k)
      selfdict[k]=d[k]

    if 'spinedit' in d:
      assert len(d['spinedit'])==0 or ('guess_fort' in d and d['guess_fort'] is not None),\
          "spinedit requires guess_fort."

  #-----------------------------------------------
  def apply_fdocc(self):
    ''' For spin setting, apply the majority and minority occupations for each spin (for FDOCCUP)
    Args: 
      spins (list): atomic spin assignments.
      basis_element (int): position of orbital in the basis. If the d is the 3rd basis set for the atom, this would be 5.
      majocc (list): ints representing relative occpuations of majority channel d or f orbitals. See FDOCCUP for more.
      minocc (list): ints representing relative occpuations of minority channel d or f orbitals. See FDOCCUP for more.
    Returns:
      str: FDOCCS input section.
    '''
    print("Warning: fdocc needs generalization to more than one spin-ful atom.")
    majocc = ' '.join([str(i) for i in self.majority_guess])
    minocc = ' '.join([str(i) for i in self.minority_guess])

    if   len(self.majority_guess) == len(self.minority_guess) == 5: lab = 3
    elif len(self.majority_guess) == len(self.minority_guess) == 5: lab = 4
    else: raise AssertionError('Majority and minority guess must be for d or f.')

    fdoccs = []
    nmod = 0
    for sdex,spin in enumerate(self.initial_spins):
      if spin == 1:
        fdoccs.append('%d %d %d'%(sdex+1,self.df_basis_element,lab))
        fdoccs.append(majocc)
        fdoccs.append(minocc)
        nmod += 1
        break # only need one per species.
      elif spin == -1:
        fdoccs.append('%d %d %d'%(sdex+1,self.df_basis_element,lab))
        fdoccs.append(minocc)
        fdoccs.append(majocc)
        nmod += 1
        break # only need one per species.
      else:
        pass
    fdoccs = ['FDOCCUP',str(nmod)] + fdoccs
    return '\n'.join(fdoccs)

  #-----------------------------------------------
  def crystal_input(self,section4=[]):

    geomlines=self.geom()
    basislines=self.basis_section()

    modisymm=[]
    if self.symmremo:
      modisymm=["SYMMREMO"]
    elif self.modisymm is None:
      if len(self.initial_spins)>0: # auto-MODISYMM
        zro=[i for i,s in enumerate(self.initial_spins) if s== 0]
        ups=[i for i,s in enumerate(self.initial_spins) if s== 1]
        dns=[i for i,s in enumerate(self.initial_spins) if s==-1]
        self.modisymm=[]
        self.modisymm+=[(i,0) for i in zro]
        self.modisymm+=[(i,1) for i in ups]
        self.modisymm+=[(i,2) for i in dns]

        modisymm+=['MODISYMM']
        modisymm+=[str(len(self.modisymm))]
        modisymm+=["%d %d"%(a+1,i+1) for a,i in self.modisymm]
    else: # Manual MODISYMM
      modisymm+=['MODISYMM']
      modisymm+=[str(len(self.modisymm))]
      modisymm+=["%d %d"%(at+1,lab+1) for at,lab in self.modisymm]

    
    outlines = ["Generated by CrystalWriter"] +\
                geomlines +\
                modisymm + \
                ["END"] +\
                basislines +\
                ["99 0"] +\
                ["CHARGED"] +\
                ["END"]

    if self.boundary=="3d": # Can also include 2d later.
      outlines+= ["SHRINK","0 %i"%self.gmesh]
      outlines+= [" ".join(map(str,self.kmesh))]

    outlines+=["DFT"]
    if self.spin_polarized:
      outlines+=["SPIN"]
    if 'predefined' in self.functional.keys() \
        and self.functional['predefined']!=None:
      outlines+=[self.functional['predefined']]
    else:
      outlines += [ 
        "EXCHANGE",
        self.functional['exchange'],
        "CORRELAT",
        self.functional['correlation'],
        "HYBRID", 
        str(self.functional['hybrid'])]
    if self.dftgrid!="":
      outlines+=[self.dftgrid]
    outlines+=["END",
      "SCFDIR",
      "SAVEPRED",
      "BIPOSIZE",
      str(self.biposize),
      "EXCHSIZE",
      str(self.exchsize),
      "TOLDEE",
      str(self.edifftol),
      "TOLINTEG",
      ' '.join(map(str,self.tolinteg)),
      "MAXCYCLE",
      str(self.maxcycle),
      "SAVEWF"
    ]
    if self.smear > 0:
      outlines+=["SMEAR",str(self.smear)]
    if self.spin_polarized:
      outlines+=['SPINLOCK','%d %d'%(self.total_spin,self.maxcycle)]

    if len(self.initial_spins)>0:
      outlines+=['ATOMSPIN',str(len(self.initial_spins))]
      outlines+=["%i %i"%(i+1,s) for i,s in enumerate(self.initial_spins)]
    if len(self.initial_spins)>0 and \
        self.majority_guess is not None and \
        self.minority_guess is not None and \
        self.df_basis_element is not None:
      outlines+=[self.apply_fdocc()]
    if len(self.spinedit)>0:
      outlines+=['SPINEDIT',str(len(self.spinedit))]
      outlines+=[str(i) for i in self.spinedit]

    outlines += [
        "FMIXING",
        str(self.fmixing)
      ]
    if not self.diis:
      outlines+=[
          "NODIIS",
        ]
      if self.anderson:
        outlines+=["ANDERSON"]
      elif len(self.broyden)!=0:
        outlines+=["BROYDEN"," ".join(map(str,self.broyden))]
      elif self.levshift!=[]:
        outlines+=["LEVSHIFT"," ".join(map(str,self.levshift))]
    else:
      outlines += self.diis_opts

    outlines+=section4
    if self.restart:
      outlines+=["GUESSP"]
    elif self.guess_fort13 is not None:
      outlines+=["GUESSYMP"]
    elif self.guess_fort is not None:
      outlines+=["GUESSP"]
    outlines+=["END"]

    return "\n".join(outlines)

  #-----------------------------------------------
  def properties_input(self):
    outlines=['NEWK']
    if self.boundary=='3d':
      outlines+=[ "0 %i"%self.gmesh,
                  " ".join(map(str,self.kmesh))]
    outlines+=["1 0"]
    if self.cryapi:
      outlines+=["CRYAPI_OUT"]
    else:
      outlines+=["67 999"]
    outlines+=["END"]
    return "\n".join(outlines)

  #-----------------------------------------------
  def write_crys_input(self,filename):
    outstr=self.crystal_input()
    with open(filename,'w') as outf:
      outf.write(outstr)
      outf.close()
    self.completed=True

  #-----------------------------------------------
  def write_prop_input(self,filename):
    outstr=self.properties_input()
    with open(filename,'w') as outf:
      outf.write(outstr)
    self.completed=True

  #-----------------------------------------------
  def check_status(self):
    # Could add consistancy check here.
    status='unknown'
    if os.path.isfile(self.cryinpfn) and os.path.isfile(self.propinpfn):
      status='ok'
    else:
      status='not_started'
    return status

########################################################
  def geom(self):
    """Generate the geometry section for CRYSTAL"""
    if self.struct_input is not None:
      # TODO make into function like geom3d?
      geomlines=[
          'CRYSTAL',
          '0 0 0',
          str(self.struct_input['symmetry']),
          ' '.join(map(str,self.struct_input['parameters'])),
          str(len(self.struct_input['coords']))
        ]
      self._elements=set()
      for coord in self.struct_input['coords']:
        geomlines+=[' '.join(map(str,coord))]
        self._elements.add(obj.crystal2qmc.periodic_table[coord[0]-200-1].capitalize()) # Assumes ECP!
      self._elements = sorted(list(self._elements)) # Standardize ordering.
      if self.supercell is not None:
        geomlines+=['SUPERCELL']
        for row in self.supercell:
          geomlines+=[' '.join(map(str,row))]
      return geomlines
    elif self.struct is not None:
      assert self.boundary in ['0d','3d'],"Invalid or not implemented boundary."
      if self.boundary=="3d":
        return self.geom3d()
      elif self.boundary=='0d': 
        return self.geom0d()
      else:
        print("Weird value of self.boundary",self.boundary)
        quit() # This shouldn't happen.
    else:
      raise AssertionError("No geometry input found; set struct or struct_input.")


########################################################

  def geom3d(self):
    lat=self.struct['lattice']
    sites=self.struct['sites']
    self._elements=set()
    for s in self.struct['sites']:
      nm=s['species'][0]['element']
      self._elements.add(nm)
    self._elements = sorted(list(self._elements)) # Standardize ordering.

    geomlines=[
        "CRYSTAL","0 0 0",
        "1",
        '%g %g %g %g %g %g'%(lat['a'],lat['b'],lat['c'],lat['alpha'],lat['beta'],lat['gamma'])
      ]

    geomlines+=["%i"%len(sites)]
    for v in sites:
      nm=v['species'][0]['element']
      nm=str(Element(nm).Z+200)
      geomlines+=[nm+" %g %g %g"%(v['abc'][0],v['abc'][1],v['abc'][2])]

    if self.supercell is not None:
      geomlines+=["SUPERCELL"]
      for row in self.supercell:
        geomlines+=[' '.join(map(str,row))]

    return geomlines

########################################################

  def geom0d(self):
    geomlines=["MOLECULE","1"]
    geomlines+=["%i"%len(self.struct['sites'])]
    self._elements=set()
    for v in self.struct['sites']:
      nm=v['species'][0]['element']
      self._elements.add(nm)
      nm=str(Element(nm).Z+200)
      geomlines+=[nm+" %g %g %g"%(v['xyz'][0],v['xyz'][1],v['xyz'][2])]
    self._elements = sorted(list(self._elements)) # Standardize ordering.

    return geomlines

########################################################

  def basis_section(self):
    if self.basislines is None:
      basislines=[]
      for e in self._elements:
        basislines+=self.generate_basis(e)
    else:
      basislines=self.basislines
      if basislines[-1]=='': basislines=basislines[:-1]
    return basislines


########################################################
  def generate_basis(self,symbol):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com> and Lucas K. Wagner
    Returns a string containing the basis section.  It is modified according to a simple recipe:
    1) The occupied atomic orbitals are kept, with exponents less than 'cutoff' removed.
    2) These atomic orbitals are augmented with uncontracted orbitals according to the formula 
        e_i = params[0]*params[2]**i, where i goes from 0 to params[1]
        These uncontracted orbitals are added for every occupied atomic orbital (s,p for most elements and s,p,d for transition metals)

    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.

    Returns:
        str: The pseudopotential and basis section.


    Uses the following member variables:
        xml_name (str): The name of the XML pseudopotential and basis
            set database.
        cutoff: smallest allowed exponent
        params: parameters for generating the augmenting uncontracted orbitals
        initial_charges    
    """
    
    maxorb=3
    basis_name="vtz"
    nangular={"s":1,"p":1,"d":1,"f":1,"g":0}
    maxcharge={"s":2,"p":6,"d":10,"f":15}
    basis_index={"s":0,"p":2,"d":3,"f":4}
    transition_metals=["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
    if symbol in transition_metals:
      maxorb=4
      nangular['s']=2
    
    tree = ElementTree()
    tree.parse(self.xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    atom_charge = int(element.find('./Effective_core_charge').text)
    if symbol in self.initial_charges.keys():
      atom_charge-=self.initial_charges[symbol]
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
    found_orbitals = []
    totcharge=0
    ret=[]
    ncontract=0
    for contraction in element.findall(basis_path):
        angular = contraction.get('Angular_momentum')
        if found_orbitals.count(angular) >= nangular[angular]:
            continue

        #Figure out which coefficients to print out based on the minimal exponent
        nterms = 0
        basis_part=[]
        for basis_term in contraction.findall('./Basis-term'):
            exp = basis_term.get('Exp')
            coeff = basis_term.get('Coeff')
            if float(exp) > self.cutoff:
                basis_part += ['  {} {}'.format(exp, coeff)]
                nterms+=1
        #now write the header 
        if nterms > 0:
          found_orbitals.append(angular)          
          charge=min(atom_charge-totcharge,maxcharge[angular])
          #put in a special case for transition metals:
          #depopulate the 4s if the atom is charged
          if symbol in transition_metals and symbol in self.initial_charges.keys() \
                  and self.initial_charges[symbol] > 0 and found_orbitals.count(angular) > 1 \
                  and angular=="s":
            charge=0
          totcharge+=charge
          ret+=["0 %i %i %g 1"%(basis_index[angular],nterms,charge)] + basis_part
          ncontract+=1

    #Add in the uncontracted basis elements
    angular_uncontracted=['s','p']
    if symbol in transition_metals:
        angular_uncontracted.append('d')

    for angular in angular_uncontracted:
        for i in range(0,self.basis_params[1]):
            exp=self.basis_params[0]*self.basis_params[2]**i
            line='{} {}'.format(exp,1.0)
            ret+=["0 %i %i %g 1"%(basis_index[angular],1,0.0),line]
            ncontract+=1

    return ["%i %i"%(Element(symbol).number+200,ncontract)] +\
            self.pseudopotential_section(symbol) +\
            ret
########################################################

  def pseudopotential_section(self,symbol):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com>
    Returns a string of the pseudopotential section which is to be written
    as part of the basis set section.

    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.
        xml_name (str): The name of the XML pseudopotential and basis
            set database.

    Returns:
        list of lines of pseudopotential section (edit by Brian Busemeyer).
    """
    tree = ElementTree()
    tree.parse(self.xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    eff_core_charge = element.find('./Effective_core_charge').text
    local_path = './Gaussian_expansion/Local_component'
    non_local_path = './Gaussian_expansion/Non-local_component'
    local_list = element.findall(local_path)
    non_local_list = element.findall(non_local_path)
    nlocal = len(local_list)
    m = [0, 0, 0, 0, 0]
    proj_path = './Gaussian_expansion/Non-local_component/Proj'
    proj_list = element.findall(proj_path)
    for projector in proj_list:
        m[int(projector.text)] += 1
    strlist = []
    strlist.append('INPUT')
    strlist.append(' '.join(map(str,[eff_core_charge,nlocal,
                                     m[0],m[1],m[2],m[3],m[4]])))
    for local_component in local_list:
        exp_gaus = local_component.find('./Exp').text
        coeff_gaus = local_component.find('./Coeff').text
        r_to_n = local_component.find('./r_to_n').text
        strlist.append(' '.join([exp_gaus, coeff_gaus,r_to_n]))
    for non_local_component in non_local_list:
        exp_gaus = non_local_component.find('./Exp').text
        coeff_gaus = non_local_component.find('./Coeff').text
        r_to_n = non_local_component.find('./r_to_n').text
        strlist.append(' '.join([exp_gaus, coeff_gaus,r_to_n]))
    return strlist
import os 


###################################################################

class CrystalReader:
  """ Extract properties of crystal run. 
  output values are stored in self.output dictionary when collect() is run. 
  """
  def __init__(self):
    self.completed=False
    self.output={}

    
#-------------------------------------------------      
  def collect(self,outfilename):
    """ Collect results from output."""
    # If the run didn't finish, then we won't find anything. 
    # In that case, we'll want to run again and collect again.
    status='killed'
    self.completed=False
    if os.path.isfile(outfilename):
      f = open(outfilename, 'r')
      try:
        lines = f.readlines()
      except UnicodeDecodeError:
        self.completed=False
        lines = []
        print(self.__class__.__name__,": Crystal output is unreadable, this usually happens when the process as been killed.")
      for li,line in enumerate(lines):
        if 'SCF ENDED - CONVERGENCE ON ENERGY' in line:
          self.output['total_energy']=float(line.split()[8])    
          print(self.__class__.__name__,": SCF ended converging on %f"%self.output['total_energy'])
          status='done'
          self.completed=True
        elif 'SCF ENDED - TOO MANY CYCLES' in line:
          last=float(line.split()[8])
          print("SCF ended at %f Ha without convergence"%last)
          status='killed'
          self.completed=False

        elif 'TOTAL ATOMIC SPINS' in line:
          moms = []
          shift = 1
          while "TTT" not in lines[li+shift]:
            moms += map(float,lines[li+shift].split())
            shift += 1
          self.output['mag_moments']=moms

        elif 'TOTAL ATOMIC CHARGES' in line:
          chgs = []
          shift = 1
          while ("SUMMED" not in lines[li+shift]) and ("TTT" not in lines[li+shift]):
            chgs += map(float,lines[li+shift].split())
            shift += 1
          self.output['atomic_charges']=chgs
    else:
      # Just to be sure/clear...
      self.completed=False
    return status


#-------------------------------------------------      
  def write_summary(self):
    print("Crystal total energy",self.output['total_energy'])


#-------------------------------------------------      
  # This can be made more efficient if it's a problem: searches whole file for
  # each query.
  def check_outputfile(outfilename,acceptable_scf=10.0):
    """ Check output file. 

    Return values:
    no_record, not_started, ok, too_many_cycles, finished (fall-back),
    scf_fail, not_enough_decrease, divergence, not_finished
    """
    if os.path.isfile(outfilename):
      outf = open(outfilename,'r',errors='ignore')
    else:
      return "not_started"

    outlines = outf.readlines()
    reslines = [line for line in outlines if "ENDED" in line]

    if len(reslines) > 0:
      if "CONVERGENCE" in reslines[0]:
        return "ok"
      elif "TOO MANY CYCLES" in reslines[0]:
        return "too_many_cycles"
      else: 
        return "finished"
      
    detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
    if len(detots) == 0:
      return "scf_fail"

    detots_net = sum(detots[1:])
    if detots_net > acceptable_scf:
      return "not_enough_decrease"

    etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
    if etots[-1] > 0:
      return "divergence"
    
    return "not_finished"
  
  
#-------------------------------------------------      
  def status(self,outfilename):
    """ Decide status of crystal run. """

    status=self.check_outputfile(outfilename)
    return status
    
if __name__=='__main__':
  print(space_group_format(136))
