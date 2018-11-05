import sys
import numpy as np
sys.path.insert(0,'..')
import qwalk_objects as obj

CRYORB = 'crys.orb'

def run_test():
  crys_writer = test_crystal_writer()
  creader = test_crystal_reader()
  orbitals,system = convert_crystal()
  var = test_variance_writer()

def test_crystal_writer():
  cwriter = obj.crystal.CrystalWriter(xml_name='../BFD_Library.xml',total_spin=5)
  cwriter.set_struct_fromcif(open('mno/MnO.cif','r').read())
  cwriter.write_crys_input('mno/test/crys.in')
  cwriter.write_prop_input('mno/test/crys.prop.in')

def test_crystal_reader():
  creader = obj.crystal.CrystalReader()
  creader.collect('mno/ref/crys.in.o')
  return creader

def convert_crystal():
  system, orbitals = obj.crystal2qmc.pack_objects('mno/ref/crystal/GRED.DAT','mno/ref/crystal/KRED.DAT',spin=5)
  orbitals[0].write_qwalk_orb('mno/test/'+CRYORB)
  return orbitals,system

# NEXT STEP: write the tests for qwalk parts.
def test_variance_writer():
  system, orbitals = obj.crystal2qmc.pack_objects('mno/ref/crystal/GRED.DAT','mno/ref/crystal/KRED.DAT',spin=5) 
  slater = obj.trialfunc.Slater(orbitals[0],CRYORB,states=[[np.arange(system.nspin[0]),np.arange(system.nspin[1])]],shift_downorb=True)
  jastrow = system.export_jast()

  varwriter = obj.variance.VarianceWriter(system,obj.trialfunc.SlaterJastrow(slater,jastrow))
  varwriter.qwalk_input('mno/test/var')

if __name__=='__main__':
  run_test()
