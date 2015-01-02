# TESTS
* refactor test suite

# Examples 
* improve examples; add more
* add one with an extended class meta bizness

# HackaMol::PeriodicTable needs to be reengineered 
* the data should be of higher quality
* access needs to be fast
* some datastructures such as KNOWN_NAMES accessible, and such mappings can be mutable.

# pdbqt files
* pull the tree of rotatable bonds groups
* pdbqt writing support 

# methods and attributes 
* add an attribute that allows a map from names to symbols/Zs to be built in reading pdb files etc

# DREAMS
* storage
* persistance
* MongoDB!

# clean-up classes
* remove the QmMolRol attributes from being loaded by default. keep role for extensions.  Make sure all current extensions pass tests. e.g. HackaMol::X::Vina depends on score attribute.  
* based on reading PDBQT may be able to simplify/rethink atom rotations about a bond.
* improve the binning in AtomGroupRole.  there is bin_this and then atom_bin 
stuff etc.
* improve use of hush_read attribute in MolReadRole.  Currently, only used for the dirty carp in PDB reading.

# Core Management
* figure out how to get Dist::Zilla to automate version numbers with some sanity
