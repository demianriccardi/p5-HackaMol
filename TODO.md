# Fix Kwalitee via CPAN

# Examples 
* improve examples; add more
* add one with an extended class meta bizness

# DREAMS
* storage
* persistance
* MongoDB!

#methods for analyzing wave functions (extensions)
*http://aim.tkgristmill.com/wfxformat.html

# clean-up classes
* remove the QmMolRol attributes from being loaded by default. keep role for extensions.  Make sure all current extensions pass tests. e.g. HackaMol::X::Vina depends on score attribute.  
* based on reading PDBQT may be able to simplify/rethink atom rotations about a bond.
* improve use of hush_read attribute in MolReadRole.  Currently, only used for the dirty carp in PDB reading.

