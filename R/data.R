###################################################################/
# Descrption: sw.data.univariate
#
###################################################################/
sw.data.univariate = function( )
{
  file = system.file( 'ABO.1df.res.rdata', package = "schmidtWorkshop" )
  load( file )
  return( as.data.table( data[ !is.na( pValue )] ) )
}

###################################################################/
# Descrption: sw.data.tree
#
###################################################################/
sw.data.tree = function( )
{
  file = system.file( 'ABO.lk.surfs.rdata', package = "schmidtWorkshop" )
  load( file )
  tree = as.data.table( tree )
  tree[  , L_minus := d[ , 1 ] ]
  tree[  , L_0     := d[ , 2 ] ]
  tree[  , L_plus  := d[ , 3 ] ]

  return( tree )
}

###################################################################/
# Descrption: Data object with likelihood surfaces for the ABO
# SNP in the ICD-10 UK Biobank data set
###################################################################/
ABO.lk.surfs = function( ) {
    file = system.file( 'ABO.lk.surfs.rdata', package = "schmidtWorkshop" )
    load( file )
    return( ABO.data )
}

ABO.lk.pars = function( ) {
    file = system.file( 'ABO.lk.surfs.rdata', package = "schmidtWorkshop" )
    load( file )
    return( ABO.pars )
}


        


