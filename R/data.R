###################################################################/
# Descrption: mi.univariate.subset
#
###################################################################/
mi.univariate.subset = function( )
{
  selection = c( "I251", "E780", "M206", "S9211", "M4792", "I422", "C837", "I461", "B24", "G35")
  return( mi.univariate()[ Code %in% selection ][ order( pValue )] )
}

###################################################################/
# Descrption: mi.univariate
#
###################################################################/
mi.univariate = function( originalCols = FALSE )
{
  file = system.file( 'MI_GRS_tree.rdata', package = "schmidtWorkshop" )
  load( file )
  data = as.data.table( MI_GRS_tree )

  if( originalCols == TRUE )
    return( data )

  setnames( data, c( "coding", "meaning", "Naffected", "BETA", "Pval" ), c( "Code", "Description", "N", "beta", "pValue") )
  return( data[, .( Code, Description, N, beta, pValue ) ][ order( pValue ) ] )
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





