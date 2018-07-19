###################################################################/
# Descrption: plot.qq
#
# assume the data is normally distributed, convert to normal and
# plot against theoretical z-score
###################################################################/
plot.qq = function( data, annot.col = 'Description' )
{
    require(plotly)

    ## remove missing data
    data <- na.omit(data)

    ## convert the pValues to a normal standard error
    data[ , zData :=-qnorm( pValue )]

    ## now add theoritical z vales
    data = data[ order( zData ) ]
    data[ , zTheor := qnorm( (1:data[,.N ] - 0.5 )/data[,.N ] ) ]


    ## make plot
    plot = plot_ly( data, x = ~zTheor, y = ~zData, type = "scatter", mode = "markers", name = "data", text = data[, get( annot.col )], hoverinfo = "text" )
    plot = layout( plot,
      xaxis = list( title = "theoreitcal z-score", range = list( floor( data[,min( zTheor ) ] ), ceiling( data[, max( zTheor ) ] ) ) ),
      yaxis = list( title = "data z-score" ),
      title = "Q-Q Plot"
    )
    plot = add_lines( plot, x = ~zData, y = ~zData, name = "null" )

  return( plot )
}
