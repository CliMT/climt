netcdf lgam.t{                                                                  

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  16;


variables:
    float lat(lat);                                                           
             lat:units = "degree";                                            
             lat:title = "LATITUDE";                                         
    float lon(lon);                                                           
             lon:units = "degree";                                            
             lon:title = "LONGITUDE";                                        
    float plev(plev);                                                         
             plev:units = "Pa";                                               
             plev:title = "PRESSURE";                                        

    float t(plev,lon,lat);                                                     
             t:units = "K";                                                   
             t:title = "TEMPERATURE";                                                        

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .100000E+03,  .300000E+03,  .100000E+04,  .300000E+04,
                      .500000E+04,  .700000E+04,  .100000E+05,  .150000E+05,
                      .200000E+05,  .250000E+05,  .300000E+05,  .400000E+05,
                      .500000E+05,  .700000E+05,  .850000E+05,  .976700E+05;
                t =   .266290E+03,  .254310E+03,  .228110E+03,  .218710E+03,
                      .212590E+03,  .208640E+03,  .206890E+03,  .211830E+03,
                      .219010E+03,  .225870E+03,  .233270E+03,  .247190E+03,
                      .258620E+03,  .274000E+03,  .282280E+03,  .287130E+03;

}                                                                               
