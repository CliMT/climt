netcdf gam3_tro.t{                                                              

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  18;


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
                      .500000E+04,  .700000E+04,  .905919E+04,  .915024E+04,
                      .100000E+05,  .150000E+05,  .200000E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .989738E+05;
                t =   .267680E+03,  .256940E+03,  .230510E+03,  .218050E+03,
                      .208800E+03,  .202290E+03,  .199940E+03,  .199850E+03,
                      .199040E+03,  .206950E+03,  .219450E+03,  .230750E+03,
                      .240120E+03,  .255180E+03,  .266620E+03,  .281900E+03,
                      .290230E+03,  .297560E+03;

}                                                                               
