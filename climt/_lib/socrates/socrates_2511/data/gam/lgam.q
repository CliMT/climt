netcdf lgam.q{                                                                  

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

    float q(plev,lon,lat);                                                     
             q:units = "None";                                                
             q:title = "MMR OF WATER VAPOUR";                                                

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .100000E+03,  .300000E+03,  .100000E+04,  .300000E+04,
                      .500000E+04,  .700000E+04,  .100000E+05,  .150000E+05,
                      .200000E+05,  .250000E+05,  .300000E+05,  .400000E+05,
                      .500000E+05,  .700000E+05,  .850000E+05,  .976700E+05;
                q =   .350778E-05,  .324034E-05,  .307863E-05,  .292937E-05,
                      .263083E-05,  .240693E-05,  .228876E-05,  .475789E-05,
                      .152999E-04,  .456508E-04,  .185340E-03,  .605153E-03,
                      .123767E-02,  .333985E-02,  .640604E-02,  .104487E-01;

}                                                                               
