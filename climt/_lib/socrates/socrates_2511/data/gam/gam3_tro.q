netcdf gam3_tro.q{                                                              

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

    float q(plev,lon,lat);                                                     
             q:units = "None";                                                
             q:title = "MMR OF WATER VAPOUR";                                                

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .100000E+03,  .300000E+03,  .100000E+04,  .300000E+04,
                      .500000E+04,  .700000E+04,  .905919E+04,  .915024E+04,
                      .100000E+05,  .150000E+05,  .200000E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .989738E+05;
                q =   .345000E-05,  .307000E-05,  .287000E-05,  .275000E-05,
                      .251000E-05,  .232000E-05,  .237000E-05,  .238000E-05,
                      .239000E-05,  .660000E-05,  .220000E-04,  .618000E-04,
                      .281000E-03,  .883000E-03,  .178000E-02,  .475000E-02,
                      .904000E-02,  .152000E-01;

}                                                                               
