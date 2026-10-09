netcdf gam3_sh.q{                                                               

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
                      .500000E+04,  .700000E+04,  .100000E+05,  .150000E+05,
                      .200000E+05,  .203471E+05,  .205516E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .983021E+05;
                q =   .360000E-05,  .341000E-05,  .325000E-05,  .308000E-05,
                      .268000E-05,  .237000E-05,  .201000E-05,  .264000E-05,
                      .783000E-05,  .932000E-05,  .102000E-04,  .271000E-04,
                      .751000E-04,  .278000E-03,  .612000E-03,  .160000E-02,
                      .321000E-02,  .584000E-02;

}                                                                               
